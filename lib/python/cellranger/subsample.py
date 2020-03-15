#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#
import cPickle
from collections import Counter
from itertools import izip

import numpy as np
import pandas as pd

import cellranger.constants as cr_constants
import cellranger.molecule_counter as cr_mc
import cellranger.rna.feature_ref as rna_feature_ref
import cellranger.rna.library as rna_library
import cellranger.utils as cr_utils
import tenkit.stats as tk_stats


def get_num_cells_per_library(library_info, filtered_barcodes_csv):
    """Get the number of cell-associated (i.e., filtered) barcodes per library.
    Note, we assume cell barcodes are assigned per GEM group, not per library.

    Args:
        library_info (dict): library_info metadata, probably from a MoleculeCounter
        filtered_barcodes_csv (str): path to filtered_barcodes.csv file
    Returns:
        np.array of int: number of cell-associated barcodes per library
    """
    # get number of cells per GEM group
    cell_bcs = cr_utils.get_cell_associated_barcode_set(filtered_barcodes_csv)
    num_cells_per_gg = Counter(cr_utils.split_barcode_seq(bc)[1] for bc in cell_bcs)
    # each library maps to a single gem group
    num_cells_per_lib = np.array(
        [num_cells_per_gg[lib['gem_group']] for lib in library_info])
    return num_cells_per_lib


def get_cell_associated_barcodes(genomes, filtered_barcodes_csv):
    """ Get cell-associated barcodes by genome.
    Args:
        genomes (list of str): Genome names.
        filtered_barcodes_csv (str): Path to CSV file.
    Returns:
        dict of (str, set): Map genome to list of cell-assoc barcodes.
            Empty-string key is for all genomes."""
    cell_bcs = {}
    for genome in genomes:
        # Get all cell-assoc barcodes (ignoring genome) for the '' (blank) genome string
        cell_bcs[genome] = cr_utils.get_cell_associated_barcode_set(
            filtered_barcodes_csv, genome)
    # All cell-associated barcodes
    cell_bcs[''] = set.union(*cell_bcs.itervalues())
    return cell_bcs


def compute_target_depths(max_target, num_targets):
    """ Construct a list of sorted, unique, integer-valued subsampling depths
    (generally corresponding to target read pairs per cell).

    Args:
        max_target (float): the largest target depth
        num_targets (int): desired number of targets, including max_target. There
            will be fewer than this many targets in case num_targets > max_target.

    Returns:
        numpy array of int: target subsampling depths (sorted, distinct, nonzero)"""
    targets_incl_zero = np.linspace(start=0, stop=max_target, num=num_targets + 1)
    distinct_targets = np.unique(np.floor(targets_incl_zero)).astype(int)
    return distinct_targets[distinct_targets > 0]


def make_subsamplings(subsample_type,
                      library_info, library_type,
                      num_cells_per_lib,
                      raw_reads_per_lib, usable_reads_per_lib,
                      fixed_depths, num_additional_depths):
    """Create metadata for subsampling jobs of a specified subsampling type
    (raw or usable reads per cell) and for a specified library type.

    Args:
        subsample_type (str): subsample based on raw or usable reads?
        library_info (dict): per-library metadata from MoleculeCounter
        library_type (str): library type to use for subsampling
        num_cells_per_lib (np.array of int): number of filtered barcodes per library
        raw_reads_per_lib (np.array of int): number of raw reads per library
        usable_reads_per_lib (np.array of int): number of usable reads per library
        fixed_depths (list of int): fixed subsampling depths (reads per cell)
            to include by default
        num_additional_depths (int): number of subsampling depths to use,
            in addition to the defaults

    Returns:
        list of dict: list of subsampling metadata, each of which is:
            {'library_type': <str>,
             'subsample_type': <str>,
             'target_read_pairs_per_cell': <int>,
             'library_subsample_rates': <np.array of floats>}
    """
    lib_indices = np.array([i for i, lib in enumerate(library_info)
                            if lib['library_type'] == library_type])

    raw_rppc_per_lib = raw_reads_per_lib.astype(float) / num_cells_per_lib
    usable_rppc_per_lib = usable_reads_per_lib.astype(float) / num_cells_per_lib
    rppc_per_lib = (raw_rppc_per_lib if subsample_type == cr_constants.RAW_SUBSAMPLE_TYPE
                    else usable_rppc_per_lib)
    # fraction of usable reads per library
    usable_frac_per_lib = usable_reads_per_lib.astype(float) / raw_reads_per_lib

    subsamplings = []
    # Pick a range of target depths that are feasible for all of the given libraries
    max_target_depth = np.min(rppc_per_lib[lib_indices])
    computed_depths = compute_target_depths(
        max_target_depth, num_additional_depths)
    # make note of the maximum after rounding
    # N.B. computed_depths is empty if max_target_depth is less than 1 read pair per cell.
    # This might happen if there are very few usable reads for this library type.
    max_computed_depth = (np.max(computed_depths)
                          if len(computed_depths) > 0 else None)
    target_depths = np.concatenate([computed_depths, fixed_depths])
    target_depths = np.unique(np.array(sorted(target_depths)))

    for i, target_depth in enumerate(target_depths):
        if subsample_type == cr_constants.MAPPED_SUBSAMPLE_TYPE:
            target_usable_reads_per_lib = target_depth * num_cells_per_lib
        else:
            # convert target raw read depth to usable read depth
            target_usable_reads_per_lib = (
                target_depth * num_cells_per_lib * usable_frac_per_lib)
        # compute subsampling rates (frac. of usable reads)
        subsample_rates = np.zeros(len(library_info), dtype=float)
        subsample_rates[lib_indices] = \
            (target_usable_reads_per_lib[lib_indices].astype(float)
             / usable_reads_per_lib[lib_indices])
        # for the largest computed (non-default) subsampling depth,
        # make sure we're subsampling the smallest library to rate=1.0
        if target_depth == max_computed_depth:
            subsample_rates = subsample_rates / np.max(subsample_rates)
        # zero out rates that are > 1
        # This can only apply to the the default subsampling targets,
        # for which we still want to run subsampling jobs with rate=0.
        subsample_rates[subsample_rates > 1.0] = 0.0

        subsamplings.append({
            'library_type': library_type,
            'subsample_type': subsample_type,
            'target_read_pairs_per_cell': target_depth,
            'library_subsample_rates': list(subsample_rates),
        })

    return subsamplings


def construct_all_subsamplings(molecule_info_h5,
                               filtered_barcodes_csv,
                               is_targeted=False):
    """Construct subsampling metadata for a range of target read depths,
    both raw reads per cell and usable reads per cell.

    Args:
        molecule_info_h5 (str): path to molecule_info.h5 file
        filtered_barcodes_csv (str): path to filtered_barcodes.csv file
        is_targeted (bool, optional): when subsampling to usable reads per cell,
            also restrict to on-target reads. Defaults to False.

    Returns:
        list of dict: metadata for subsampling job, produced by make_subsamplings
            and consumed by run_subsampling
    """
    # Get required info from the mol info
    with cr_mc.MoleculeCounter.open(molecule_info_h5, 'r') as mc:
        library_info = mc.library_info
        num_cells_per_lib = get_num_cells_per_library(
            library_info, filtered_barcodes_csv)

        raw_reads_per_lib = np.array(mc.get_raw_read_pairs_per_library())
        if not is_targeted:
            usable_reads_per_lib = np.array(mc.get_usable_read_pairs_per_library())
        else:
            usable_reads_per_lib = np.array(mc.get_on_target_usable_read_pairs_per_library())

    if num_cells_per_lib.sum() == 0:
        return []

    subsamplings = []

    for library_type in rna_library.sorted_library_types(library_info):
        libraries = [l for l in library_info
                     if l['library_type'] == library_type]
        is_targeted_gex = (library_type == rna_library.GENE_EXPRESSION_LIBRARY_TYPE
                           and any(rna_library.has_target_set(lib) for lib in libraries))
        if is_targeted_gex:
            fixed_depths = cr_constants.SUBSAMPLE_TARGETED_FIXED_DEPTHS
        else:
            fixed_depths = cr_constants.SUBSAMPLE_FIXED_DEPTHS

        for subsample_type in cr_constants.ALL_SUBSAMPLE_TYPES:
            subsamplings.extend(make_subsamplings(
                subsample_type,
                library_info,
                library_type,
                num_cells_per_lib,
                raw_reads_per_lib,
                usable_reads_per_lib,
                fixed_depths,
                cr_constants.SUBSAMPLE_NUM_ADDITIONAL_DEPTHS))

    return subsamplings


def run_subsampling(molecule_info_h5, subsample_info, filtered_barcodes_csv, feature_indices,
                    chunk_start, chunk_len):
    """
    Runs a subsampling chunk.
    :param molecule_info_h5: Path to a MoleculeCounter file.
    :param subsample_info: A subsampling produced by make_subsamplings
    :param filtered_barcodes_csv: A CSV of filtered (cell) barcodes
    :param feature indices: indices of filtered features
    :param chunk_start: integer chunk start
    :param chunk_len: integer chunk len
    :return: data dictionary
    """
    np.random.seed(0)
    mc = cr_mc.MoleculeCounter.open(molecule_info_h5, 'r')
    # Get cell-associated barcodes
    genomes = mc.feature_reference.get_genomes()
    cell_bcs_by_genome = get_cell_associated_barcodes(genomes, filtered_barcodes_csv)

    # Load chunk of relevant data from the mol_info
    chunk = slice(int(chunk_start), int(chunk_start) + int(chunk_len))
    mol_library_idx = mc.get_column_lazy('library_idx')[chunk]
    mol_read_pairs = mc.get_column_lazy('count')[chunk]
    mol_gem_group = mc.get_column_lazy('gem_group')[chunk]
    mol_barcode_idx = mc.get_column_lazy('barcode_idx')[chunk]
    mol_feature_idx = mc.get_column_lazy('feature_idx')[chunk]

    barcodes = mc.get_ref_column('barcodes')

    if feature_indices is not None:
        # subset molecules to targeted panel for each library
        feature_indices = np.array(sorted(feature_indices))
        mask = np.isin(mol_feature_idx, feature_indices)
        mol_library_idx = mol_library_idx[:][mask]
        mol_read_pairs = mol_read_pairs[:][mask]
        mol_gem_group = mol_gem_group[:][mask]
        mol_barcode_idx = mol_barcode_idx[:][mask]
        mol_feature_idx = mol_feature_idx[:][mask]

    # Give each cell-associated barcode an integer index
    cell_bcs = sorted(list(cell_bcs_by_genome['']))
    cell_bc_to_int = {bc: i for i, bc in enumerate(cell_bcs)}

    # Give each genome an integer index
    genome_to_int = {g: i for i, g in enumerate(genomes)}
    feature_int_to_genome_int = np.fromiter(
        (genome_to_int[f.tags.get(rna_feature_ref.GENOME_FEATURE_TAG, '')] for f in mc.feature_reference.feature_defs),
        dtype=int)
    mol_genome_idx = feature_int_to_genome_int[mol_feature_idx]

    # determine which (library type, genome) pairs have any associated reads
    lib_types = rna_library.sorted_library_types(mc.library_info)
    lib_type_to_int = {l: i for i, l in enumerate(lib_types)}
    lib_idx_to_lib_type_idx = np.fromiter((lib_type_to_int[lib['library_type']] for lib in mc.library_info),
                                          dtype=np.int)

    lib_type_genome_any_reads = np.zeros((len(lib_types), len(genomes)), dtype=np.bool)
    lib_genome_idx_pairs = set(izip(mol_library_idx[mol_read_pairs > 0],
                                    mol_genome_idx[mol_read_pairs > 0]))
    for (lib_idx, genome_idx) in lib_genome_idx_pairs:
        lib_type_idx = lib_idx_to_lib_type_idx[lib_idx]
        lib_type_genome_any_reads[lib_type_idx, genome_idx] = True

    # Run each subsampling task on this chunk of data
    n_tasks = len(subsample_info)
    n_genomes = len(genomes)
    n_cells = len(cell_bcs)

    umis_per_bc = np.zeros((n_tasks, n_genomes, n_cells))
    read_pairs_per_bc = np.zeros((n_tasks, n_genomes, n_cells))
    features_det_per_bc = np.zeros((n_tasks, n_genomes, n_cells))
    read_pairs_per_task = np.zeros((n_tasks, n_genomes))
    umis_per_task = np.zeros((n_tasks, n_genomes))

    for task_idx, task in enumerate(subsample_info):
        print 'subsampling task: {}'.format(task)
        # Per-library subsampling rates
        rates_per_library = np.array(task['library_subsample_rates'], dtype=float)

        if np.count_nonzero(rates_per_library) == 0:
            continue

        mol_rate = rates_per_library[mol_library_idx]

        # Subsampled read pairs per molecule
        new_read_pairs = np.random.binomial(mol_read_pairs, mol_rate)

        # Compute tallies for each barcode
        group_keys = (mol_gem_group, mol_barcode_idx)
        group_values = (mol_feature_idx, mol_genome_idx, new_read_pairs)
        for (gg, bc_idx), (feature_idx, genome_idx, read_pairs) in \
                cr_utils.numpy_groupby(group_values, group_keys):

            barcode = cr_utils.format_barcode_seq(barcodes[bc_idx], gg)

            cell_idx = cell_bc_to_int.get(barcode)

            for this_genome_idx in xrange(len(genomes)):
                umis = np.flatnonzero((read_pairs > 0) & (genome_idx == this_genome_idx))
                this_genome_read_pairs = np.sum(read_pairs[genome_idx == this_genome_idx])

                # Tally UMIs and median features detected
                if barcode in cell_bcs_by_genome[genomes[this_genome_idx]]:
                    # This is a cell-associated barcode for this genome
                    umis_per_bc[task_idx, this_genome_idx, cell_idx] = len(umis)
                    read_pairs_per_bc[task_idx, this_genome_idx, cell_idx] = this_genome_read_pairs
                    features_det_per_bc[task_idx, this_genome_idx, cell_idx] = np.count_nonzero(
                        np.bincount(feature_idx[umis]))

                # Tally numbers for duplicate fraction
                read_pairs_per_task[task_idx, this_genome_idx] += this_genome_read_pairs
                umis_per_task[task_idx, this_genome_idx] += len(umis)

    data = {
        'umis_per_bc': umis_per_bc,
        'features_det_per_bc': features_det_per_bc,
        'read_pairs_per_bc': read_pairs_per_bc,
        'read_pairs': read_pairs_per_task,
        'umis': umis_per_task,
        'lib_type_genome_any_reads': lib_type_genome_any_reads,
    }
    mc.close()
    return data


def make_metric_name(name, library_type, genome, ss_type, ss_depth, target_mode):
    lt_prefix = rna_library.get_library_type_metric_prefix(library_type)
    if target_mode is None:
        return '%s%s_%s_%s_%s' % (lt_prefix, genome, ss_type, ss_depth, name)
    return '%s%s_%s_%s_%s_%s' % (lt_prefix, genome, ss_type, ss_depth, name, target_mode)


def compute_dup_frac(read_pairs, umis):
    return tk_stats.robust_divide(read_pairs - umis, read_pairs) if read_pairs > 0 else 0.0


def join_metrics(metrics):
    """
    Join together a list of metric dicts (each value is a numpy vector)
    :param metrics: list of metrics
    :return: joined dictionary
    """
    if len(metrics) == 0:
        return dict()
    with open(metrics[0]) as f:
        data = cPickle.load(f)
    for m in metrics[1:]:
        with open(m) as f:
            chunk_data = cPickle.load(f)
            for k, v in chunk_data.iteritems():
                data[k] += v
    return data


def calculate_subsampling_metrics(data,
                                  molecule_info_h5,
                                  filtered_barcodes_csv,
                                  subsample_info,
                                  target_mode):
    """
    Calculate subsampling metrics (summary) from a joined data structure from run_subsampling
    :param data: A merged dictionary of data from run_subsampling
    :param molecule_info_h5: path to a MoleculeInfo file
    :param filtered_barcodes_csv: path to a list of cell barcodes
    :param subsample_info: subsampling info produced by construct_all_subsamplings
    :param target_mode: String of target mode for metrics suffix. One of ['ontarget', 'offtarget', None].
    :return: dict (JSON) metrics
    """
    # Compute metrics for each subsampling rate
    summary = {}

    with cr_mc.MoleculeCounter.open(molecule_info_h5, 'r') as mc:
        genomes = mc.feature_reference.get_genomes()
        lib_types = rna_library.sorted_library_types(mc.library_info)
        lib_type_map = dict((lt, idx) for (idx, lt) in enumerate(lib_types))
    cell_bcs_by_genome = get_cell_associated_barcodes(genomes, filtered_barcodes_csv)

    # Give each cell-associated barcode an integer index
    cell_bcs = sorted(list(cell_bcs_by_genome['']))
    cell_bc_to_int = {bc: i for i, bc in enumerate(cell_bcs)}

    for i, task in enumerate(subsample_info):
        lib_type = task['library_type']
        lib_type_idx = lib_type_map[lib_type]
        ss_type = task['subsample_type']
        ss_depth = task['target_read_pairs_per_cell']

        if rna_library.has_genomes(lib_type):
            genome_ints = list(range(data['umis_per_bc'].shape[1]))
        else:
            genome_ints = [0]

        # Per-genome metrics
        for g in genome_ints:
            if not data['lib_type_genome_any_reads'][lib_type_idx, g]:
                continue
            genome = genomes[g]

            # Only compute on cell-associated barcodes for this genome.
            # This only matters when there are multiple genomes present.
            cell_inds = np.array(sorted(cell_bc_to_int[bc] for bc in cell_bcs_by_genome[genome]))

            mean_reads_per_cell = np.mean(data['read_pairs_per_bc'][i,g,cell_inds])
            summary[make_metric_name('subsampled_filtered_bcs_mean_read_counts',
                                     lib_type, genome, ss_type, ss_depth, target_mode)] = mean_reads_per_cell

            median_reads_per_cell = np.median(data['read_pairs_per_bc'][i,g,cell_inds])
            summary[make_metric_name('subsampled_filtered_bcs_median_read_counts',
                                     lib_type, genome, ss_type, ss_depth, target_mode)] = median_reads_per_cell

            median_umis_per_cell = np.median(data['umis_per_bc'][i,g,cell_inds])
            summary[make_metric_name('subsampled_filtered_bcs_median_counts',
                                     lib_type, genome, ss_type, ss_depth, target_mode)] = median_umis_per_cell

            median_features_per_cell = np.median(data['features_det_per_bc'][i,g,cell_inds])
            summary[make_metric_name('subsampled_filtered_bcs_median_unique_genes_detected',
                                     lib_type, genome, ss_type, ss_depth, target_mode)] = median_features_per_cell

            dup_frac = compute_dup_frac(data['read_pairs'][i,g],  data['umis'][i,g])
            summary[make_metric_name('subsampled_duplication_frac',
                                     lib_type, genome, ss_type, ss_depth, target_mode)] = dup_frac

        # Whole-dataset duplication frac
        all_read_pairs = np.sum(data['read_pairs'][i,:])
        all_umis = np.sum(data['umis'][i,:])
        dup_frac = compute_dup_frac(all_read_pairs, all_umis)
        summary[make_metric_name('subsampled_duplication_frac',
                                 lib_type, rna_library.MULTI_REFS_PREFIX, ss_type, ss_depth, target_mode)] = dup_frac
    return summary


def parse_subsample_summary(summary, suffix=None, prefix=None):
    """
    Convenience function for extracting a subsampling results table from a summary JSON.
    :param summary: a loaded JSON, as produced by calculate_subsampling_metrics
    :return: pandas DataFrame
    """
    munged_data = []
    names = ['Mean reads per cell',
             'Median reads per cell',
             'Median UMIs per cell',
             'Median genes per cell']
    metrics = ['subsampled_filtered_bcs_mean_read_counts',
               'subsampled_filtered_bcs_median_read_counts',
               'subsampled_filtered_bcs_median_counts',
               'subsampled_filtered_bcs_median_unique_genes_detected']
    if suffix is not None:
        metrics = ['{}_{}'.format(x, suffix) for x in metrics]
    for name, metric in zip(names, metrics):
        keys = [x for x in summary.keys() if metric in x and cr_constants.MAPPED_SUBSAMPLE_TYPE in x]
        if prefix is not None:
            keys = [x for x in keys if x.startswith(prefix)]
        ss_target = [int(x.split(metric)[0].rsplit('_', 2)[1]) for x in keys]
        data = [summary[x] for x in keys]
        for target, val in zip(ss_target, data):
            munged_data.append([target, val, name])
    df = pd.DataFrame(munged_data, columns=['target', 'val', 'name'])
    df = df.pivot(index='target', columns='name', values='val').astype(int)
    # remove rows that are all 0s
    df = df[(df != 0).any(axis=1)]
    # handle the case where all subsamplings were skipped due to super low coverage
    if df.shape == (0, 0):
        df = pd.DataFrame(columns=names)
    df.index.name = 'Target mean reads per cell'
    return df
