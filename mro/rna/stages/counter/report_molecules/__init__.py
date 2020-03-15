#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
""" Report info on each detected molecule
"""

from collections import defaultdict, OrderedDict, Counter
import itertools
import math
import numpy as np
import tenkit.bam as tk_bam
import tenkit.stats as tk_stats
from tenkit.safe_json import safe_jsonify
import cellranger.constants as cr_constants
import cellranger.molecule_counter as cr_mol_counter
from cellranger.molecule_counter import MoleculeCounter
import cellranger.report as cr_report
import cellranger.io as cr_io
import cellranger.rna.feature_ref as rna_feature_ref
import cellranger.rna.library as rna_library
import cellranger.utils as cr_utils

__MRO__ = '''
stage REPORT_MOLECULES(
    in  bam[]  inputs,
    in  path   reference_path,
    in  csv[]  target_features,
    in  csv    feature_reference,
    in  map    align,
    in  string barcode_whitelist,
    in  json   extract_reads_summary,
    in  json   attach_bcs_and_umis_summary,
    in  json   mark_duplicates_summary,
    in  csv    filtered_barcodes,
    in  int    recovered_cells,
    in  int    force_cells,
    in  bool   disable_targeted,
    out h5     output,
    out json   summary,
    src py     "../rna/stages/counter/report_molecules",
) split (
    in  bam    chunk_input,
) using (
    mem_gb   = 2,
    volatile = strict,
)
'''

MAX_MEM_GB = 64

def split(args):
    # Estimate the total number of rows in the final molecule info. Worst case.
    if cr_utils.is_metric_in_json(args.extract_reads_summary, 'total_reads'):
        total_reads = cr_utils.get_metric_from_json(args.extract_reads_summary, 'total_reads')
    else: # there is no GEX library, query the antibody library total reads instead
        total_reads = cr_utils.get_metric_from_json(args.extract_reads_summary, rna_library.get_library_type_metric_prefix(rna_library.ANTIBODY_METRIC_PREFIX) + 'total_reads')
    mol_info_rows = total_reads

    # Memory for chunk
    if len(args.inputs) > 0:
        avg_rows_per_chunk = int(total_reads / len(args.inputs))
        avg_chunk_mem_gb = int(math.ceil((32 * avg_rows_per_chunk)/2.5e8))
        chunk_mem_gb = min(MAX_MEM_GB, max(8, avg_chunk_mem_gb))
    else:
        chunk_mem_gb = 1

    # Memory for concatenating molecule info
    # N = total number of rows
    # 8*N bytes to store the sort indices
    # (8+8+8)*N bytes to load, concatenate, and index into a 64-bit data column
    mol_info_mem_gb = int(math.ceil((32 * mol_info_rows)/2.5e8))
    join_mem_gb = min(MAX_MEM_GB, max(4, mol_info_mem_gb))

    chunks = []
    for chunk_input in args.inputs:
        chunks.append({
            'chunk_input': chunk_input,
            '__mem_gb': chunk_mem_gb,
        })
    join_def = {
        '__mem_gb': join_mem_gb,
    }
    return {'chunks': chunks, 'join': join_def}

def main(args, outs):
    outs.coerce_strings()

    # Load whitelist
    whitelist = cr_utils.load_barcode_whitelist(args.barcode_whitelist)
    barcode_to_idx = OrderedDict((k, i) for i,k in enumerate(whitelist))

    # Load library info from BAM
    in_bam = tk_bam.create_bam_infile(args.chunk_input)
    library_info = rna_library.get_bam_library_info(in_bam)

    if not args.disable_targeted and args.target_features is not None:
        # Load target features
        target_features = cr_io.load_target_feature_dict(library_info, args.target_features)
        # load library_idx order
        by_library_features = [set(cr_io.load_target_features_csv(x)) for x in args.target_features]
    else:
        target_features = None
    # Load feature reference
    feature_ref = rna_feature_ref.from_transcriptome_and_csv(args.reference_path,
                                                             args.feature_reference,
                                                             target_features)

    # Get cell-associated barcodes by genome
    filtered_bcs_by_genome = cr_utils.load_barcode_csv(args.filtered_barcodes)
    filtered_bc_union = cr_utils.get_cell_associated_barcode_set(args.filtered_barcodes)

    # Create the barcode info
    barcode_info = MoleculeCounter.build_barcode_info(filtered_bcs_by_genome,
                                                      library_info,
                                                      whitelist)

    # Create the molecule info file
    mc = MoleculeCounter.open(outs.output, mode='w',
                              feature_ref=feature_ref,
                              barcodes=whitelist,
                              library_info=library_info,
                              barcode_info=barcode_info)

    # Initialize per-library metrics
    lib_metrics = {}
    for lib_idx in xrange(len(library_info)):
        lib_metrics[str(lib_idx)] = {}
        lib_metrics[str(lib_idx)][cr_mol_counter.USABLE_READS_METRIC] = 0
        if not args.disable_targeted:
            lib_metrics[str(lib_idx)][cr_mol_counter.ON_TARGET_USABLE_READS_METRIC] = 0

    # Record read-counts per molecule. Note that UMIs are not contiguous
    # in the input because no sorting was done after UMI correction.

    prev_gem_group = None
    prev_barcode_idx = None

    for (gem_group, barcode_seq), reads_iter in \
        itertools.groupby(in_bam, key=cr_utils.barcode_sort_key_no_umi):
        if barcode_seq is None:
            continue

        barcode_idx = barcode_to_idx[barcode_seq]

        # Assert expected sort order of input BAM
        assert gem_group >= prev_gem_group
        if gem_group == prev_gem_group:
            assert barcode_idx >= prev_barcode_idx

        is_cell_barcode = cr_utils.format_barcode_seq(barcode_seq, gem_group) in filtered_bc_union

        counts = defaultdict(int)

        for read in reads_iter:
            # ignore read2 to avoid double-counting. the mapping + annotation should be equivalent.
            if read.is_secondary or \
               read.is_read2 or \
               cr_utils.is_read_low_support_umi(read) or \
               not cr_utils.is_read_conf_mapped_to_feature(read):
                continue

            umi_seq = cr_utils.get_read_umi(read)
            if umi_seq is None:
                continue

            umi_int = MoleculeCounter.compress_umi_seq(umi_seq, MoleculeCounter.get_column_dtype('umi').itemsize*8)

            feature_ids = cr_utils.get_read_gene_ids(read)
            assert len(feature_ids) == 1
            feature_int = feature_ref.id_map[feature_ids[0]].index

            library_idx = cr_utils.get_read_library_index(read)

            counts[(umi_int, library_idx, feature_int)] += 1

            if is_cell_barcode:
                lib_metrics[str(library_idx)][cr_mol_counter.USABLE_READS_METRIC] += 1
                if not args.disable_targeted and feature_int in by_library_features[library_idx]:
                    lib_metrics[str(library_idx)][cr_mol_counter.ON_TARGET_USABLE_READS_METRIC] += 1

            prev_gem_group = gem_group
            prev_barcode_idx = barcode_idx

        # Record data for this barcode
        gg_int = MoleculeCounter.get_column_dtype('gem_group').type(gem_group)
        mc.append_column('gem_group', np.repeat(gg_int, len(counts)))
        bc_int = MoleculeCounter.get_column_dtype('barcode_idx').type(barcode_idx)
        mc.append_column('barcode_idx', np.repeat(bc_int, len(counts)))

        feature_ints = np.fromiter((k[2] for k in counts.iterkeys()),
                                   dtype=MoleculeCounter.get_column_dtype('feature_idx'),
                                   count=len(counts))
        # Sort by feature for fast matrix construction
        order = np.argsort(feature_ints)
        feature_ints = feature_ints[order]
        mc.append_column('feature_idx', feature_ints)
        del feature_ints

        li_ints = np.fromiter((k[1] for k in counts.iterkeys()),
                               dtype=MoleculeCounter.get_column_dtype('library_idx'),
                               count=len(counts))[order]
        mc.append_column('library_idx', li_ints)
        del li_ints

        umi_ints = np.fromiter((k[0] for k in counts.iterkeys()),
                               dtype=MoleculeCounter.get_column_dtype('umi'),
                               count=len(counts))[order]
        mc.append_column('umi', umi_ints)
        del umi_ints

        count_ints = np.fromiter(counts.itervalues(),
                                 dtype=MoleculeCounter.get_column_dtype('count'),
                                 count=len(counts))[order]
        mc.append_column('count', count_ints)
        del count_ints

    in_bam.close()

    mc.set_metric(cr_mol_counter.LIBRARIES_METRIC, dict(lib_metrics))

    mc.save()

def join(args, outs, chunk_defs, chunk_outs):
    summary = cr_utils.merge_jsons_as_dict([
        args.extract_reads_summary,
        args.attach_bcs_and_umis_summary,
        args.mark_duplicates_summary,
    ])

    # Hack for getting reference metadata -
    # this used to be computed in prior stages.
    # This is needed for storage in the molecule_info HDF5.
    tmp_reporter = cr_report.Reporter()
    tmp_reporter.store_reference_metadata(args.reference_path,
                                          cr_constants.REFERENCE_TYPE,
                                          cr_constants.REFERENCE_METRIC_PREFIX)
    ref_metadata = tmp_reporter.report(cr_constants.DEFAULT_REPORT_TYPE)
    summary.update(ref_metadata)

    # Load library info from BAM
    in_bam = tk_bam.create_bam_infile(args.inputs[0])
    library_info = rna_library.get_bam_library_info(in_bam)

    metrics = MoleculeCounter.get_metrics_from_summary(summary, library_info,
                                                       args.recovered_cells, args.force_cells)

    input_h5_filenames = [chunk_out.output for chunk_out in chunk_outs]
    # update with metrics that were computed in the chunks
    chunk_metrics = [cr_mol_counter.USABLE_READS_METRIC]
    if not args.disable_targeted:
        chunk_metrics.append(cr_mol_counter.ON_TARGET_USABLE_READS_METRIC)
    for chunk_metric in chunk_metrics:
        summed_lib_metrics = MoleculeCounter.sum_library_metric(input_h5_filenames, chunk_metric)
        for lib_key, value in summed_lib_metrics.iteritems():
            metrics[cr_mol_counter.LIBRARIES_METRIC][lib_key][chunk_metric] = value

    MoleculeCounter.concatenate(outs.output, input_h5_filenames, metrics=metrics)

    # write out targeting specific metrics
    with MoleculeCounter.open(outs.output, 'r') as mc:
        genomes = mc.feature_reference.get_genomes()
        lib_types = rna_library.sorted_library_types(mc.library_info)

    if not args.disable_targeted:
        # group the on-target metric by genome
        json_metrics = Counter()
        # keep track of totals for fractions
        read_counter = Counter()
        for lib_idx, value in metrics[cr_mol_counter.LIBRARIES_METRIC].iteritems():
            genome = genomes[int(lib_idx)]
            lib_type = lib_types[int(lib_idx)]
            genome_key = make_metric_name('on_target_conf_mapped_reads', lib_type, genome)
            multi_key = make_metric_name('on_target_conf_mapped_reads', lib_type, 'multi')
            for key in [genome_key, multi_key]:
                json_metrics[key] += value['on_target_usable_read_pairs']
                read_counter[key] += value['raw_read_pairs']

        # calculate per-genome and combined fractions
        for key, val in read_counter.iteritems():
            json_metrics[key + '_frac'] = tk_stats.robust_divide(json_metrics[key], val)

        with open(outs.summary, 'w') as fh:
            fh.write(safe_jsonify(json_metrics))


def make_metric_name(name, library_type, genome):
    lt_prefix = rna_library.get_library_type_metric_prefix(library_type)
    return '%s%s_%s' % (lt_prefix, genome, name)
