import enum
import json
import numpy as np
import pandas as pd
import martian
from collections import OrderedDict


import cellranger.rna.library as rna_library
import cellranger.cell_calling as cr_cell
import cellranger.feature.antibody.analysis as ab_utils

import cellranger.constants as cr_constants

from cellranger.metrics import BarcodeFilterResults
from cellranger.stats import determine_max_filtered_bcs, find_within_ordmag, summarize_bootstrapped_top_n
import tenkit.stats as tk_stats

class FilterMethod(enum.Enum):
    # Caller-provided list of barcodes
    MANUAL = 0
    # Take the top N barcodes by count
    TOP_N_BARCODES = 1
    # Take barcodes within an order of magnitude of the max by count
    ORDMAG = 2
    # The above, then find barcodes that differ from the ambient profile
    ORDMAG_NONAMBIENT = 3

def get_filter_method_name(fm):
    if fm == FilterMethod.MANUAL:
        return 'manual'
    elif fm == FilterMethod.TOP_N_BARCODES:
        return 'topn'
    elif fm == FilterMethod.ORDMAG:
        return 'ordmag'
    elif fm == FilterMethod.ORDMAG_NONAMBIENT:
        return 'ordmag_nonambient'
    else:
        raise ValueError('Unsupported filter method value %d' % fm)


###############################################################################
def remove_bcs_with_high_umi_corrected_reads(correction_data, matrix):
    """ Given a CountMatrix and and csv file containing information about umi corrected reads,
        detect all barcodes with unusually high fraction of corrected reads (probably aggregates),
        and remove them from the CoutMatrix """

    bcs_to_remove, reads_lost, removed_bcs_df = ab_utils.detect_aggregate_bcs(correction_data)
    bcs_to_remove = set(matrix.bc_to_int(bc) for bc in bcs_to_remove)
    # make sure filtered_bcs is in deterministic order or any later bootstrap sampling will not be deterministic
    filtered_bcs = [i for i in xrange(matrix.bcs_dim) if i not in bcs_to_remove]
    cleaned_matrix = matrix.select_barcodes(filtered_bcs)

    ### report how many aggregates were found, and the fraction of reads those accounted for
    metrics_to_report = {}
    report_prefix  = rna_library.get_library_type_metric_prefix(rna_library.ANTIBODY_LIBRARY_TYPE)
    metrics_to_report[report_prefix + 'number_highly_corrected_GEMs'] = len(bcs_to_remove)
    metrics_to_report[report_prefix + 'reads_lost_to_highly_corrected_GEMs'] = reads_lost

    return cleaned_matrix, metrics_to_report, removed_bcs_df

################################################################################
def call_initial_cells(matrix, genomes, is_antibody_only, unique_gem_groups, method, gg_recovered_cells, cell_barcodes, force_cells):
    # (gem_group, genome) => dict
    filtered_metrics_groups = OrderedDict()
    # (gem_group, genome) => list of barcode strings
    filtered_bcs_groups = OrderedDict()

    for genome in genomes:
        if is_antibody_only: # this is the case when there in only Antibody Capture library present
            genome_matrix = matrix
        else:
            genome_matrix = matrix.select_features_by_genome(genome)

        # Make initial cell calls for each gem group individually
        for gem_group in unique_gem_groups:
            gg_matrix = genome_matrix.select_barcodes_by_gem_group(gem_group)
            gg_filtered_metrics, gg_filtered_bcs = _call_cells_by_gem_group(gg_matrix, unique_gem_groups, method, gg_recovered_cells, cell_barcodes, force_cells)
            filtered_metrics_groups[(gem_group, genome)] = gg_filtered_metrics
            filtered_bcs_groups[(gem_group, genome)] = gg_filtered_bcs

    return filtered_metrics_groups, filtered_bcs_groups


def filter_cellular_barcodes_manual(matrix, cell_barcodes):
    """ Take all barcodes that were given as cell barcodes """
    subset = matrix.select_barcodes_by_seq(cell_barcodes)
    bc_idx = np.where(subset.bc_mask)[0]
    metrics = BarcodeFilterResults.init_with_constant_call(len(bc_idx))
    return bc_idx, metrics, None



def _call_cells_by_gem_group(gg_matrix, unique_gem_groups, method, gg_recovered_cells, cell_barcodes, force_cells):

    if method == FilterMethod.ORDMAG or method == FilterMethod.ORDMAG_NONAMBIENT:
        gg_total_diversity = gg_matrix.bcs_dim
        gg_bc_counts = gg_matrix.get_counts_per_bc()
        gg_filtered_indices, gg_filtered_metrics, msg = filter_cellular_barcodes_ordmag(
            gg_bc_counts, gg_recovered_cells, gg_total_diversity)
        gg_filtered_bcs = gg_matrix.ints_to_bcs(gg_filtered_indices)

    elif method == FilterMethod.MANUAL:
        with(open(cell_barcodes)) as f:
            cell_barcodes = json.load(f)
        gg_filtered_indices, gg_filtered_metrics, msg = filter_cellular_barcodes_manual(
            gg_matrix, cell_barcodes)
        gg_filtered_bcs = gg_matrix.ints_to_bcs(gg_filtered_indices)

    elif method == FilterMethod.TOP_N_BARCODES:
        gg_bc_counts = gg_matrix.get_counts_per_bc()
        if force_cells is not None:
            gg_force_cells = int(float(force_cells) / float(len(unique_gem_groups)))
        gg_filtered_indices, gg_filtered_metrics, msg = filter_cellular_barcodes_fixed_cutoff(
            gg_bc_counts, gg_force_cells)
        gg_filtered_bcs = gg_matrix.ints_to_bcs(gg_filtered_indices)

    else:
        martian.exit("Unsupported BC filtering method: %s" % method)

    if msg is not None:
        martian.log_info(msg)

    return gg_filtered_metrics, gg_filtered_bcs

def call_additional_cells(matrix, unique_gem_groups, genomes, filtered_bcs_groups):
    # Track these for recordkeeping
    eval_bcs_arrays = []
    umis_per_bc_arrays = []
    loglk_arrays = []
    pvalue_arrays = []
    pvalue_adj_arrays = []
    nonambient_arrays = []
    genome_call_arrays = []

    # Do it by gem group, but agnostic to genome
    for gg in unique_gem_groups:
        gg_matrix = matrix.select_barcodes_by_gem_group(gg)

        # Take union of initial cell calls across genomes
        gg_bcs = sorted(list(reduce(set.union,
                                    [set(bcs) for group, bcs in filtered_bcs_groups.iteritems() if group[0] == gg])))

        result = cr_cell.find_nonambient_barcodes(gg_matrix, gg_bcs)
        if result is None:
            print 'Failed at attempt to call non-ambient barcodes in GEM well %s' % gg
            continue

        # Assign a genome to the cell calls by argmax genome counts
        genome_counts = []
        for genome in genomes:
            genome_counts.append(gg_matrix.view() \
                                    .select_features_by_genome(genome) \
                                    .select_barcodes(result.eval_bcs) \
                                    .get_counts_per_bc())
        genome_counts = np.column_stack(genome_counts)
        genome_calls = np.array(genomes)[np.argmax(genome_counts, axis=1)]

        umis_per_bc = gg_matrix.get_counts_per_bc()

        eval_bcs_arrays.append(np.array(gg_matrix.bcs)[result.eval_bcs])
        umis_per_bc_arrays.append(umis_per_bc[result.eval_bcs])
        loglk_arrays.append(result.log_likelihood)
        pvalue_arrays.append(result.pvalues)
        pvalue_adj_arrays.append(result.pvalues_adj)
        nonambient_arrays.append(result.is_nonambient)
        genome_call_arrays.append(genome_calls)

        # Update the lists of cell-associated barcodes
        for genome in genomes:
            eval_bc_strs = np.array(gg_matrix.bcs)[result.eval_bcs]
            filtered_bcs_groups[(gg, genome)].extend(eval_bc_strs[(genome_calls == genome) & (result.is_nonambient)])

    if len(eval_bcs_arrays) > 0:
        nonambient_summary = pd.DataFrame(OrderedDict([
            ('barcode', np.concatenate(eval_bcs_arrays)),
            ('umis', np.concatenate(umis_per_bc_arrays)),
            ('ambient_loglk', np.concatenate(loglk_arrays)),
            ('pvalue', np.concatenate(pvalue_arrays)),
            ('pvalue_adj', np.concatenate(pvalue_adj_arrays)),
            ('nonambient', np.concatenate(nonambient_arrays)),
            ('genome', np.concatenate(genome_call_arrays)),
        ]))
    else:
        nonambient_summary = pd.DataFrame()

    return filtered_bcs_groups, nonambient_summary




################################################################################
def merge_filtered_metrics(filtered_metrics):
    """ Merge all the barcode filter results and return them as a dictionary """
    result = BarcodeFilterResults(0)
    dresult = {}
    for i, fm in enumerate(filtered_metrics):
        dresult.update(fm.to_dict_with_prefix(i+1))
        # Compute metrics over all gem groups
        result.filtered_bcs += fm.filtered_bcs
        result.filtered_bcs_lb += fm.filtered_bcs_lb
        result.filtered_bcs_ub += fm.filtered_bcs_ub
        result.max_filtered_bcs += fm.max_filtered_bcs
        result.filtered_bcs_var += fm.filtered_bcs_var

    # Estimate CV based on sum of variances and means
    result.filtered_bcs_cv = tk_stats.robust_divide(
        np.sqrt(result.filtered_bcs_var), result.filtered_bcs)

    dresult.update(result.__dict__)
    return dresult


def combine_initial_metrics(genomes, filtered_metrics_groups, genome_filtered_bcs, method, summary):
    # Combine initial-cell-calling metrics
    for genome in genomes:
        # Merge metrics over all gem groups for this genome
        txome_metrics = [v for k,v in filtered_metrics_groups.iteritems() if k[1] == genome]
        txome_summary = merge_filtered_metrics(txome_metrics)

        prefix = genome + '_' if genome != None else ''
        # Append method name to metrics
        summary.update({
            ('%s%s_%s' % (prefix,
                           key,
                           get_filter_method_name(method))): txome_summary[key] \
            for (key,_) in txome_summary.iteritems()})

        summary['%sfiltered_bcs' % prefix] = len(genome_filtered_bcs[genome])

        # NOTE: This metric only applies to the initial cell calls
        summary['%sfiltered_bcs_cv' % prefix] = txome_summary['filtered_bcs_cv']

    return summary


def filter_cellular_barcodes_ordmag(bc_counts, recovered_cells, total_diversity):
    """ Simply take all barcodes that are within an order of magnitude of a top barcode
        that likely represents a cell
    """
    np.random.seed(0)

    if recovered_cells is None:
        recovered_cells = cr_constants.DEFAULT_RECOVERED_CELLS_PER_GEM_GROUP

    metrics = BarcodeFilterResults(0)
    max_filtered_bcs = determine_max_filtered_bcs(total_diversity, recovered_cells)
    metrics.max_filtered_bcs = max_filtered_bcs

    nonzero_bc_counts = bc_counts[bc_counts > 0]
    if len(nonzero_bc_counts) == 0:
        msg = "WARNING: All barcodes do not have enough reads for ordmag, allowing no bcs through"
        return [], metrics, msg

    baseline_bc_idx = int(round(float(recovered_cells) * (1 - cr_constants.ORDMAG_RECOVERED_CELLS_QUANTILE)))
    baseline_bc_idx = min(baseline_bc_idx, len(nonzero_bc_counts) - 1)
    assert baseline_bc_idx < max_filtered_bcs

    # Bootstrap sampling; run algo with many random samples of the data
    top_n_boot = np.array([
        find_within_ordmag(np.random.choice(nonzero_bc_counts, len(nonzero_bc_counts)), baseline_bc_idx)
        for _ in xrange(cr_constants.ORDMAG_NUM_BOOTSTRAP_SAMPLES)
    ])

    metrics.update(summarize_bootstrapped_top_n(top_n_boot))

    # Get the filtered barcodes
    top_n = metrics.filtered_bcs
    top_bc_idx = np.sort(np.argsort(bc_counts)[::-1][0:top_n])
    return top_bc_idx, metrics, None


def filter_cellular_barcodes_fixed_cutoff(bc_counts, cutoff):
    nonzero_bcs = len(bc_counts[bc_counts > 0])
    top_n = min(cutoff, nonzero_bcs)
    top_bc_idx = np.sort(np.argsort(bc_counts)[::-1][0:top_n])
    metrics = BarcodeFilterResults.init_with_constant_call(top_n)
    return top_bc_idx, metrics, None
