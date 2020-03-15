#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
from collections import defaultdict
import json
import numpy as np
import pandas as pd
import random
import martian
import tenkit.safe_json as tk_safe_json
import cellranger.chemistry as cr_chem
import cellranger.matrix as cr_matrix
import cellranger.constants as cr_constants
import cellranger.rna.matrix as rna_matrix
import cellranger.rna.library as rna_library
import cellranger.rna.report_matrix as rna_report_mat
import cellranger.utils as cr_utils
from cellranger.csv_utils import write_filtered_barcodes

import cellranger.cell_calling_helpers as helpers
from cellranger.cell_calling_helpers import FilterMethod

FILTER_BARCODES_MIN_MEM_GB = 2.0

__MRO__ = """
stage FILTER_BARCODES(
    in  string sample_id,
    in  h5     matrices_h5,
    in  json   raw_fastq_summary,
    in  json   attach_bcs_summary,
    in  int    recovered_cells,
    in  int    force_cells,
    in  h5     barcode_summary,
    in  csv    barcode_correction_csv,
    in  string barcode_whitelist,
    in  bool   is_antibody_only,
    in  path   reference_path,
    in  int[]  gem_groups,
    in  map    chemistry_def,
    in  json   cell_barcodes          "Cell barcode override",
    out json   summary,
    out csv    filtered_barcodes,
    out csv    aggregate_barcodes,
    out h5     filtered_matrices_h5,
    out path   filtered_matrices_mex,
    out csv    nonambient_calls,
    src py     "stages/counter/filter_barcodes",
) split using (
)
"""


def split(args):
    # We need to store one full copy of the matrix.
    mem_gb = 2 * \
        cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.matrices_h5)
    mem_gb = max(mem_gb, FILTER_BARCODES_MIN_MEM_GB)

    return {
        'chunks': [],
        'join': {
            '__mem_gb': mem_gb,
        }
    }


def main(_args, _outs):
    martian.throw('main is not supposed to run.')


def join(args, outs, _chunk_defs, _chunk_outs):
    filtered_matrix = filter_barcodes(args, outs)

    matrix_attrs = cr_matrix.make_matrix_attrs_count(args.sample_id, args.gem_groups, cr_chem.get_description(args.chemistry_def))
    filtered_matrix.save_h5_file(outs.filtered_matrices_h5, extra_attrs=matrix_attrs, sw_version=martian.get_pipelines_version())

    rna_matrix.save_mex(filtered_matrix,
                        outs.filtered_matrices_mex,
                        martian.get_pipelines_version())


def filter_barcodes(args, outs):
    random.seed(0)
    np.random.seed(0)

    correction_data = pd.read_csv(args.barcode_correction_csv)
    raw_matrix = cr_matrix.CountMatrix.load_h5_file(args.matrices_h5)
    if np.isin(rna_library.ANTIBODY_LIBRARY_TYPE, correction_data.library_type):
        matrix, metrics_to_report, removed_bcs_df = helpers.remove_bcs_with_high_umi_corrected_reads(
            correction_data, raw_matrix)
        # report all idenitified aggregate barcodes, together with their reads, umi corrected reads, fraction of corrected reads, and fraction of total reads
        removed_bcs_df.to_csv(outs.aggregate_barcodes)
        summary = metrics_to_report
    else:
        matrix = raw_matrix
        summary = {}

    if args.cell_barcodes is not None:
        method = FilterMethod.MANUAL
    elif args.force_cells is not None:
        method = FilterMethod.TOP_N_BARCODES
    else:
        if args.is_antibody_only:
            method = FilterMethod.ORDMAG
        else:
            method = FilterMethod.ORDMAG_NONAMBIENT

    summary['total_diversity'] = matrix.bcs_dim
    summary['filter_barcodes_method'] = helpers.get_filter_method_name(method)

    # Get unique gem groups
    unique_gem_groups = sorted(list(set(args.gem_groups)))

    # Get per-gem group cell load
    if args.recovered_cells is not None:
        gg_recovered_cells = int(
            float(args.recovered_cells) / float(len(unique_gem_groups)))
    else:
        gg_recovered_cells = cr_constants.DEFAULT_RECOVERED_CELLS_PER_GEM_GROUP

    if args.is_antibody_only:
        # Only use antibody capture matrix for cell calling
        ab_matrix = matrix.view().select_features_by_type(
            rna_library.ANTIBODY_LIBRARY_TYPE)
        genomes = cr_utils.get_reference_genomes(args.reference_path)
        filtered_metrics_groups, filtered_bcs_groups = helpers.call_initial_cells(ab_matrix, genomes, args.is_antibody_only, unique_gem_groups, method,
                                                                                  gg_recovered_cells, args.cell_barcodes, args.force_cells)
        # Do not do additional cell calling if antibody only
        outs.nonambient_calls = None

        # Record all filtered barcodes
        genome_filtered_bcs = defaultdict(set)
        filtered_bcs = set()
        for (_, genome), bcs in filtered_bcs_groups.iteritems():
            genome_filtered_bcs[genome].update(bcs)
            filtered_bcs.update(bcs)

        # Combine initial-cell-calling metrics
        summary = helpers.combine_initial_metrics(
            genomes, filtered_metrics_groups, genome_filtered_bcs, method, summary)

    else:
        # Only use gene expression matrix for cell calling
        gex_matrix = matrix.view().select_features_by_type(
            rna_library.GENE_EXPRESSION_LIBRARY_TYPE)

        # Make initial cell calls for each genome separately
        genomes = gex_matrix.get_genomes()

        filtered_metrics_groups, filtered_bcs_groups = helpers.call_initial_cells(gex_matrix, genomes, args.is_antibody_only, unique_gem_groups, method,
                                                                                  gg_recovered_cells, args.cell_barcodes, args.force_cells)

        # Do additional cell calling
        outs.nonambient_calls = None
        if method == FilterMethod.ORDMAG_NONAMBIENT:
            # We need the full gene expression matrix instead of just a view
            full_gex_matrix = matrix.select_features_by_type(
                rna_library.GENE_EXPRESSION_LIBRARY_TYPE)
            filtered_bcs_groups, nonambient_summary = helpers.call_additional_cells(full_gex_matrix,
                                                                                    unique_gem_groups,
                                                                                    genomes, filtered_bcs_groups)
            nonambient_summary.to_csv(outs.nonambient_calls)

        # Record all filtered barcodes
        genome_filtered_bcs = defaultdict(set)
        filtered_bcs = set()
        for (_, genome), bcs in filtered_bcs_groups.iteritems():
            genome_filtered_bcs[genome].update(bcs)
            filtered_bcs.update(bcs)

        # Combine initial-cell-calling metrics
        summary = helpers.combine_initial_metrics(
            genomes, filtered_metrics_groups, genome_filtered_bcs, method, summary)

    # Deduplicate and sort filtered barcode sequences
    # Sort by (gem_group, barcode_sequence)
    def barcode_sort_key(x):
        return cr_utils.split_barcode_seq(x)[::-1]

    for genome, bcs in genome_filtered_bcs.iteritems():
        genome_filtered_bcs[genome] = sorted(
            list(set(bcs)), key=barcode_sort_key)
    filtered_bcs = sorted(list(set(filtered_bcs)), key=barcode_sort_key)

    # Re-compute various metrics on the filtered matrix
    reads_summary = cr_utils.merge_jsons_as_dict(
        [args.raw_fastq_summary, args.attach_bcs_summary])
    matrix_summary = rna_report_mat.report_genomes(matrix,
                                                   reads_summary=reads_summary,
                                                   barcode_summary_h5_path=args.barcode_summary,
                                                   recovered_cells=args.recovered_cells,
                                                   cell_bc_seqs=genome_filtered_bcs)

    # Write metrics json
    combined_summary = matrix_summary.copy()
    combined_summary.update(summary)
    with open(outs.summary, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(
            combined_summary), f, indent=4, sort_keys=True)

    # Write the filtered barcodes file
    write_filtered_barcodes(outs.filtered_barcodes, genome_filtered_bcs)

    # Select cell-associated barcodes
    filtered_matrix = matrix.select_barcodes_by_seq(filtered_bcs)

    return filtered_matrix
