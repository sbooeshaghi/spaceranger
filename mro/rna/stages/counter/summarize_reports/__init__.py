#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import martian

import cellranger.io as cr_io
import cellranger.report as cr_report
import cellranger.utils as cr_utils
import cellranger.webshim.common as cr_webshim
import cellranger.websummary.sample_properties as sp
from cellranger.webshim.constants.shared import PIPELINE_COUNT
from cellranger.websummary.web_summary_builder import build_web_summary_html_sc

__MRO__ = """
stage SUMMARIZE_REPORTS(
    in  json[] summaries,
    in  string sample_id,
    in  string sample_desc,
    in  path   reference_path,
    in  path   analysis,
    in  h5     barcode_summary_h5,
    in  h5     filtered_gene_bc_matrices_h5,
    in  csv    filtered_barcodes,
    in  csv    feature_reference,
    in  string barcode_whitelist,
    in  int[]  gem_groups,
    out json   metrics_summary_json,
    out csv    metrics_summary_csv,
    out html   web_summary,
    out csv    feature_reference,
    src py     "stages/counter/summarize_reports",
) using (
    mem_gb = 4,
)
"""

def main(args, outs):
    id_dict = {
        "sample_id" : args.sample_id,
        "sample_desc" : args.sample_desc
    }
    cr_report.merge_jsons(args.summaries, outs.metrics_summary_json, [id_dict])

    sample_data_paths = sp.SampleDataPaths(
        summary_path=outs.metrics_summary_json,
        barcode_summary_path=args.barcode_summary_h5,
        analysis_path=args.analysis,
        filtered_barcodes_path=args.filtered_barcodes,
    )

    genomes = cr_utils.get_reference_genomes(args.reference_path)
    sample_properties = sp.ExtendedCountSampleProperties(sample_id=args.sample_id,
                                                         sample_desc=args.sample_desc,
                                                         barcode_whitelist=args.barcode_whitelist,
                                                         reference_path=args.reference_path,
                                                         genomes=genomes)

    # TODO: Move metrics CSV somewhere else
    sample_data = cr_webshim.load_sample_data(sample_properties, sample_data_paths)
    cr_webshim.build_metrics_summary_csv(outs.metrics_summary_csv, sample_properties, sample_data, PIPELINE_COUNT)

    build_web_summary_html_sc(outs.web_summary, sample_properties, sample_data_paths, PIPELINE_COUNT)

    if args.feature_reference is not None:
        cr_io.copy(args.feature_reference, outs.feature_reference)
    return
