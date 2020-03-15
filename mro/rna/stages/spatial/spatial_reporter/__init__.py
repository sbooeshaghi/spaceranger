#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
import os
import cellranger.io as cr_io
import cellranger.websummary.spatial_utils as sp_utils

__MRO__ = """
stage SPATIAL_REPORTER(
    in  json[] summaries,
    in  string sample_id,
    in  string sample_desc,
    in  string slide_serial_info,
    in  path   reference_path,
    in  path   analysis,
    in  h5     barcode_summary_h5,
    in  csv    filtered_barcodes,
    in  string barcode_whitelist,
    in  int[]  gem_groups,
    in  h5     matrix,
    in  json   scalefactors,
    in  float  fraction_under_tissue,
    in  png    tissue_hires_image,
    in  png    tissue_lowres_image,
    in  jpg    aligned_fiducials,
    in  jpg    detected_tissue_image,
    in  csv    tissue_positions_list,
    out path   spatial,
    out html   web_summary,
    out json   metrics_summary_json,
    out csv    metrics_summary_csv,
    src py     "stages/spatial/spatial_reporter",
)
"""

def main(args, outs):
    ### combine all spatial outs into an `spatial` path
    cr_io.makedirs(outs.spatial, allow_existing=True)

    ### copy input arguments destined for outs.spatial
    spatial_file_list = [ args.tissue_hires_image, args.tissue_lowres_image, args.scalefactors,
                  args.detected_tissue_image, args.tissue_positions_list, args.aligned_fiducials ]

    for src_path in spatial_file_list:
        if src_path is not None and os.path.exists(src_path):
            _, src_base = os.path.split(src_path)   # peel off the dirname
            cr_io.copy(src_path, os.path.join(outs.spatial, src_base))

    ### make a web summary
    web_sum_data = sp_utils.create_common_spatial_summaries(args, outs,
                                                            outs.metrics_summary_csv)
    sp_utils.ws_builder.write_html_file(outs.web_summary, web_sum_data)
    return
