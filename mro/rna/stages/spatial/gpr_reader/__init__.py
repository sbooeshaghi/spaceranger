#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
import os
import json
import martian
import cellranger.spatial.utils as spatial_utils

__MRO__ = """
stage GPR_READER(
    in  string  slide_serial_capture_area,
    in  path    loupe_alignment_file,
    in  gpr     gpr_file,
    in  string  barcode_whitelist,
    out json    gpr_data,
    src py      "stages/spatial/gpr_reader",
)
"""

def main(args, outs):
    gpr_reader(args.slide_serial_capture_area, args.loupe_alignment_file, args.gpr_file, args.barcode_whitelist,
               outs.gpr_data)

def gpr_reader(slide_sample_area_id, loupe_alignment_file, gpr_file_in, barcode_whitelist,
               gpr_data_out):
    """
    gpr_reader - entry point for stage to run the Go-based gpr file reader and provide JSON output

    slide_sample_area_id - input slide id of the form <batch_number>-<serial_number>-<area_id>
    loupe_alignment_file - None for automatic pathway or a loupe alignment file pathname
    gpr_file_in          - None (for automatic retrieval or manual pathway) or gpr file pathname
    gpr_data_out         - pathname for JSONified GPR data output to the pipeline
    """

    # files path for this stage
    out_path = os.path.dirname(gpr_data_out)

    # return if no galfile for specified barcode whitelist
    galfile = spatial_utils.get_galfile_path(barcode_whitelist)
    if not os.path.exists(galfile):
        martian.log_info("gpr_reader - no GAL file for whitelist {}".format(barcode_whitelist))
        return

    # Run gprreader with appropriate parameters for various scenarios
    if loupe_alignment_file:
        return

    if slide_sample_area_id:
        slide_sample_id, area_id = spatial_utils.parse_slide_sample_area_id(slide_sample_area_id)
        if gpr_file_in:
            spatial_utils.call_gprreader("read", gpr_file_in, area_id, out_path)
            _, basename_ext = os.path.split(gpr_file_in)
            basename, _ = os.path.splitext(basename_ext)
            gpr_data = spatial_utils.read_from_json(os.path.join(out_path, basename + '_' + area_id + '.json'))
        else:
            spatial_utils.call_gprreader("fetch", slide_sample_id, area_id, out_path)
            gpr_data = spatial_utils.read_from_json(slide_sample_id + '_' + area_id + '.json')
    else: # default case
        galfile = spatial_utils.get_galfile_path(barcode_whitelist)
        spatial_utils.call_gprreader("default", galfile, "B1", out_path)
        gpr_data = spatial_utils.read_from_json('default_B1.json')

    # Regardless of where the GPR is coming from, check it has the required data
    if not gpr_data['spots']['oligo'] or not gpr_data['spots']['fiducial']:
        martian.exit('The slidefile is missing either oligo or fiducial spot information.')

    # Output GPR json
    with open(gpr_data_out, 'w') as json_out:
        json.dump(gpr_data, json_out)
