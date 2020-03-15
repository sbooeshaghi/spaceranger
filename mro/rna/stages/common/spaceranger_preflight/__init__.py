#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
import martian
import cellranger.preflight as cr_preflight

__MRO__ = """
stage SPACERANGER_PREFLIGHT_LOCAL(
    in map[]  sample_def,
    in string chemistry_name,
    in map    custom_chemistry_def,
    in path   reference_path,
    in bool   check_executables,
    in int    recovered_cells,
    in int    force_cells,
    in string[] allowed_chems,
    in int    r1_length,
    in int    r2_length,
    in path   tissue_image_path,
    in path   loupe_alignment_file,
    in gpr    gpr_file,
    in string slide_serial_capture_area,
    src py    "stages/common/spaceranger_preflight",
)
"""

def run_preflight_checks(args):
    cr_preflight.check_sample_info(args.sample_def, args.reference_path)
    cr_preflight.check_common_preflights(args.reference_path,
                            args.chemistry_name, args.custom_chemistry_def, args.allowed_chems,
                            args.r1_length, args.r2_length, args.check_executables,
                            args.recovered_cells, args.force_cells)
    cr_preflight.check_spatial_arguments(args.tissue_image_path,
                                         args.loupe_alignment_file,
                                         args.gpr_file,
                                         args.slide_serial_capture_area
                                         )
    cr_preflight.check_image_open(args.tissue_image_path, args.loupe_alignment_file)

def main(args, outs):
    try:
        run_preflight_checks(args)
    except cr_preflight.PreflightException as e:
        martian.exit(e.msg)

    cr_preflight.record_package_versions()
