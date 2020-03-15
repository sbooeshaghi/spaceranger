#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

import cellranger.io as cr_io
import martian
import tempfile
import subprocess
import os
import tenkit.log_subprocess as tk_subproc

import cellranger.spatial.utils as spatial_utils

__MRO__ = """
stage CLOUPE_TILE_IMAGES(
    in  path tissue_image_path,
    in  int  tile_size,
    in  bool skip_stage,
    in  bool no_secondary_analysis,
    out json dzi_info,
    out path dzi_tile_path,
    src py   "stages/cloupe/cloupe_tile_images",
) split using (
)
"""
TILE_DATASET_NAME = "tiles"

def do_not_make_tiles(args):
    """
    Returns True if there is a reason why this stage should not make tiles
    """
    if args.no_secondary_analysis:
        martian.log_info("Skipping tile generation by instruction (--no-secondary-analysis)")
        return True
    if args.skip_stage:
        martian.log_info("Skipping tile generation as stage disabled")
        return True
    if args.tissue_image_path is None:
        martian.log_info("Skipping tile generation due to missing tissue image path")
        return True
    if not os.path.exists(args.tissue_image_path):
        martian.log_info("Skipping tile generation due to incorrect tissue image path")
        return True
    return False

def split(args):
    # low mem usage if skipped
    if do_not_make_tiles(args):
        return {'chunks': [], 'join': {'__mem_gb': 1}}

    mem_gb = spatial_utils.call_tiffer_mem_estimate_gb(args.tissue_image_path, "tile") + 1 # account for possible overhead
    # address space requirement is a bit more than usual relative to
    # memory usage because of the golang process.
    vmem_gb = mem_gb + 5
    return {
        'chunks': [],
        'join': {
            '__mem_gb': int(mem_gb),
            '__vmem_gb': int(vmem_gb),
        },
    }

def join(args, outs, chunk_defs, chunk_outs):
    if do_not_make_tiles(args):
        outs.dzi_info = None
        outs.dzi_tiles_path = None
        return

    tmp_outs = tempfile.mkdtemp()

    call = ["cloupetiler",
            args.tissue_image_path,
            tmp_outs,
            TILE_DATASET_NAME,
            "--tilesize", str(args.tile_size)]

    # this is copied from crconverter; some paths may require
    # encoding as a utf-8 byte string for use in check_output.
    call_bytes = [arg.encode('utf-8') for arg in call]

    # but keep the arg 'call' here because log_info attempts to
    # decode the message if we give it a byte string rather than
    # a unicode string.
    martian.log_info(u"Running cloupetiler: %s" % u" ".join(call))
    try:
        results = tk_subproc.check_output(call_bytes, stderr=subprocess.STDOUT)
        martian.log_info("cloupetiler output: %s" % results)
        generated_dzi_path = os.path.join(tmp_outs, "%s.dzi" % TILE_DATASET_NAME)
        cr_io.move(generated_dzi_path, outs.dzi_info)
        cr_io.move(os.path.join(tmp_outs, "%s_files" % TILE_DATASET_NAME), outs.dzi_tiles_path)
    except subprocess.CalledProcessError, e:
        outs.dzi_info = None
        outs.dzi_tiles_path = None
        martian.throw("Could not generate image tiles: \n%s" % e.output)
