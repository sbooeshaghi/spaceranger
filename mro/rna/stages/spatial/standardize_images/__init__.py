#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
import cellranger.feature.utils as feature_utils
import cellranger.spatial.utils as spatial_utils
import martian

__MRO__ = """
stage STANDARDIZE_IMAGES(
    in path     tissue_image_path,
    in path     spot_image_path,
    out json    scalefactors_json,
    out png     tissue_hires_image,
    out png     tissue_lowres_image,
    out png     spot_hires_image,
    src py      "stages/spatial/standardize_images",
) split using (
)
"""
HIRES_MAX_DIM = 2000
LORES_MAX_DIM = 600

def split(args):
    # estimate downsampling bytes
    num_tissue_gb = 1
    num_spots_gb = 1
    if args.tissue_image_path:
        num_tissue_gb = spatial_utils.call_tiffer_mem_estimate_gb(args.tissue_image_path, HIRES_MAX_DIM)
    if args.spot_image_path:
        num_spots_gb = spatial_utils.call_tiffer_mem_estimate_gb(args.spot_image_path, LORES_MAX_DIM)
    num_gbs = num_tissue_gb
    if num_spots_gb > num_tissue_gb:
        num_gbs = num_spots_gb

    # add one for overhead
    num_gbs = num_gbs + 3
    return {
        'chunks': [],
        'join': {
            '__mem_gb': int(num_gbs),
            '__vmem_gb': int(num_gbs) + 6, # martian default is +3
        }
    }

def main(_args, _outs):
    martian.throw('main is not supposed to run.')

def join(args, outs, _chunk_defs, _chunk_outs):
    standardize_images( args.tissue_image_path, args.spot_image_path,
            		outs.scalefactors_json,
           	 	outs.tissue_hires_image,
			outs.tissue_lowres_image,
            		outs.spot_hires_image)

def standardize_images( tissue_image_path, spot_image_path,
        		scalefactors_json,
        		tissue_hires_image,
			tissue_lowres_image,
        		spot_hires_image):
    """ Read tissue and (optionally) spot images, and return their downsampled copies """

    # Downsample input tissue image by calling tiffer
    tissue_hires_json = spatial_utils.call_tiffer_resample(tissue_image_path, HIRES_MAX_DIM, tissue_hires_image)
    tissue_lowres_json = spatial_utils.call_tiffer_resample(tissue_image_path, LORES_MAX_DIM, tissue_lowres_image)

    # Create a dictionary to store scaling factors
    scalefactors = {}
    scalefactors.update( { "tissue_hires_scalef": tissue_hires_json['scale'],
                           "tissue_lowres_scalef": tissue_lowres_json['scale']} )

    # If there is a spot image, downsample that as well
    if spot_image_path is not None:
        spot_hires_json = spatial_utils.call_tiffer_resample(spot_image_path, HIRES_MAX_DIM, spot_hires_image)
        scalefactors.update( {"spot_hires_scalef": spot_hires_json['scale']} )

    feature_utils.write_json_from_dict(scalefactors, scalefactors_json)
