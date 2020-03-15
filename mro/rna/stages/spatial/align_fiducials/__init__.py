#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
import json
from cellranger.spatial.fiducial_alignment import align_fiducials_5k, align_fiducials_from_override
import cellranger.spatial.utils as spatial_utils
import cellranger.metrics_names as metrics_names
import martian
import numpy as np

__MRO__ = """
stage ALIGN_FIDUCIALS(
    in  jpg  downsampled_tissue_image,
    in  json scalefactors_json,
    in  json gpr_data,
    in  path loupe_alignment_file,
    in  string transform_method,
    out txt  spot_positions_list,
    out txt  fiducial_positions_list,
    out jpg  aligned_fiducials,
    out jpg  detected_keypoints,
    out json scalefactors_json,
    out json alignment_metrics,
    src py   "stages/spatial/align_fiducials",
)
"""

def main(args, outs):
    metrics = align_fiducials(args.loupe_alignment_file, args.downsampled_tissue_image,
                    args.scalefactors_json, args.gpr_data, args.transform_method,
                    outs.spot_positions_list, outs.fiducial_positions_list,
                    outs.aligned_fiducials, outs.detected_keypoints, outs.scalefactors_json)

    with open(outs.alignment_metrics, 'w') as f:
        json.dump(metrics, f, indent=4)

def align_fiducials(loupe_alignment_file, tissue_hires_image, scalefactors_json, gpr_data,
                    transform_method, spot_positions_list, fiducial_positions_list,
                    aligned_fiducials, detected_keypoints, scalefactors_json_out):

    ### Read the downsampled tissue image
    img = spatial_utils.cv_read_image_standard(tissue_hires_image)
    scalefactors = json.load(open(scalefactors_json, 'r'))
    iscale = scalefactors["tissue_hires_scalef"]

    ### Write out metrics for automated alignment
    metrics = {}
    if loupe_alignment_file:
        with open(loupe_alignment_file, 'r') as gprfile:
            gpr_data = json.load(gprfile)
        fid_dia, spot_dia = align_fiducials_from_override(img, iscale, gpr_data,
                                      spot_positions_list, fiducial_positions_list,
                                      aligned_fiducials)
        metrics["alignment_method"] = 'Manual Alignment'
    else:
        with open(gpr_data,'r') as gprfile:
            gpr_data = json.load(gprfile)
        spots = gpr_data['spots']

        try:
            reg_metrics, (fid_dia, spot_dia) = align_fiducials_5k(img, iscale, spots, transform_method,
                                   spot_positions_list, fiducial_positions_list,
                                   aligned_fiducials, detected_keypoints)
            metrics["num_fiducials_detected"] = reg_metrics.nfids
            metrics["alignment_method"] = 'Automatic Alignment'
            metrics["residual_spot_error"] = reg_metrics.mse if not np.isnan(reg_metrics.mse) else -1
            metrics["alignment_sigma2"] = reg_metrics.sigma2
            msg = reg_metrics.test()
            metrics[metrics_names.SUSPECT_ALIGNMENT] = ( msg is not None )
            if msg:
                martian.log_info(msg)
        except RuntimeError, e:
            martian.exit(e.message)

    # scale diameters to full resolution and write out
    fid_dia /= iscale
    spot_dia /= iscale

    spatial_utils.update_scalefactors_for_diameters(scalefactors_json, fid_dia, spot_dia, scalefactors_json_out)

    return metrics
