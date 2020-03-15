#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
import os
import json
import cv2
import numpy as np
import martian
import subprocess
import tenkit.log_subprocess as tk_subproc
import cellranger.constants as cr_constants

def shrink_to_max(shape, target):
    """Given an image shape (x,y), if the largest dimension is larger than
    target, return a scalefactor to shrink the image to approximately target in that dimension,
    otherwise return 1.0"""

    scale = 1.0 * target / max(shape)
    return min(1.0, scale)

def cv_read_image_standard(filename, maxval=255):
    """ Canonical function to read image files """

    img = cv2.imread(filename, cv2.IMREAD_UNCHANGED)
    if img is not None:
        if len(img.shape) == 3:
            img = cv2.cvtColor(img,cv2.COLOR_RGB2GRAY)
        cv2.normalize(img, img, maxval, 0, cv2.NORM_MINMAX)
        if img.dtype=='uint16':
            img = img.astype('uint8')

    return img

def cv_write_downsampled_image(cvimg, filename, maxres):
    """ Scale, resize, and return the downsampled image """

    scalef = shrink_to_max( cvimg.shape, maxres )
    rcvimg = cv2.resize(cvimg, (0,0), fx=scalef, fy=scalef)
    params = [ cv2.IMWRITE_JPEG_QUALITY, 80 ]
    cv2.imwrite(filename, rcvimg, params)

    return scalef

def cv_composite_labelmap(background, labelmap, colormap, alpha):
    """ Take a background image (generally grayscale) and composite with
        a color image based on a labelmap

        alpha value is [0,1.0)  (note open - can't be 1.0)

        example:
        cv_composite_labelmap(gray, mask, {0:(255,0,0), 1:(0,0,255), 0.75})

        takes a binary labelmap (mask) and composites it with the gray image
        to return a color image with 0.25 (25%) background and 0.75 (75%) labelmap
        where the zero mask values are red (255,0,0) and the one mask values are
        blue (0,0,255)
        """
    if background.shape[:2] != labelmap.shape[:2]:
        raise RuntimeError("cv_composite_labelmap called with incompatible images")

    if alpha < 0.0 or alpha >= 1.0:
        raise ValueError("cv_composite_labelmap passed illegal alpha value - should be [0.,1.)")

    if background.ndim < 3:
        output = cv2.cvtColor(background,cv2.COLOR_GRAY2RGB)
    else:
        output = background.copy()

    for label,color in colormap.iteritems():
        mask = (labelmap == label)
        output[mask,:] = (1.0-alpha)*output[mask,:] +  np.multiply(alpha,color)

    return output

def write_detected_keypoints_image(img, keypoints, diameter, filename):
    """ Draw the detected keypoints over a brightfield image """

    output = img.copy()
    for kp in keypoints:
        col, row = int(kp.pt[0]), int(kp.pt[1])
        cv2.circle(output, (col,row), int( diameter / 2.0 ), (0,255,0), 3 )
    cv2.imwrite(filename, output)

def write_aligned_fiducials_image(img, aligned_fiducials, diameter, filename):
    """ Draw the aligned fiducials over a brightfield image """
    output = img.copy()
    for spot in aligned_fiducials:
        x, y = spot
        col, row = int(x), int(y)
        embiggen = 1.2  # add 20% to spot size to be able to visualize underlying fiducial
        cv2.circle(output, (col,row), int( diameter / 2.0 * embiggen + 0.5), (0,0,255), 2 )
    cv2.imwrite(filename, output)

def update_scalefactors_for_diameters(in_json_path, fid_dia, spot_dia, out_json_path):
    """serialize (optional) fiducial diameter and spot diameter to scalefactors JSON"""
    scalefactors=json.load(open(in_json_path,'r'))

    if fid_dia is not None:
        scalefactors["fiducial_diameter_fullres"] = fid_dia
    if spot_dia is not None:
        scalefactors["spot_diameter_fullres"] = spot_dia

    json.dump(scalefactors, open(out_json_path,'w'))


def transform_oligo_locations(spots, destination_points, source_points, inverse_rotation_mtx, inverse_t_matrix, scale):
    """ Given oligo locations from a GPR file, transform their x/y coordinates using the outputs of fiducial alignment """

    oligo_coordinates = np.array([ [ oligo['row'], oligo['col'], oligo['x'], oligo['y'],  ] for oligo in spots['oligo'] ])
    spots_xy = oligo_coordinates[:,2:4].astype('double')
    spots_centered = spots_xy - np.median(destination_points, 0)

    transformed_xy = np.dot(spots_centered + np.tile(inverse_t_matrix, (np.shape(spots_centered)[0], 1)), inverse_rotation_mtx)
    uncentered_xy = transformed_xy + np.median(source_points, 0)
    uncentered_xy = np.vstack([uncentered_xy[:,1], uncentered_xy[:,0]]).T
    scaled_xy = uncentered_xy / scale
    transformed_coordinates = np.hstack([oligo_coordinates[:,0:2], scaled_xy]).astype(int)

    return transformed_coordinates

def write_to_list(lst, coordinates):
    """ Write given coordinates to a txt file """

    with open(lst, 'w') as fd:
        for spot in coordinates:
            fd.write( ",".join( map(str, spot)))
            fd.write( "\n")

def call_tiffer_mem_estimate(image_path, downsample_size):
    """ Run tiffer memreq with given parameters """
    call = ["tiffer",
        "memreq",
        image_path,
        "--size",str(downsample_size)]
    proc=tk_subproc.Popen(call,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    (proc_stdout, proc_stderr) = proc.communicate(None)
    if proc.returncode != 0:
        msgs=["Could not generate image downsample memory estimate:", proc_stdout, proc_stderr ]
        martian.exit("\n".join([x for x in msgs if x]))
    output_bytes = int(proc_stdout.strip())
    return output_bytes

def call_tiffer_mem_estimate_gb(image_path, downsample_size):
    """ call_tiffer_mem_estimate but in GB """
    output_bytes = call_tiffer_mem_estimate(image_path, downsample_size)
    bytes_gb = np.ceil(output_bytes/(1024.0*1024.0*1024.0))
    return bytes_gb

def call_tiffer_resample(image_path, downsample_size, full_output_path):
    """ Run tiffer with given parameters """

    full_output_path = os.path.abspath(full_output_path)
    output_path = os.path.dirname(full_output_path)
    output_file_name = os.path.basename(full_output_path)
    (output_file_name_noext, ext) = os.path.splitext(output_file_name)
    ext = ext[1:]

    call = ["tiffer",
            "resample",
            image_path,
            output_path,
            "--jsonfile",
            "--size", str(downsample_size),
            "--format", ext,
            "--outname", output_file_name_noext]

    try:
        _ = tk_subproc.check_output(call, stderr=subprocess.STDOUT)
        # ^ the unused returned variable is warning/error messages.
        output_json_wpath = os.path.join(output_path, '{}.json'.format(output_file_name_noext))
        output_json = json.load(open(output_json_wpath, 'r'))
    except subprocess.CalledProcessError, e:
        martian.exit("Could not generate downsampled image: \n%s" % e.output)
    return output_json

def call_gprreader(mode, input_file, well, output_path):
    """ Run gprreader with given parameters """

    call = ["gprreader",
             mode,
	     input_file,
             output_path,
             "--area",well]

    try:
        tk_subproc.check_output(call, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError, e:
        if mode == "fetch":
            martian.exit("Could not retrieve spot layout data: \n%s" % e.output )
        else:
            martian.exit("Could not read spot layout file: \n%s" % e.output )

def parse_slide_sample_area_id(slide_sample_area_id):
    """ Given an input to the pipeline like V19L01-006-B1,
        parse out slide sample id and area id """

    slide_sample_id, area_id = slide_sample_area_id[:-3], slide_sample_area_id[-2:]
    return slide_sample_id, area_id

def get_galfile_path(barcode_whitelist):
    """ Given a barcode whitelist, return the path to the corresponding GAL file """

    path_to_galfile = os.path.join(cr_constants.BARCODE_WHITELIST_PATH, barcode_whitelist + '.gal')
    return path_to_galfile

def read_from_json(filename):
    """ Read from a given json file """

    with open(filename) as json_file:
        gpr_data = json.load(json_file)
    return gpr_data
