#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

import numpy as np
import cv2
import scipy
import skimage
import skimage.morphology
import skimage.feature
import skimage.filters
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import cellranger.barcodes.utils as bc_utils
import cellranger.spatial.utils as spatial_utils
import martian

def get_mask(original, plot=False, plot_title=None, plotfname=None, bounding_box=None):
    """
    The get_mask function takes an image (array) in grayscale and uses the opencv grabcut
    algorithm to detect tissue section(s) (foreground) on the glass slide surface (background).
    Markers for initialization of the grabcut algorithm are based on otsu tresholding of the
    grayscale and the gradient (local max-min) of the image. Returns a binary mask (where 1=tissue,
    0=background) as well as a qc image with the mask overlayed on the input image.
    """

    if len(original.shape) != 2:
        raise RuntimeError("non-2D image (color?) passed to get_mask nD={}".format(len(original.shape)))

    if bounding_box is None:
        height, width = original.shape[:2]
        # order matters below - should trace a bounding box by adjacent edges
        bounding_box = np.array([ (0,0), (width-1,0), (width-1, height-1), (0, height-1) ])

    gray = box(original,bounding_box)

    longest_side = max(original.shape[:2])
    # set sizes of objects and holes (assumes approx square input image)
    small_holes_size = int(longest_side/2.0)
    small_objects_size = int(longest_side/2.0)
    large_holes_size = int(longest_side*50.0)

    # Calculate grayscale intensity Otsu treshold, often good proxy for where tissue is
    # (not as good when background and tissue has similar gray scale histograms)
    otsu_treshed = gray <= skimage.filters.threshold_otsu(gray)
    # Remove holes and tissue debris
    otsu_treshed = skimage.morphology.remove_small_objects(otsu_treshed, small_objects_size)
    otsu_treshed = skimage.morphology.remove_small_holes(otsu_treshed, small_holes_size)

    # Get the gradient (local max - local min) of the gray scale image,
    # high gradient usually indicates tissue (not as good for when larger
    # areas are out of focus or tissue is really smooth potentially
    # affected by low resolution images)
    gradient = skimage.filters.rank.gradient(gray, skimage.morphology.disk(5))

    # Binarize the gradient into two classes 1=FG and 0=BG using Otsu treshold
    inverted_grad = skimage.util.invert(gradient, signed_float=False)
    otsu_of_gradient = inverted_grad <= skimage.filters.threshold_otsu(inverted_grad)
    otsu_of_gradient = skimage.morphology.remove_small_objects(otsu_of_gradient, small_objects_size)
    otsu_of_gradient = skimage.morphology.remove_small_holes(otsu_of_gradient, small_holes_size)

    # Detect canny edges on the grayscale image (many edges usually indicate tissue)
    canny_edges = skimage.feature.canny(gray)
    closed_canny = skimage.morphology.closing(canny_edges)
    closed_canny = scipy.ndimage.distance_transform_edt(~closed_canny) <= longest_side*0.01

    # Sum upp the two estimates of tissue placement
    # (gradient based and Outsu on grayscale intensity)
    BACKGROUND = 0
    DETECTED_BY_ONE_METHOD = 1
    DETECTED_BY_TWO = 2
    DETECTED_BY_ALL = 3
    otsu_sum = np.add(np.add(otsu_of_gradient.astype('uint8'),
                             otsu_treshed.astype('uint8')),
                      closed_canny.astype('uint8'))

    # Start making markers for the grabcut
    markers_gc = np.zeros(gray.shape).astype('uint8')
    # to keep track of not yet classed vs obvious background
    classed = np.zeros(otsu_sum.shape).astype('uint8')

    ##### below is order dependent based on priority, pixels may be assign GC_BGD early and then
    ##### the same pixel may get GC_FGD later

    # If classed as background by both methods add a margin of 1% image longest side (in pixels) and
    # set to an obvious background pixels
    background = np.zeros(otsu_sum.shape).astype('uint8')
    background[otsu_sum == BACKGROUND] = 1
    background = scipy.ndimage.distance_transform_edt(background) >= longest_side*0.01
    markers_gc[background == 1] = cv2.GC_BGD
    classed[background == 1] += 1

    # Take the two estimates (otsu_sum) fill all holes and set everything detected by at least one
    # method to be probable Background (instead of obvious background)
    # This is done so no holes will be classed as obvious background
    no_holes = np.zeros(otsu_sum.shape).astype('bool')
    no_holes[otsu_sum >= DETECTED_BY_ONE_METHOD] = True
    # remove_small_holes treats 0/1 mask different than false/true - use boolean
    no_holes = skimage.morphology.remove_small_holes(no_holes, large_holes_size)
    markers_gc[no_holes >= 1] = cv2.GC_PR_BGD
    classed[no_holes >= 1] += 1

    # If detected by at least one method set to be a possible foreground pixel
    markers_gc[otsu_sum >= DETECTED_BY_ONE_METHOD] = cv2.GC_PR_FGD
    classed[otsu_sum >= DETECTED_BY_ONE_METHOD] += 1

    # If detected by two methods add a margin of 5% (inward) image longest side (in pixels)
    # basically make the estimate smaller by some amount around the boundaries
    # set as an obvious foreground (object) pixel
    foreground = np.zeros(otsu_sum.shape).astype('uint8')
    foreground[otsu_sum == DETECTED_BY_TWO] = 1
    foreground = scipy.ndimage.distance_transform_edt(foreground) >= longest_side*0.05
    markers_gc[foreground == 1] = cv2.GC_FGD
    classed[foreground == 1] += 1

    # If detected by all methods add a margin of 2.5% image longest side (in pixels)
    # same as above, but with smaller margin of safety due to greater certainty
    # set as an obvious foreground (object) pixel
    foreground = np.zeros(otsu_sum.shape).astype('uint8')
    foreground[otsu_sum == DETECTED_BY_ALL] = 1
    foreground = scipy.ndimage.distance_transform_edt(foreground) >= longest_side*0.025
    markers_gc[foreground == 1] = cv2.GC_FGD
    classed[foreground == 1] += 1

    # Within tissue estimates (no_holes) but zero in outsu sum should be probable background
    # Essentially encourage interior holes when no method has indicated foreground.
    otsu_background = np.zeros(otsu_sum.shape).astype('uint8')
    otsu_background[otsu_sum == BACKGROUND] = 1
    probable_foreground_hole = np.add(otsu_background.astype('uint8'), no_holes.astype('uint8'))
    temp = np.ones(otsu_sum.shape).astype('uint8')
    temp[probable_foreground_hole == 2] = 0
    temp = scipy.ndimage.distance_transform_edt(temp) >= longest_side*0.015
    temp = skimage.util.invert(temp, signed_float=False)
    probable_foreground_hole = temp.astype('uint8')
    markers_gc[temp == 1] = cv2.GC_PR_BGD

    # Set any unclassed pixels to be possible background
    markers_gc[classed == 0] = cv2.GC_PR_BGD

    # if statement for switching creation of plots for debugging and evaluation on and off
    if plot:
        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=[15, 15.5])
        if plot_title:
            fig.suptitle(plot_title)

        axes[0][0].imshow(unbox(original,
                                bounding_box,
                                gray,
                                border_thickness=5),
                          cmap='gray',
                          interpolation='nearest')
        axes[0][0].set_title('Input Image')
        if bounding_box is not None:
            bbox = Polygon(bounding_box, closed=True, alpha=0.25, facecolor='red')
            axes[0][0].add_patch(bbox)

        axes[0][1].imshow(canny_edges,
                          cmap=plt.cm.gray,
                          interpolation='nearest')
        axes[0][1].set_title('Canny Edges')

        axes[0][2].imshow(otsu_treshed,
                          cmap=plt.cm.gray,
                          interpolation='nearest')
        axes[0][2].set_title('Otsu treshold')

        axes[1][0].imshow(closed_canny,
                          cmap=plt.cm.gray,
                          interpolation='nearest')
        axes[1][0].set_title('Closed Canny')

        axes[1][1].imshow(gradient,
                          cmap=plt.cm.gray,
                          interpolation='nearest')
        axes[1][1].set_title('Local Gradient')

        axes[1][2].imshow(otsu_of_gradient,
                          cmap=plt.cm.gray,
                          interpolation='nearest')
        axes[1][2].set_title('Otsu of Gradient')

        axes[2][0].imshow(otsu_sum,
                          cmap=plt.cm.tab10,
                          interpolation='nearest')
        axes[2][0].set_title('Sum of Otsus')

        axes[2][1].imshow(markers_gc,
                          cmap=plt.cm.tab10,
                          interpolation='nearest')
        axes[2][1].set_title('Grabcut Markers')

    # make the image compatible with the grabcut (color)
    im = cv2.cvtColor(gray.astype('uint8'), cv2.COLOR_GRAY2RGB)

    # run the opencv grabcut
    bgmodel = np.zeros((1, 65), dtype='float64')
    fgmodel = np.zeros((1, 65), dtype='float64')
    mask, bgmodel, fgmodel = cv2.grabCut(im, markers_gc, None,
                                         bgmodel, fgmodel, 6, cv2.GC_INIT_WITH_MASK)
    mask = np.where((mask == 2)|(mask == 0), 0, 1)

    # generate qc image as 24-bit RGB instead of crazy float64 RGB
    qc = im
#    qc = spatial_utils.cv_composite_labelmap(qc, mask, {0:(255,0,0), 1:(0,0,255)}, alpha=0.25)
    qc = spatial_utils.cv_composite_labelmap(qc, mask, {1:(0,0,255)}, alpha=0.25)   # try skipping background coloring

    # if enabled add the qc to the debugging/evaluation plot and save/show
    if plot:
        axes[2][2].imshow(qc)
        axes[2][2].set_title('GrabCut')

        plt.subplots_adjust(left=None, bottom=0, right=None, top=0.95, wspace=0.1, hspace=0.1)

        if plotfname:
            plt.savefig(plotfname)
        else:
            plt.show()

        plt.close()

    # reembed the mask in the original image shape
    mask = unbox(np.zeros(original.shape[:2],
                          dtype=np.uint8),
                 bounding_box,
                 mask,
                 border_thickness=0, interp=cv2.INTER_NEAREST)
    # reembed the qc image in the original
    qc = unbox(original,
               bounding_box,
               qc,
               border_thickness=2, interp=cv2.INTER_NEAREST)

    return mask, qc

def get_manual_qc_image(image, plot=False, bounding_box=None):
    """
    Prepare a QC image to indicate where the user specified tissue vs background,
    instead of doing automatic tissue detection.
    """
    if bounding_box:
        y_pos, x_pos, height, width = bounding_box
    else:
        y_pos = 0
        x_pos = 0
        height, width = image.shape[:2]

    original = image
    image = original[y_pos:y_pos+height, x_pos:x_pos+width]

    # just show the border
    qc = unbox(original,
               [y_pos, x_pos, height, width],
               image,
               border_thickness=2)
    return qc


def _get_bbox_transform(bounding_box):
    len1=np.linalg.norm( bounding_box[1] - bounding_box[0] )
    len2=np.linalg.norm( bounding_box[2] - bounding_box[1] )

    cols,rows = int(len1+0.5), int(len2+0.5)
    srcpts=np.float32(bounding_box[:3])
    dstpts=np.float32(np.array( [ [0,0], [cols-1,0], [cols-1,rows-1] ] ))

    trans=cv2.getAffineTransform(src=srcpts,dst=dstpts)

    return trans, cols, rows

def box(original, bounding_box, interp=cv2.INTER_LINEAR):
    """
    Given a bounding box comprising coordinates of vertices (NOT x,y,w,h), extract a possibly rotated
    rectangle from the original image and return the subimage
    """

    trans, cols, rows = _get_bbox_transform(bounding_box)
    return cv2.warpAffine(original,trans,(cols,rows),borderMode=cv2.BORDER_REPLICATE, flags=interp)

def _unbox1(original, bounding_box, boxed_image, interp=cv2.INTER_LINEAR):
    sentinel = np.iinfo(boxed_image.dtype).max
    img_max = np.max(boxed_image)
    if img_max == sentinel:
        cv2.normalize(boxed_image, boxed_image, sentinel-1, 0, cv2.NORM_MINMAX)

    trans, _, _ = _get_bbox_transform(bounding_box)
    output=cv2.warpAffine(src=boxed_image,M=trans,flags=cv2.WARP_INVERSE_MAP|cv2.BORDER_CONSTANT|interp,dsize=(original.shape[1],original.shape[0]),borderValue=sentinel)
    output[output==sentinel] = original[output==sentinel]
    return output

def unbox(original, bounding_box, boxed_image, border_thickness=None, interp=cv2.INTER_LINEAR):
    """
    Function for re-embedding a cropped bounding boxed image back in the original image (or shape)

    If border thickness is > 0, then we return a color version of what otherwise MUST be grayscale inputs.
    """

    if len(boxed_image.shape) == 3 and len(original.shape) == 2:
        original = cv2.cvtColor(original, cv2.COLOR_GRAY2RGB)
        output = np.zeros( original.shape, dtype='uint8' )
        for i in range(3):
            output[:,:,i] = _unbox1(original[:,:,i], bounding_box, boxed_image[:,:,i], interp)
    else:
        output = _unbox1(original, bounding_box, boxed_image, interp)

    # add a blue border
    if border_thickness:
        if border_thickness < 0:
            raise ValueError("negative border thickness passed to unbox")

        if len(output.shape) == 2:
            output=cv2.cvtColor(output,cv2.COLOR_GRAY2RGB)

        cv2.drawContours(output, [np.intp(bounding_box)], 0, color=(0,0,255), thickness=2)

    return output

def _pad_rotated_rectangle(rrect, padding):
    """
    given an opencv rotated rectangle, add "padding" pixels to all sides and return a reconstructed
    RotatedRect
    """
    (centerx,centery),(width,height),angle = rrect
    width += 2*padding
    height += 2*padding

    return ((centerx,centery),(width,height),angle)

def get_bounding_box(spots, scale, padding):
    """
    Generates a rectangular bounding box based on a set of points WITHOUT the assumption that
    the box is axis-aligned.  The box is padded on each side by the amount padding specified in
    original (unscaled) image units.  The box returned is scaled.  Returns the (x,y) coordinates of
    four corners.  Makes heavy use of OpenCV convenience functions.

    Args:
        spots (iterable): A list of spots centre coordinates (row, col, y_pixel, x_pixel).
        scale (float): The float used to scale the spot coordinates to the (usu reduced-size) input image
        padding (int): The padding to be used (generally spot diameter)

    Returns:
        bounding_box (list): box coordinates [ (x1,y1)...(x4,y4) ] as float
    """
    x=np.array([x for _,_,_,x in spots])    # spots are row,col,y,x so we're taking x,y for openCV
    y=np.array([y for _,_,y,_ in spots])

    # compute the "rotated rectangle" surrounding the (padded) array - no assumption of axis alignment
    rrect=cv2.minAreaRect(np.column_stack((x,y)))
    rrect=_pad_rotated_rectangle(rrect,padding)

    # compute the corners of the rotated rectangle
    #
    # Given points 0, 1, 2, 3, these define four sides 0:1, 1:2, 2:3, and 3:0.
    # OpenCV code guarantees that 0:1 is opposite 2:3, 1:2 is opposite 3:0, and
    # 0:1 is adjacent to 1:2.  This implies that the lengths defined by 0:1 and 1:2
    # give the two side lengths of the rectangle.  Also, drawing the sides "in order"
    # will trace out a continuous contour.
    bbox=cv2.boxPoints(rrect)

    # scale the corners
    sbox=np.round(np.multiply(scale, bbox))

    return sbox



def read_barcode_coordinates(barcode_whitelist):
    """
    Function that reads the barcode_whitelist file and returns a dict with spatial coordinates
    (row and column as integers) as keys and barcode sequnce as values.

        Args:
            barcode_whitelist (str): Name of the barcode whitelist.

        Returns:
            barcodes (dict): Dictionary with (row, column) as keys and barcode sequences as values.
    """

    # read barcode coordinates
    barcodes = dict()
    fn = bc_utils.barcode_whitelist_path(barcode_whitelist + "_coordinates")
    for line in open(fn, 'r'):
        if '#' in line:
            continue
        fields = line.strip().split()
        c, r = int(fields[1])-1, int(fields[2])-1
        barcodes[(r, c)] = fields[0]           # note row, col
    return barcodes

def generate_tissue_positions_list(tissue_positions_list,
                                   spots,
                                   scale,
                                   mask,
                                   spot_diameter,
                                   tspots,
                                   barcodes):
    """
    Function for generating the tissue_positions_list output file. Determines which spots overlap
    the tissue section based on the spot positions and the tissue mask. The overlap is determined
    as the fraction of spot pixels that overlap the mask. This calculation uses the spot radius and
    image scale to estimate the ovelap.

        Args:
            tissue_positions_list (str):
            spots (iterable): A list of spots centre coordinates (row, col, y_pixel, x_pixel).
            scale (float): The float used to scale the input image (0.0 - 1.0).
            mask (numpy array): Binary mask defining tissue covered pixels in scaled image.
            spot_diameter (int): The estimated spot diameter in pixels.
            tspots (pandas DataFrame): Dataframe with row and column coords for spots under tissue.
            barcodes (dict): Dictionary with (row, column) as keys and barcode sequences as values.

        Returns:
            tissue_barcodes (list): List of barcode sequences overlapping the tissue section.
            img_tissue (set): Set with tuples of scaled spot coords (pixels) overlapping tissue.
            img_no_tissue (set): Set with tuples of scaled spot coords (pixels) outside tissue.
    """

    # generate tissue positions list
    tissue_barcodes = []
    img_tissue = set()
    img_no_tissue = set()
    missing=[]

    with open(tissue_positions_list, 'w') as fd:
        for row, col, rowpos, colpos in spots:
            srowpos = int(rowpos*scale)
            scolpos = int(colpos*scale)
            if mask is not None:
                height, width = mask.shape[:2]
                if 0 <= srowpos < height and 0 <= scolpos < width:
                    tissue = int(mask[srowpos, scolpos] > 0)
                else:
                    tissue = 0
#                # get the spot pixels and determine if at least 50% of them is tissue (1 in mask)
#                y_axis, x_axis = np.ogrid[:mask.shape[0], :mask.shape[1]]
#                spot_centre_dist = np.sqrt((x_axis-scolpos)**2+(y_axis-srowpos)**2)
#                spot_area_mask = spot_centre_dist <= spot_diameter/2.0
#                if np.mean(mask[spot_area_mask].astype('float64')) >= np.float64(0.5):
#                    tissue = 1
#                else:
#                    tissue = 0
            else:
                tissue = int(any((tspots['x'] == col)&(tspots['y'] == row).values))

            if (row, col) in barcodes:
                barcode = barcodes[(row, col)]
                barcode += "-1"         # CELLRANGER-1761 needs to address this
                if tissue > 0:
                    tissue_barcodes.append(barcode)
                    img_tissue.add((srowpos, scolpos))
                else:
                    img_no_tissue.add((srowpos, scolpos))
                fd.write("{},{},{},{},{},{}\n".format(barcode, tissue, row, col, rowpos, colpos))
            else:
                missing.append((row,col))

    # This is temporary below as I don't want to kill the pipeline, but this should be converted
    # to an error or, better, a preflight check.  This should never happen.
    if len(missing):
        martian.alarm("barcodes were missing from the coordinates file (up to ten printed): {}".format(missing[:10]))

    return tissue_barcodes, img_tissue, img_no_tissue

def generate_manual_tissue_positions_list(tissue_positions_list,
                                          spots,
                                          scale,
                                          tissue_binary_mask,
                                          barcodes,
                                          gemgroup):
    """
    Function for generating a tissue positions list file for downstream compatibility
    from a predetermined list of spots and a tissue binary mask.

        Args:
            tissue_positions_list (str): The path to write the tissue positions list to.
            spots (iterable): A list of spot centers (row, col, y_pixel, x_pixel)
            scale (float): The float used to scale the input image (0.0 - 1.0)
            tissue_binary_mask (iterable): A binary list of whether the spots are in tissue.
            barcodes (dict): Dicionary with (row, column) as keys and barcode sequences as values.
            gemgroup (int): The GEM group suffix to add to each barcode.

        Returns:
            tissue_barcodes (list): List of all barcode sequences overlapping the tissue section
            img_tissue (set): Set with tuples of scaled spot coords (pixels) overlapping tissue.
            img_no_tissue (set): Set with tuples of scaled spot coords (pixels) outside tissue.
    """
    # generate tissue positions list
    tissue_barcodes = []
    img_tissue = set()
    img_no_tissue = set()
    with open(tissue_positions_list, 'w') as fd:
        for ((row, col, rowpos, colpos), in_tissue) in zip(spots, tissue_binary_mask):
            srowpos = int(rowpos*scale)
            scolpos = int(colpos*scale)
            if (row, col) in barcodes:
                barcode = "%s-%s" % (barcodes[(row, col)], str(gemgroup))
                if in_tissue > 0:
                    tissue_barcodes.append(barcode)
                    img_tissue.add((srowpos, scolpos))
                else:
                    img_no_tissue.add((srowpos, scolpos))
                fd.write("{},{},{},{},{},{}\n".format(barcode, in_tissue, row, col, rowpos, colpos))

    return tissue_barcodes, img_tissue, img_no_tissue


def write_QC_image(qc, scale, img_tissue, img_no_tissue, spot_diameter, fiducial_coords,
                    fiducial_diameter, detected_tissue_image, text_annotation=None, spots_alpha=0.25):
    """
    Writes the QC image to disk. Spots overlapping with the tissue section in blue.

        Args:
            qc (numpy array): QC image from the get_mask function to plot the spots on
            img_tissue (iterable): list or similar of image coordinates for spots under tissue
            img_no_tissue (iterable): list or similar of image coordinates for spots outside tissue
            spot_diameter (integer): the diameter of a spot (in fullres pixels)
            fiducial_coords (iterable): list of image coordinates row,col for the fiducials OR None
            fiducial_diameter (integer): the diameter of a fiducial (in fullres pixels)
            detected_tissue_image (str): path to the detected_tissue_image output
            text_annotation (str): Annotation text to be added to QC image
            spots_alpha (float): Tissue spots are colored according to (1-spots_alpha)*background + spots_alpha*blue.

        Returns:
             None
    """

    # write the QC image

    height, width  = qc.shape[:2]
    fig, ax = plt.subplots(figsize=(width/72.0,height/72.0),dpi=72)
    ax.imshow(qc)
    if spot_diameter == 1:
        spot_diameter = matplotlib.rcParams['lines.markersize']/scale
    if fiducial_diameter == 1:
        fiducial_diameter = matplotlib.rcParams['lines.markersize']/scale

    ax.scatter(np.array([int(_c) for _r, _c in img_tissue]),
               np.array([int(_r) for _r, _c in img_tissue]),
               c='blue', s=int((spot_diameter*scale)**2), edgecolor='black', alpha=spots_alpha,
               linewidth=3)
    ax.scatter(np.array([int(_c) for _r, _c in img_no_tissue]),
               np.array([int(_r) for _r, _c in img_no_tissue]),
               c='None', s=int((spot_diameter*scale)**2), edgecolor='black', alpha=0.25,
               linewidth=3)
    if fiducial_coords is not None:
        ax.scatter(
               np.array([int(_c) for _, _c in fiducial_coords]),
               np.array([int(_r) for _r, _ in fiducial_coords]),
               c='None', s=int((fiducial_diameter*scale)**2), edgecolor='red', alpha=0.75,
                    linewidth=3)
    ax.set_axis_off()
    if text_annotation:
        ax.text(0.98, 0.02, text_annotation, horizontalalignment='right',
                verticalalignment='center', transform=ax.transAxes, fontsize=36)
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    fig.savefig(detected_tissue_image, format='jpg', dpi=72, pad_inches=0)
    plt.close(fig)

def estimate_spot_diameter(spots):
    """
    Function that takes a list of spots centre coordinates (row, col, y_pixel, x_pixel).
    The spot radius (in pixels) is estimated as a fixed proportion of the average spot to spot
    distance between spots with even column coordinates.
    This estimate will only be valid for 5k arrays.
    If input spot coordinates is from a 1k array (both even and odd cordinates on same row) the
    function will return an estimated spot radius of 1px (falling back to spot centers).

        Args:
            spots (iterable): A list of spots centre coordinates (row, col, y_pixel, x_pixel).

        Returns:
             spot_diameter (int): The estimated spot radius in full-resolution pixels.
    """

    # make a dictionary to sort spots by row and column coordinates
    positions_by_row_col = {}
    for spot in spots:
        row, col, rowpos, colpos = spot
        try:
            positions_by_row_col[row][col] = (colpos, rowpos)
        except KeyError:
            positions_by_row_col[row] = {col: (colpos, rowpos)}

    # calculate center to center distances of spots, for each spot calculate the distance to the
    # neighbour straight to the right i.e. the spot on the same row 2 columns later
    distances = []
    for row,cols in positions_by_row_col.items():
        for col in cols:
                try:
                    x1, y1 = positions_by_row_col[row][col]
                    x2, y2 = positions_by_row_col[row][col+2]
                    distances.append(((x2-x1)**2+(y2-y1)**2)**(1/2.0))
                except KeyError:
                    pass # End of row or missing spots

                # if both even and odd columns are present on the same row this is a old array
                # i.e fall back to no spot areas (estimated radius = 1 pixel)
                if col+1 in positions_by_row_col[row]:
                    return 1

    average_spot_to_spot_distance = np.mean(np.array(distances))

    spot_radius = average_spot_to_spot_distance*0.65*0.5

    if spot_radius > 1:
        return 2.0*spot_radius
    else:
        return 1.0
