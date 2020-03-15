import os
import json
import pandas as pd
import cv2
import cellranger.spatial.manual_alignment as ma_utils
import cellranger.spatial.utils as spatial_utils
import cellranger.spatial.tissue_detection as tissue_detection
import martian

__MRO__ = """
stage DETECT_TISSUE(
    in  jpg    downsampled_tissue_image,
    in  json   scalefactors_json,
    in  path   spot_data_sel,
    in  bool   skip_tissue_detection,
    in  json   gpr_data,
    in  path   loupe_alignment_file,
    in  txt    spot_positions_list,
    int txt    fiducial_positions_list,
    in  string barcode_whitelist,
    out jpg    detected_tissue_image,
    out csv    tissue_positions_list,
    out json   barcodes_under_tissue,
    out float  fraction_under_tissue,
    src py     "stages/spatial/detect_tissue",
)
"""

def main(args, outs):
    if args.skip_tissue_detection:
        return

    if args.loupe_alignment_file:
        fraction_under_tissue = generate_manual_tissue_files(args.loupe_alignment_file,
                                     args.downsampled_tissue_image,
                                     args.scalefactors_json,
                                     args.barcode_whitelist,
                                     outs.detected_tissue_image,
                                     outs.tissue_positions_list,
                                     outs.barcodes_under_tissue)
    else:
        fraction_under_tissue = detect_tissue(args.downsampled_tissue_image,
                                              args.scalefactors_json,
                                              args.spot_positions_list,
                                              args.fiducial_positions_list,
                                              args.barcode_whitelist,
                                              args.spot_data_sel,
                                              args.gpr_data,
                                              outs.detected_tissue_image,
                                              outs.tissue_positions_list,
                                              outs.barcodes_under_tissue)
    outs.fraction_under_tissue = fraction_under_tissue

def detect_tissue(tissue_hires_image,
                  scalefactors_json,
                  spot_positions_list,
                  fiducial_positions_list,
                  barcode_whitelist,
                  spot_data_sel,
                  gpr_data,
                  detected_tissue_image,
                  tissue_positions_list,
                  barcodes_under_tissue):

    #
    # a bunch of this should move into helper
    #
    ### Read the downsampled tissue image
    rbf = spatial_utils.cv_read_image_standard(tissue_hires_image)
    scalefactors = json.load(open(scalefactors_json, 'r'))
    scale = scalefactors["tissue_hires_scalef"]

    # read spot positions list
    if spot_positions_list is not None:
        spots = [tuple(map(int, x.strip().split(','))) for x in open(spot_positions_list, 'r')]
    else:
        raise RuntimeError('No spots found among input parameters to tissue detection.')

    if fiducial_positions_list is not None and os.path.exists(fiducial_positions_list):
        fiducials =  [tuple(map(int, x.strip().split(','))) for x in open(fiducial_positions_list, 'r')]
        fiducials =  [ (int(scale*r+0.5),int(scale*c+0.5)) for _,_,r,c in fiducials ]
    else:
        fiducials = None

    # get the spot radius from 5K fiducial array or use a sentinel value
    spot_diameter=scalefactors.get("spot_diameter_fullres",1.0)
    fiducial_diameter=scalefactors.get("fiducial_diameter_fullres", 1.0)

    if spot_data_sel is None:
        bounding_box = tissue_detection.get_bounding_box(spots, scale, spot_diameter)
        mask, qc = tissue_detection.get_mask(rbf,
                                             plot=True,
                                             plot_title=None,
                                             plotfname='initialization.jpg',
                                             bounding_box=bounding_box)
        tspots = None

    else:
        # read selected spot positions list
        tspots = pd.read_csv(spot_data_sel, sep='\t')
        tspots['y'] -= 1
        tspots['x'] -= 1
        mask = None
        qc = cv2.cvtColor(rbf.astype('uint8'), cv2.COLOR_GRAY2RGB)

    barcodes = tissue_detection.read_barcode_coordinates(barcode_whitelist)

    tissue_barcodes, img_tissue, img_no_tissue = tissue_detection.generate_tissue_positions_list(
                                                                            tissue_positions_list,
                                                                            spots,
                                                                            scale,
                                                                            mask,
                                                                            spot_diameter,
                                                                            tspots,
                                                                            barcodes)
    fraction_under_tissue = len(tissue_barcodes) / float(len(barcodes))
    if fraction_under_tissue == 0:
        martian.exit("No tissue is detected on the spots by automatic alignment. Please use manual alignment.")

    # write json
    json.dump(tissue_barcodes, open(barcodes_under_tissue, 'w'))

    tissue_detection.write_QC_image(qc, scale,
                                    img_tissue,
                                    img_no_tissue,
                                    spot_diameter,
                                    fiducials,
                                    fiducial_diameter,
                                    detected_tissue_image)
    return fraction_under_tissue

def generate_manual_tissue_files(loupe_alignment_file,
                                 tissue_hires_image,
                                 scalefactors_json,
                                 barcode_whitelist,
                                 detected_tissue_image,
                                 tissue_positions_list,
                                 barcodes_under_tissue):
    # use alignment file + tissue image to label spots
    rbf = spatial_utils.cv_read_image_standard(tissue_hires_image)
    scalefactors = json.load(open(scalefactors_json, 'r'))
    scale = scalefactors["tissue_hires_scalef"]

    # read spot locations from loupe alignment file
    if loupe_alignment_file is None:
        raise RuntimeError('No manual alignment file supplied.')

    with open(loupe_alignment_file, 'r') as infile:
        spot_data = json.load(infile)
    spots = ma_utils.read_spot_positions_list_from_alignment(spot_data)

    # read also fiducial locations
    loupe_fiducials = spot_data['fiducial']
    fiducials = list()
    for fiducial_dict in loupe_fiducials:
        fiducials.append((int(scale * fiducial_dict['imageY']), int(scale * fiducial_dict['imageX'])))
    fiducial_diameter = scalefactors.get("fiducial_diameter_fullres", 1.0)

    # Estimate the spot radius
    spot_diameter = tissue_detection.estimate_spot_diameter(spots)

    # have this in our back pocket if we need to show bounding box on manual alignment
    # (may need color conversion to run first) - TALK TO NEIL FIRST - this has
    # changed as the previous bounding box assumed axis aligned images
    #bounding_box = tissue_detection.get_bounding_box(spots, scale, spot_radius)
    #qc = tissue_detection.get_manual_qc_image(rbf, plot=True, bounding_box=bounding_box)

    # this is eventually done somewhere in get_mask but we need it here explicitly
    qc = cv2.cvtColor(rbf.astype('uint8'), cv2.COLOR_GRAY2RGB)
    barcodes = tissue_detection.read_barcode_coordinates(barcode_whitelist)
    tissue_mask = ma_utils.read_spot_tissue_mask_from_alignment(spot_data)

    tissue_barcodes, img_tissue, img_no_tissue = \
        tissue_detection.generate_manual_tissue_positions_list(
            tissue_positions_list,
            spots,
            scale,
            tissue_mask,
            barcodes,
            1) # gemgroup 1 -- see CELLRANGER-1761 to be more flexible
    fraction_under_tissue = len(tissue_barcodes) / float(len(barcodes))

    json.dump(tissue_barcodes, open(barcodes_under_tissue, 'w'))

    # There isn't a tissue mask in the manual alignment path; only tissue spots are drawn in blue, therefore the
    # coloring needs adjusted so both paths would produce images with similar shades of blue.
    # The automatic path blends background and blue of tissue pixels in this way:
    # (1-tissue_alpha)*background + (tissue_alpha)*blue, where tissue_alpha is 0.25, see utils.cv_composite_labelmap;
    # then write_QC_image blends tissue dots in a similar way, but using spots_alpha. Simple arithmetic shows the
    # following transformation renders the tissue spots in similar shades of blue for both paths.
    tissue_alpha = 0.25
    spots_alpha = 0.25
    manual_align_alpha = spots_alpha + tissue_alpha - spots_alpha * tissue_alpha
    tissue_detection.write_QC_image(qc, scale,
                                    img_tissue,
                                    img_no_tissue,
                                    spot_diameter,
                                    fiducials,
                                    fiducial_diameter,
                                    detected_tissue_image,
                                    None,
                                    manual_align_alpha)

    return fraction_under_tissue
