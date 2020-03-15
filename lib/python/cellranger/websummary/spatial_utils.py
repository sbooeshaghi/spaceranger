#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
""" Utilities designed to help create the web summary used in Spatial."""


from collections import namedtuple
import json

import pandas as pd
import cellranger.metrics_names as metrics_names

import cellranger.utils as cr_utils
import cellranger.report as cr_report #pylint: disable=no-name-in-module
import cellranger.websummary.sample_properties as samp_props
from cellranger.spatial.image import WebImage
import cellranger.webshim.common as cr_webshim
from cellranger.websummary.metrics import SpatialMetricAnnotations, output_metrics_csv_from_annotations
from cellranger.websummary.react_components import SharedCoordinatePlotCollection

import cellranger.websummary.plotly_tools as pltly
import cellranger.websummary.analysis_tab as cr_at
import cellranger.websummary.summary_tab as cr_st
import cellranger.websummary.web_summary_builder as ws_builder

SPATIAL_PIPELINE_NAME = "count"
SPATIAL_COMMAND_NAME = "Space Ranger"

SpatialReportArgs = namedtuple('SpatialReportArgs',
                                ['sample_id', 'sample_desc', 'scalefactors', 'matrix',
                                'tissue_lowres_image', 'detected_tissue_image',
                                'tissue_positions_list', 'analysis'])


SENSITIVITY_PLOT_HELP = {
    "data": [["",
              ["(left) Total UMI counts for each spot overlayed on the tissue image. "
               "Spots with greater UMI counts likely have higher RNA content than spots "
                "with fewer UMI counts. ",
               "(right) Total UMI counts for spots displayed by a 2-dimensional embedding produced "
               "by the t-SNE algorithm. In this space, pairs of spots that are close to each other "
               "have more similar gene expression profiles than spots that are distant from each "
               "other."
               ]]],
    "title": "UMIs Detected"
}


def _make_array_plot_layout(title, encoded_image):
    xlo, xhi, ylo, yhi = encoded_image.cropbox

    layout = {
        "title": {
            "text": title,
            "yanchor": "top",
        },
        "autosize": True,
        "hovermode": "closest",
        "legend": {"y": 0.9, "yanchor": "top"},
        "xaxis": {
            "type": "linear",
            "fixedrange": True,
            "showgrid": False,
            "zeroline": False,
            "showticklabels": False,
            "range": [xlo, xhi],
            "scaleanchor": "y",
            "scaleratio" : 1.0,
        },
        "yaxis": {
            "type": "linear",
            "fixedrange": True,
            "showgrid": False,
            "zeroline": False,
            "showticklabels": False,
            "range": [yhi, ylo],
            "domain": [0.0,0.90],
        },
        "margin": {
            "b": 10,
            "t": 10,
        },
        "images": [
            {
                "source": encoded_image.base64_encoded_str,
                "xref": "x",
                "yref": "y",
                "x": 0,
                "y": 0,
                "sizex": encoded_image.width,
                "sizey": encoded_image.height,
                "sizing": "stretch",
                "opacity": 1,
                "layer": "below",
            },
        ],
        "sliders": [
            {
                "pad": {"t": 10,
                        "r": 0,
                        "b": 10},
                "yanchor": 'center',
                "xanchor": 'left',
                "active": 4,
                "showactive": True,
                "tickwidth": 1,
                "currentvalue": {
                    "xanchor": 'center',
                    "prefix": 'Spot Opacity: ',
                    "font": {
                        "color": '#888',
                        "size": 12
                    },
                    "visible": False,
                },
                "steps": [{
                    "label": '0.0',
                    "method": 'restyle',
                    "args": ['opacity', 0.00]
                }, {
                    "label": '',
                    "method": 'restyle',
                    "args": ['opacity', 0.25]
                }, {
                    "label": '',
                    "method": 'restyle',
                    "args": ['opacity', 0.50]
                }, {
                    "label": '',
                    "method": 'restyle',
                    "args": ['opacity', 0.75]
                }, {
                    "label": '1.0',
                    "method": 'restyle',
                    "args": ['opacity', 1.00]
                }]
            }
        ],
    }
    return layout

def make_array_plot(title, encoded_image, data):
    """Generate plot layout for displaying spatial images overlaid with spots"""
    assert isinstance(encoded_image, WebImage)
    return {
        "config": pltly.SPATIAL_PLOT_CONFIG,
        "layout": _make_array_plot_layout(title, encoded_image),
        "data": data
    }

def umi_on_spatial_plot(sample_data, coords, encoded_image):
    """

    :param sample_data: Instance of webshim.data.SampleData
    :param encoded_image: Instance of WebImage
    :return:
    """
    analysis = sample_data.get_analysis(cr_at.SingleGenomeAnalysis)
    if analysis is None:
        return None

    # make a df with the barcodes as the index
    bc_umi_counts = pd.DataFrame(index=analysis.matrix.bcs)
    # Add the UMI column by getting the sum UMI for each barcode
    reads_per_bc = analysis.matrix.get_counts_per_bc()
    bc_umi_counts['UMI'] = reads_per_bc
    cbd = coords.join(bc_umi_counts)
    cbd = cbd.dropna()
    color = cr_at.get_umi_color(cbd['UMI'].values)
    title = "Tissue Plot with Spots Colored by UMI Count"
    data = [{
        "name": "Spots",
        "x": list(cbd["imagex"].values),
        "y": list(cbd["imagey"].values),
        "type": "scatter",
        "mode": "markers",
        "marker": {
            "opacity": 0.9,
            "size": encoded_image.markersize,
            "sizemode": "diameter",
            "color": color,
            "colorscale": "Jet",
            "colorbar": {
                "title": "UMI counts",
                "yanchor": "top",
                "y": 0.9,
                "len": 0.9,
                "xpad": 20,
         },
        },
        "text": ['UMI counts: {:,d}'.format(reads) for reads in reads_per_bc],
    }]
    return make_array_plot(title, encoded_image, data)

def tissue_plots_by_clustering_spatial(sample_data, coords, encoded_image):
    """
    Get the data structure that represents the plots for all the tissues.

    :param sample_data: Instance of webshim.data.SampleData
    :param coords:
    :param encoded_image: Instance of WebImage
    :return: A SharedCoordinatePlotCollection object
    """
    if sample_data is None:
        return None

    analysis = sample_data.get_analysis(cr_at.SingleGenomeAnalysis)
    if analysis is None:
        return None

   # order clustering by order: graph, kmeans=2, =3 etc
    clusterings = cr_at.sort_analysis_clusterings(analysis.clusterings)
    # align clustering label with image by barcode
    clusters_dict = {clustering_key:clustering.clusters for
                     (clustering_key, clustering) in clusterings}
    clusters_dict['barcode'] = analysis.matrix.bcs
    clusters_df = pd.DataFrame.from_dict(clusters_dict)
    clusters_df.set_index('barcode', inplace=True)
    cluster_tissue_df = clusters_df.join(coords[['imagex', 'imagey']])

    title = "Tissue Plot with Spots Colored by Clustering"
    # pylint: disable=too-many-function-args
    tissue_plot_data = SharedCoordinatePlotCollection(list(cluster_tissue_df['imagex']),
                                                     list(cluster_tissue_df['imagey']),
                                                     pltly.SPATIAL_PLOT_CONFIG,
                                                     _make_array_plot_layout(title, encoded_image),
                                                     "scatter",
                                                     {"size": encoded_image.markersize, "sizemode": "diameter"})
    for (clustering_key, clustering) in clusterings:
        tissue_plot_data.add_clustering(clustering_key, clustering)
    return tissue_plot_data

def get_scalefactors(scalefactors_fn):
    return json.load(open(scalefactors_fn, "r"))


def get_scaled_coordinates(tissue_positions_list_fn, scalefactors_fn):
    coords = pd.read_table(tissue_positions_list_fn,
                           names=["barcode", "tissue", "y", "x", "imagey", "imagex"],
                           index_col="barcode", sep=",")

    # read in scalefactors json and adjust coords for downsampled image
    scalef = get_scalefactors(scalefactors_fn)
    coords["imagey"] *= scalef["tissue_lowres_scalef"]
    coords["imagex"] *= scalef["tissue_lowres_scalef"]
    return coords


def create_common_spatial_summaries(args, outs, metrics_csv_out=None):
    """This method generates the CS portion of the summarize stage outputs
    it is shared by both the PD and CS code."""

    # Both manual and automatic alignment pathways should generate the
    # fraction of spots under tissue metric
    assert args.fraction_under_tissue is not None
    alignment_dict = {
        'aligned_fiducials': args.aligned_fiducials,
        'fraction_under_tissue': args.fraction_under_tissue,
    }
    id_dict = {
        "sample_id": args.sample_id,
        "sample_desc": args.sample_desc,
        "spatial_slide_info": args.slide_serial_info,
    }
    cr_report.merge_jsons(args.summaries, outs.metrics_summary_json, [alignment_dict, id_dict])
    ref_genomes = cr_utils.get_reference_genomes(args.reference_path)
    sample_properties = samp_props.ExtendedCountSampleProperties(sample_id=args.sample_id,
                                                                 sample_desc=args.sample_desc,
                                                                 barcode_whitelist=args.barcode_whitelist,
                                                                 reference_path=args.reference_path,
                                                                 genomes=ref_genomes)
    sample_data_paths = samp_props.SampleDataPaths(
        summary_path=outs.metrics_summary_json,
        barcode_summary_path=args.barcode_summary_h5,
        analysis_path=args.analysis,
        filtered_barcodes_path=args.filtered_barcodes
    )
    sample_data = cr_webshim.load_sample_data(sample_properties, sample_data_paths)

    spatial_args = SpatialReportArgs(sample_id=args.sample_id,
                                     sample_desc=args.sample_desc,
                                     tissue_lowres_image=args.tissue_lowres_image,
                                     scalefactors=args.scalefactors,
                                     matrix=args.matrix,
                                     tissue_positions_list=args.tissue_positions_list,
                                     detected_tissue_image=args.detected_tissue_image,
                                     analysis=args.analysis)
    web_sum_data = build_web_summary_data_spatial(sample_properties, sample_data, spatial_args)
    ## Make a metrics.csv file
    if metrics_csv_out: # not made in PD
        spatial_mets = SpatialMetricAnnotations()
        output_metrics_csv_from_annotations(spatial_mets, sample_data, metrics_csv_out, sample_properties.genomes)
    return web_sum_data

def _add_image_alignment_alarm(web_sum_data, sample_data, metadata):
    """
    Adds a warning if issues detected with image registration
    :param web_sum_data:
    :param sample_data:
    :param metadata: A SpatialMetricsAnnotation instance
    :return:
    """
    assert isinstance(metadata, SpatialMetricAnnotations)
    alarms = metadata.gen_metric_list(sample_data.summary, [metrics_names.SUSPECT_ALIGNMENT], [])
    alarm_dicts = [metric.alarm_dict for metric in alarms if metric.alarm_dict]
    if alarm_dicts:
        web_sum_data.alarms.extend(alarm_dicts)

def build_web_summary_data_spatial(sample_properties, sample_data, spatial_args):
    metadata = SpatialMetricAnnotations()
    species_list = sample_properties.genomes
    web_sum_data = ws_builder.build_web_summary_data_common(sample_properties,
                                                 sample_data, SPATIAL_PIPELINE_NAME, metadata,
                                                 SPATIAL_COMMAND_NAME)

    _add_image_alignment_alarm(web_sum_data, sample_data, metadata)

    detected_tissue_image = WebImage(spatial_args.detected_tissue_image)
    small_img = detected_tissue_image.resize_and_encode_image(new_width=470)
    zoom_images = {
        "small_image" : small_img.base64_encoded_str,
        "big_image" :detected_tissue_image.base64_encoded_str,
        "sizes" : {
            "width": detected_tissue_image.width,
            "height": detected_tissue_image.height
        }
    }
    cr_st.add_data(web_sum_data.summary_tab, web_sum_data.alarms,
             cr_st.spot_calling_table(metadata, sample_data, species_list, zoom_images,
                                cr_st.CELL_CALLING_METRIC_KEYS, cr_st.CELL_CALLING_ALARM_KEYS))


    # Add in the clustering over tissue plots to the left pane of the clustering selector
    scalef = get_scalefactors(spatial_args.scalefactors)
    coords = get_scaled_coordinates(spatial_args.tissue_positions_list, spatial_args.scalefactors)
    l,r=min(coords["imagex"]), max(coords["imagex"])
    t,b=min(coords["imagey"]), max(coords["imagey"])
    hoffset = 0.13 * (r - l + 1) # arbitrary looking numbers due to difference
    voffset = 0.10 * (b - t + 1) # in vertical vs. horizontal packing

    # ensure that the crop box is within the image
    # (plotly does weird things otherwise)
    cb_x0 = max(l-hoffset, 0)
    cb_x1 = min(r+hoffset, detected_tissue_image.width-1)
    cb_y0 = max(t-voffset, 0)
    cb_y1 = min(b+voffset, detected_tissue_image.height-1)

    # compute marker size for this image
    plot_height = 265.5     # this is fixed for our design, but in a better world we'd
                            # get this from Plotly.  If this works, we can consider
                            # just updating markerref in javascript instead of doing this
    spot_diameter_image = scalef["tissue_lowres_scalef"]*scalef["spot_diameter_fullres"]
    scaled_spot_diameter = spot_diameter_image * plot_height / ( cb_y1 - cb_y0 )

    lowres_tissue = WebImage(spatial_args.tissue_lowres_image, \
                                    cropbox=[cb_x0, cb_x1, cb_y0, cb_y1], \
                                    markersize=scaled_spot_diameter)

    if web_sum_data.clustering_selector:
        web_sum_data.clustering_selector.left_plots = \
            tissue_plots_by_clustering_spatial(sample_data,
                                               coords,
                                               lowres_tissue)

    # Make UMI Count plots
    umi_spatial_plot = umi_on_spatial_plot(sample_data, coords, lowres_tissue)
    umi_tsne_plot = cr_at.umi_on_tsne_plot(sample_data, spatial=True)
    if umi_spatial_plot and umi_tsne_plot:
        web_sum_data.analysis_tab['umi_plots'] = {
            "help_txt": SENSITIVITY_PLOT_HELP,
            "umi_spatial_plot": umi_spatial_plot,
            cr_at.UMI_TSNE_PLOT: umi_tsne_plot
        }
    return web_sum_data
