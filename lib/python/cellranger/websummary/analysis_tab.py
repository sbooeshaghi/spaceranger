#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

# Here, you could find everything shown in the summary tab

import numpy as np
from copy import deepcopy
import cellranger.webshim.common as cr_webshim
import cellranger.webshim.constants.gex as ws_gex_constants
from cellranger.analysis.singlegenome import SingleGenomeAnalysis
import cellranger.websummary.plotly_tools as pltly
from cellranger.websummary.react_components import (Clusterings,
                                                    SharedCoordinatePlotCollection,
                                                    ClusteringData,
                                                    ClusteringSelector)

TSNE_LAYOUT_CONFIG = {
    "xaxis": {
        "type": "linear",
        "title": "t-SNE1",
        "showline": False,
        "zeroline": True,
        "fixedrange": False,
    },
    "yaxis": {
        "type": "linear",
        "title": "t-SNE2",
        "showline": False,
        "zeroline": True,
        "fixedrange": False,
    },
    "margin": {'t': 30},    # needed to keep screenshots from hitting the top
    'hovermode': 'closest',
}

def _make_tsne_umi_layout(title):
    new_layout = deepcopy(TSNE_LAYOUT_CONFIG)
    new_layout["title"] = title
    new_layout["hover_mode"] = "closest"
    return new_layout

DIFFEXP_TABLE_HELP = {
    "helpText": "The differential expression analysis seeks to find, for each cluster, features that are more highly expressed in that cluster relative to the rest of the sample. "
    "Here a differential expression test was performed between each cluster and the rest of the sample for each feature. "
    "The Log2 fold-change (L2FC) is an estimate of the log2 ratio of expression in a cluster to that in all other cells. "
    "A value of 1.0 indicates 2-fold greater expression in the cluster of interest. "
    "The p-value is a measure of the statistical significance of the expression difference and is based on a negative binomial test. "
    "The p-value reported here has been adjusted for multiple testing via the Benjamini-Hochberg procedure. "
    "In this table you can click on a column to sort by that value. "
    "Also, in this table features were filtered by (Mean UMI counts > 1.0) and the top N features by L2FC for each cluster were retained. "
    "Features with L2FC < 0 or adjusted p-value >= 0.10 were grayed out. "
    "The number of top features shown per cluster, N, is set to limit the number of table entries shown to 10,000; N=%10,000/K^2 where K is the number of clusters. "
    "N can range from 1 to 50. "
    "For the full table, please refer to the 'differential_expression.csv' files produced by the pipeline.",
    "title": "Top Features by Cluster (Log2 fold-change, p-value)"
}

SPATIAL_DIFFEXP_TABLE_HELP = {
    "helpText": "The differential expression analysis seeks to find, for each cluster, features that are more highly expressed in that cluster relative to the rest of the sample. "
    "Here a differential expression test was performed between each cluster and the rest of the sample for each feature. "
    "The Log2 fold-change (L2FC) is an estimate of the log2 ratio of expression in a cluster to that in all other spots. "
    "A value of 1.0 indicates 2-fold greater expression in the cluster of interest. "
    "The p-value is a measure of the statistical significance of the expression difference and is based on a negative binomial test. "
    "The p-value reported here has been adjusted for multiple testing via the Benjamini-Hochberg procedure. "
    "In this table you can click on a column to sort by that value. "
    "Also, in this table features were filtered by (Mean UMI counts > 1.0) and the top N features by L2FC for each cluster were retained. "
    "Features with L2FC < 0 or adjusted p-value >= 0.10 were grayed out. "
    "The number of top features shown per cluster, N, is set to limit the number of table entries shown to 10,000; N=%10,000/K^2 where K is the number of clusters. "
    "N can range from 1 to 50. "
    "For the full table, please refer to the 'differential_expression.csv' files produced by the pipeline.",
    "title": "Top Features by Cluster (Log2 fold-change, p-value)"
}

TSNE_CLUSTERING_PLOT_HELP = {
    "data": [["",
              ["(left) Shown here are the total UMI counts for each cell-barcode. "
               "Cells with greater UMI counts likely have higher RNA content than cells with fewer UMI counts. "
               "The axes correspond to the 2-dimensional embedding produced by the t-SNE algorithm. "
               "In this space, pairs of cells that are close to each other have more similar gene expression profiles than cells that are distant from each other. "
               "The display is limited to a random subset of cells.",
               "(right) These are the assignments of each cell-barcode a clusters by an automated clustering algorithm. "
               "The clustering groups together cells that have similar expression profiles. "
               "The axes correspond to the 2-dimensional embedding produced by the t-SNE algorithm. "
               "In this space, pairs of cells that are close to each other have more similar gene expression profiles than cells that are distant from each other. "
               "The display is limited to a random subset of cells. "
               "Please use Loupe Browser to view the entire dataset."]]],
    "title": "t-SNE Projection"
}

SPATIAL_TSNE_CLUSTERING_PLOT_HELP = {
    "data": [["",
              ["(left) Shown here are the total UMI counts for each spot-barcode. "
               "Spots with greater UMI counts likely have higher RNA content than spots with fewer UMI counts. "
               "The axes correspond to the 2-dimensional embedding produced by the t-SNE algorithm. "
               "In this space, pairs of spots that are close to each other have more similar gene expression profiles than spots that are distant from each other. "
               "The display is limited to a random subset of spots.",
               "(right) These are the assignments of each spot-barcode a clusters by an automated clustering algorithm. "
               "The clustering groups together spots that have similar expression profiles. "
               "The axes correspond to the 2-dimensional embedding produced by the t-SNE algorithm. "
               "In this space, pairs of spots that are close to each other have more similar gene expression profiles than spots that are distant from each other. "
               "The display is limited to a random subset of spots. "
               "Please use Loupe Browser to view the entire dataset."]]],
    "title": "Clustering"
}

SEQ_SATURATION_PLOT_HELP = {
    "helpText" : "This plot shows the Sequencing Saturation metric as a function of downsampled sequencing depth (measured in mean reads per cell), up to the observed sequencing depth. "
    "Sequencing Saturation is a measure of the observed library complexity, and approaches 1.0 (100%) when all converted mRNA transcripts have been sequenced. "
    "The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point. "
    "The dotted line is drawn at a value reasonably approximating the saturation point.",
    "title" : "Sequencing Saturation"
}

SPATIAL_SEQ_SATURATION_PLOT_HELP = {
    "helpText" : "This plot shows the Sequencing Saturation metric as a function of downsampled sequencing depth (measured in mean reads per spot), up to the observed sequencing depth. "
    "Sequencing Saturation is a measure of the observed library complexity, and approaches 1.0 (100%) when all converted mRNA transcripts have been sequenced. "
    "The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point. "
    "The dotted line is drawn at a value reasonably approximating the saturation point.",
    "title" : "Sequencing Saturation"

}

MEDIAN_GENE_PLOT_HELP = {
    "helpText": "This plot shows the Median Genes per Cell as a function of downsampled sequencing depth in mean reads per cell, up to the observed sequencing depth. "
    "The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point.",
    "title": "Median Genes per Cell"
}

SPATIAL_MEDIAN_GENE_PLOT_HELP = {
    "helpText": "This plot shows the Median Genes per Spot as a function of downsampled sequencing depth in mean reads per spot, up to the observed sequencing depth. "
    "The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point.",
    "title": "Median Genes per Spot"
}

BARNYARD_PLOT_HELP = [["Cell UMI Counts Plot",
                       ["Each point represents a cell-barcode. "
                        "The axes measure the total UMI counts in each cell-barcode that mapped to each transcriptome reference. "
                        "The points are colored by the number of inferred cells in the GEM associated with each barcode. "
                        "A multiplet represents either a GEM inferred to have encapsulated >1 cell or a barcode sequence that was shared by multiple single-cell GEMs."]]]

UMI_TSNE_PLOT = "umi_tsne_plot"



def barnyard_table(metadata, sample_data, sample_properties, species_list):
    """ Barnyard table and barnyard plot  """
    def barnyard_plot(sample_data, sample_properties):
        chart = {
            'config': pltly.PLOT_CONFIG,
            'layout': {
                'title': 'Cell UMI Counts',
                'showlegend': True,
                'hovermode': 'closest',
            },
            'data': [
                {
                    'x': [],
                    'y': [],
                    'mode': 'markers',
                    'type': 'scattergl',
                },
            ],
        }

        return cr_webshim.plot_barnyard_barcode_counts(chart, sample_properties, sample_data)

    if len(species_list) <= 1 or sample_data is None or sample_data.summary is None:
        return None

    metric_keys = [
        'filtered_bcs_observed_all',
        'filtered_bcs_inferred_multiplets',
        'filtered_bcs_inferred_multiplet_rate',
        'filtered_bcs_inferred_multiplet_rate_lb',
        'filtered_bcs_inferred_multiplet_rate_ub',
        'multi_filtered_bcs_mean_count_purity'
    ]

    gems = {}
    metrics = metadata.gen_metric_list(
        sample_data.summary, metric_keys, species_list)
    gems['table'] = {"rows": [[metric.name, metric.value_string]
                              for metric in metrics]}

    helptext = metadata.gen_metric_helptext(metric_keys) + BARNYARD_PLOT_HELP
    gems['help'] = {'title': 'GEM Partitions', 'data': helptext}

    barnyard = {
        'plot': barnyard_plot(sample_data, sample_properties),
        'gems': gems,
    }
    return {"barnyard": barnyard}

def cluster_on_tsne(clustering):
    clustering_labels = clustering.clusters
    num_cells = len(clustering_labels)
    data = []

    for i in range(max(clustering_labels)):
        index = i + 1
        name = "Cluster {}".format(index)
        indices = np.where(clustering_labels == index)[0]
        prop = len(indices)*1.0 / num_cells
        data.append({
            "name": "{} ({:.1%})".format(name, prop),
            "indices": list(indices),
            "type": "scattergl",
            "mode": "markers",
            "marker": {"opacity": 0.9, "size": 4},
            "text": "{}: {:.1%}".format(name, prop)
        })
    return data

def diffexp_table(diffexp, clustering, analysis):
    """ Show diffexp table for the given clustering """
    n_clusters = clustering.clusters.max()

    # Limit the number of entries in the DE table
    n_genes = int(
        np.floor(float(ws_gex_constants.MAX_DE_TABLE_ENTRIES) / (n_clusters**2)))
    if n_genes < 1:
        n_genes = 1
    elif n_genes > ws_gex_constants.MAX_TOP_N_GENES:
        n_genes = ws_gex_constants.MAX_TOP_N_GENES

    columns = [{
        "Header": "Feature",
        "columns": [{"Header": "ID", "accessor": "feature.id"},
                    {"Header": "Name", "accessor": "feature.fn"}]
    }]

    # Get the union of top DE genes
    top_genes = set()
    for i in xrange(n_clusters):
        # Filter genes by mean count and sort by log2 fold-change, descending
        means = diffexp.data[:, 0+3*i]
        log2fcs = diffexp.data[:, 1+3*i]

        keep_indices = np.flatnonzero(
            means >= ws_gex_constants.TOP_DE_GENES_MIN_MEAN)
        top_gene_indices = keep_indices[log2fcs[keep_indices].argsort()[
            ::-1]][:n_genes]

        for j in top_gene_indices:
            top_genes.add(analysis.matrix.int_to_feature_id(j))

        columns.append({"Header": "Cluster {}".format(i+1),
                        "columns": [{
                            "Header": "L2FC",
                            "accessor": "c{}.l".format(i+1),
                            "greyedout": "c{}.g".format(i+1)},
            {
                            "Header": "p-value",
                            "accessor": "c{}.p".format(i+1),
                            "greyedout": "c{}.g".format(i+1)}
        ]})

    table = []
    for gene_id in top_genes:
        i = analysis.matrix.feature_id_to_int(gene_id)
        gene_name = analysis.matrix.feature_id_to_name(gene_id)

        row = {"feature": {"fn": gene_name, "id": gene_id}}

        for j in xrange(n_clusters):
            log2fc = diffexp.data[i, 1+(3*j)]
            adj_p_value = diffexp.data[i, 2+(3*j)]

            if log2fc <= 0 or adj_p_value >= ws_gex_constants.PVALUE_DEEMPHASIS_CUTOFF:
                greyed = True
            else:
                greyed = False

            cn = "c{}".format(j+1)
            row[cn] = {
                "p": adj_p_value,
                "l": log2fc,
            }
            # Only output this if True
            # And set it as 1 to use less space than the boolean in JSON
            # And be interpreted correctly by the web summary code (copied below)
            # const gcval = get(row.original, gc);
            # const col = gcval ? '#DDD' : default_col;
            if greyed:
                row[cn]["g"] = 1

        table.append(row)

    # Sort by log2fc, descending, in first cluster
    if n_clusters > 0:
        table = sorted(table, key=lambda row: row['c1']['l'], reverse=True)

    return {
        "columns": columns,
        "data": table,
    }

def sort_analysis_clusterings(clusterings):
    return sorted(clusterings.items(), key=lambda (k, v): v.global_sort_key)

def _get_unit_and_plt_type(is_spatial=False):
    """ Get the plot type and unit to use"""
    if is_spatial:
        unit = "Spots"
        plt_type = "scatter"
    else:
        unit = "Cells"
        plt_type = "scattergl"
    return unit, plt_type

def analysis_by_clustering(sample_data, spatial=False):
    """
    Get the tSNE (colored by clustering) and diffexp table for each clustering
    """
    if sample_data is None:
        return None

    analysis = sample_data.get_analysis(SingleGenomeAnalysis)
    if analysis is None:
        return None

    tsne_coordinates = analysis.get_tsne().transformed_tsne_matrix
    unit, plt_type = _get_unit_and_plt_type(spatial)


    title = "t-SNE Projection of {} Colored by Clustering".format(unit)
    #pylint: disable=too-many-function-args
    tsne_data = SharedCoordinatePlotCollection(list(tsne_coordinates[:, 0]),
                                                     list(tsne_coordinates[:, 1]),
                                                     pltly.PLOT_CONFIG,
                                                     _make_tsne_umi_layout(title),
                                                     plt_type,
                                                     {"opacity": 0.9, "size": 4})

    diff_exp_tables = Clusterings()
    # order clustering by order: graph, kmeans=2, =3 etc
    clusterings = sort_analysis_clusterings(analysis.clusterings)
    for (clustering_key, clustering) in clusterings:
        # Add t-sne plot for clustering
        tsne_data.add_clustering(clustering_key, clustering)
        # Add differential expression table for clustering
        diffexp = analysis.differential_expression[clustering_key]
        diffexp_tbl = diffexp_table(diffexp, clustering, analysis)
        diftblclust = ClusteringData(clustering_key, clustering, diffexp_tbl)
        diff_exp_tables.add_clustering(diftblclust)
    if not spatial:
        clust_select = ClusteringSelector(TSNE_CLUSTERING_PLOT_HELP, DIFFEXP_TABLE_HELP)
    else:
        clust_select = ClusteringSelector(SPATIAL_TSNE_CLUSTERING_PLOT_HELP,
                                          SPATIAL_DIFFEXP_TABLE_HELP)
    clust_select.right_plots = tsne_data
    clust_select.tables = diff_exp_tables
    return clust_select


def get_umi_color(reads_per_bc):
    vmin, vmax = np.percentile(reads_per_bc, ws_gex_constants.TSNE_TOTALCOUNTS_PRCT_CLIP)
    color = [min(vmax, max(vmin, int(v))) for v in reads_per_bc]
    return color


def umi_on_tsne_plot(sample_data, spatial=False):
    """ UMI count on tSNE plot """
    if sample_data is None:
        return None

    analysis = sample_data.get_analysis(SingleGenomeAnalysis)
    if analysis is None:
        return None

    tsne_coordinates = analysis.get_tsne().transformed_tsne_matrix
    reads_per_bc = analysis.matrix.get_counts_per_bc()
    color = get_umi_color(reads_per_bc)
    unit, plt_type = _get_unit_and_plt_type(spatial)
    title = "t-SNE Projection of {} Colored by UMI Counts".format(unit)
    data = [{
        "name": unit,
        "x": list(tsne_coordinates[:, 0]),
        "y": list(tsne_coordinates[:, 1]),
        "type": plt_type,
        "mode": "markers",
        "marker": {
            "opacity": 0.9,
            "size": 4,
            "color": color,
            "colorscale": "Jet",
            "colorbar": {"title": "UMI counts"}
        },
        "text": ['UMI counts: {:,d}'.format(reads) for reads in reads_per_bc],
    }]

    # Note: the help text has been included in tsne_cluster plot
    layout = TSNE_LAYOUT_CONFIG.copy()
    layout['title'] = title
    umi_tsne_plot = {
            "config": pltly.PLOT_CONFIG,
            "layout": layout,
            "data": data,
    }
    return umi_tsne_plot

def seq_saturation_plot(sample_data, sample_properties, spatial=False):
    chart = {
        'config': pltly.PLOT_CONFIG,
        'layout': {
            'showlegend': False,
            'hovermode': 'closest',
            'xaxis': {
                'title': 'Mean Reads per {}'.format("Cell" if spatial is False else "Spot"),
                "fixedrange": False,
            },
            'yaxis': {
                'title': 'Sequencing Saturation',
                'range': [0, 1],
                "fixedrange": False,
            },
            'shapes': [
                {
                    'type': 'line',
                    'x0': 0,
                    'y0': 0.9,
                    'x1': 0,
                    'y1': 0.9,
                    'line': {
                        'color': 'rgb(128, 128, 128)',
                        'width': 4,
                        'dash': 'dot',
                    },
                },
            ],
        },
        'data': [],  # data entries are built in the function
    }

    kwargs = {
        'metric_suffix': 'subsampled_duplication_frac',
        'show_multi_genome_only': True,
    }

    plot = cr_webshim.plot_subsampled_scatterplot_metric(chart, sample_properties, sample_data, **kwargs)

    if plot:
        return {"seq_saturation_plot" : {
                'plot' : plot,
                'help' : SEQ_SATURATION_PLOT_HELP if spatial is False else SPATIAL_SEQ_SATURATION_PLOT_HELP
                }
            }
    else:
        return None

def median_gene_plot(sample_data, sample_properties, species_list, spatial=False):
    unit = "Cell" if spatial is False else "Spot"
    chart = {
        'config': pltly.PLOT_CONFIG,
        'layout': {
            'showlegend': True,
            'hovermode': 'closest',
            'xaxis': {'title': 'Mean Reads per {}'.format(unit),
                      "fixedrange": False,},
            'yaxis': {'title': 'Median Genes per {}'.format(unit),
                      "fixedrange": False,},
        },
        'data': [],  # data entries are built in the function
    }

    kwargs = {
        'metric_suffix': 'subsampled_filtered_bcs_median_unique_genes_detected',
        'references': species_list,
    }

    plot = cr_webshim.plot_subsampled_scatterplot_metric(chart, sample_properties, sample_data, **kwargs)

    if plot:
        return {"median_gene_plot" : {
                'plot' : plot,
                'help' : MEDIAN_GENE_PLOT_HELP if spatial is False else SPATIAL_MEDIAN_GENE_PLOT_HELP
                }
            }
    else:
        return None
