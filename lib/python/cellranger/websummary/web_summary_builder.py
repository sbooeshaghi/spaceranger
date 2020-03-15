#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
# pylint: disable=too-few-public-methods

from cellranger.websummary.react_components import (WebSummaryData,
                                                    write_html_file)
from cellranger.websummary.analysis_tab import (umi_on_tsne_plot, analysis_by_clustering,
                                                median_gene_plot, seq_saturation_plot,
                                                barnyard_table)
from cellranger.websummary.metrics import MetricAnnotations
from cellranger.webshim.common import load_sample_data

from cellranger.websummary.summary_tab import (FEATURE_BARCODINGS, pipeline_info_table,
                                               CELL_CALLING_METRIC_KEYS, ANTIBODY_CELL_CALLING_METRIC_KEYS,
                                               CELL_CALLING_ALARM_KEYS, ANTIBODY_CELL_CALLING_ALARM_KEYS,
                                               sequencing_table, feature_barcode_sequencing_table, add_data,
                                               hero_metrics,
                                               feature_barcode_application_table, mapping_table,
                                               cell_calling_table, batch_correction_table,
                                               aggregation_table, feature_barcode_aggregation_table)

from cellranger.websummary.sample_properties import CountSampleProperties, SampleDataPaths


CELLRANGER_COMMAND_NAME = "Cell Ranger"

def _build_summary_tab_common(metadata, sample_data, species_list, sample_properties, pipeline, ws_data):
    """ Code to create summary tab shared by spatial and single cell, updates ws_data """
    alarm_list = ws_data.alarms
    summary_tab = ws_data.summary_tab
    add_data(summary_tab, alarm_list, hero_metrics(metadata, sample_data, species_list))
    add_data(summary_tab, alarm_list, pipeline_info_table(sample_data, sample_properties, pipeline))
    add_data(summary_tab, alarm_list, sequencing_table(metadata, sample_data, species_list))
    add_data(summary_tab, alarm_list, mapping_table(metadata, sample_data, species_list))
    return

def _build_analysis_tab_common(sample_data, species_list, sample_properties, ws_data,
                               is_spatial=False):
    """ Code to create analysis tab shared by spatial and single cell """
    alarm_list = ws_data.alarms
    analysis_tab = ws_data.analysis_tab
    add_data(analysis_tab, alarm_list,
             seq_saturation_plot(sample_data, sample_properties, spatial=is_spatial))
    add_data(analysis_tab, alarm_list,
             median_gene_plot(sample_data, sample_properties, species_list, spatial=is_spatial))
    ws_data.clustering_selector = analysis_by_clustering(sample_data, spatial=is_spatial)
    return

def _build_diagnostic_values(sample_data, ws_data):
    """ Some metrics are useful for debugging but not yet ready to be displayed to the customer,
    we hide these in the web summary json inside a "diagnostics" field if we encounter them. """
    # assert isinstance(sample_data, SampleData)
    assert isinstance(ws_data, WebSummaryData)
    diagnostic_metrics = ["tso_frac", "sample_index_bases_with_q30_frac"]
    mets = {}
    for met in diagnostic_metrics:
        if met in sample_data.summary:
            mets[met] = sample_data.summary[met]
    ws_data.diagnostics = mets

def build_web_summary_data_common(sample_properties, sample_data, pipeline, metadata, command):
    """ Produce common data shared by both spatial and cell ranger, VDJ is currently independent """
    is_spatial = command != CELLRANGER_COMMAND_NAME
    wsd = WebSummaryData(sample_properties, command, pipeline)
    species_list = sample_properties.genomes
    _build_summary_tab_common(metadata, sample_data,
                              species_list, sample_properties,
                              pipeline, wsd)
    _build_analysis_tab_common(sample_data,
                               species_list, sample_properties,
                               wsd, is_spatial=is_spatial)
    _build_diagnostic_values(sample_data, wsd)
    return wsd


def build_web_summary_html_sc(filename, sample_properties, sample_data_paths, pipeline):
    assert isinstance(sample_properties, CountSampleProperties)
    assert isinstance(sample_data_paths, SampleDataPaths)
    sample_data = load_sample_data(sample_properties, sample_data_paths)
    metadata = MetricAnnotations()
    species_list = sample_properties.genomes
    web_sum_data = build_web_summary_data_common(sample_properties,
                                                 sample_data, pipeline, metadata,
                                                 CELLRANGER_COMMAND_NAME)
    # Single cell specific stuff
    add_data(web_sum_data.summary_tab, web_sum_data.alarms,
             aggregation_table(metadata, sample_data, sample_properties))
    add_data(web_sum_data.summary_tab, web_sum_data.alarms,
             batch_correction_table(metadata, sample_data, species_list))


    # TODO: Detecting antibody-only/GEX-Less was done by two different conditions in two different
    # stages before.  In `count` the check was if 'cells' not in web_sum_data.summary_tab
    # and in the check was `aggr` was if len(species_list) == 0:
    # currently trying the condition below which should count for both.
    if len(sample_properties.genomes) ==0: # This is the case when there was no GEX input library
        add_data(web_sum_data.summary_tab, web_sum_data.alarms,
                  cell_calling_table(metadata, sample_data, sample_properties, species_list,
                         ANTIBODY_CELL_CALLING_METRIC_KEYS, ANTIBODY_CELL_CALLING_ALARM_KEYS))
    else:
        add_data(web_sum_data.summary_tab, web_sum_data.alarms,
                 cell_calling_table(metadata, sample_data,
                                    sample_properties, species_list, CELL_CALLING_METRIC_KEYS,
                                    CELL_CALLING_ALARM_KEYS))

    for fb in FEATURE_BARCODINGS:
        # Not all three of theses will be present

        # feature barcoding sequencing info (Count)
        add_data(web_sum_data.summary_tab, web_sum_data.alarms,
                 feature_barcode_sequencing_table(metadata, sample_data, species_list, fb))
        # feature barcoding application metric (Count)
        add_data(web_sum_data.summary_tab, web_sum_data.alarms,
                 feature_barcode_application_table(metadata, sample_data, species_list, fb))
        # aggregation metrics (AGGR)
        add_data(web_sum_data.summary_tab, web_sum_data.alarms,
                 feature_barcode_aggregation_table(metadata, sample_data, sample_properties,
                                                   fb))

    # t-SNE plot appears as a constant left plot in the Cell Ranger
    # clustering selector
    if web_sum_data.clustering_selector:
        web_sum_data.clustering_selector.left_plots = umi_on_tsne_plot(sample_data)
    # Barnyard
    add_data(web_sum_data.analysis_tab, web_sum_data.alarms,
             barnyard_table(metadata, sample_data, sample_properties, species_list))
    write_html_file(filename, web_sum_data)
    return
