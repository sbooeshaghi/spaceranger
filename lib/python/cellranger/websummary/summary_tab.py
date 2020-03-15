#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

# Here, you could find everything shown in the summary tab

import cellranger.webshim.common as cr_webshim
import cellranger.constants as cr_constants
import cellranger.rna.library as rna_library
import cellranger.reference as cr_reference
import cellranger.webshim.constants.shared as shared_constants
import cellranger.vdj.constants as vdj_constants
import cellranger.websummary.plotly_tools as pltly
import cellranger.websummary.sample_properties as wsp

ALARMS = shared_constants.ALARMS

# Feature barcoding internel name <-> display name
FEATURE_BARCODINGS = [rna_library.CRISPR_METRIC_PREFIX, rna_library.ANTIBODY_METRIC_PREFIX, rna_library.CUSTOM_METRIC_PREFIX]

FB_DISPLAY_NAME = {
    rna_library.CRISPR_METRIC_PREFIX: "CRISPR",
    rna_library.ANTIBODY_METRIC_PREFIX: "Antibody",
    rna_library.CUSTOM_METRIC_PREFIX: "Custom Feature"
}

RANK_PLOT_HELP = [["Barcode Rank Plot",
                   ["The plot shows the count of filtered UMIs mapped to each barcode.  As barcodes are not determined to be cell-associated strictly based on their UMI count, but instead are determined by their expression profiles, some regions of the graph contain both cell-associated and background-associated barcodes.  The color of the graph in these regions is based on the local density of barcodes that are cell-associated."]]]


# metric and alarm keys for cell calling (GEX and feature barcoding)
CELL_CALLING_METRIC_KEYS = [
    'filtered_bcs',
    'filtered_bcs_conf_mapped_barcoded_reads_cum_frac',
    'multi_transcriptome_total_raw_reads_per_filtered_bc',
    'filtered_bcs_median_unique_genes_detected',
    'filtered_bcs_total_unique_genes_detected',
    'filtered_bcs_median_counts',
]

ANTIBODY_CELL_CALLING_METRIC_KEYS = [
    'ANTIBODY_filtered_bcs_transcriptome_union',
    'ANTIBODY_multi_transcriptome_total_raw_reads_per_filtered_bc',
]

CELL_CALLING_ALARM_KEYS = [
    'filtered_bcs_conf_mapped_barcoded_reads_cum_frac'
]

# should become ANTIBODY_multi_filtered_bcs_conf_mapped_barcoded_reads_cum_frac when we figure out the appropriate alarm thresholds
ANTIBODY_CELL_CALLING_ALARM_KEYS = None


# metric keys for sequencing (GEX and feature barcoding)
SEQUENCING_METRIC_KEYS = [
    'total_reads',
    'good_bc_frac',
    'good_umi_frac',
    'multi_cdna_pcr_dupe_reads_frac',
    'bc_bases_with_q30_frac',
    'read_bases_with_q30_frac',
    'read2_bases_with_q30_frac',
    'sample_index_bases_with_q30_frac',
    'umi_bases_with_q30_frac',
]

SEQUENCING_ALARM_KEYS = [
    'good_bc_frac',
    'good_umi_frac',
    'bc_bases_with_q30_frac',
    'read_bases_with_q30_frac',
    'sample_index_bases_with_q30_frac',
    'umi_bases_with_q30_frac',
]

AGGREGATION_METRIC_KEYS = [
    "frac_reads_kept",
    "pre_normalization_raw_reads_per_filtered_bc",
    "pre_normalization_cmb_reads_per_filtered_bc",
]

# metric keys for feature barcoding application
FB_APP_METRIC_KEYS = {
    'CRISPR': [
        'CRISPR_feature_bc_extracted_frac',
        'CRISPR_multi_transcriptome_conf_mapped_reads_frac',
        'CRISPR_frac_feature_reads_usable',
        'CRISPR_feature_reads_usable_per_cell',
        'CRISPR_unrecognized_feature_bc_frac',
        'CRISPR_feature_reads_in_cells',
        'CRISPR_frac_cells_with_protospacer',
        'CRISPR_frac_cells_with_multiple_protospacer',
        'CRISPR_multi_filtered_bcs_median_counts'
    ],
    'ANTIBODY': [
        'ANTIBODY_multi_transcriptome_conf_mapped_reads_frac',
        'ANTIBODY_frac_feature_reads_usable',
        'ANTIBODY_feature_reads_usable_per_cell',
        'ANTIBODY_reads_lost_to_highly_corrected_GEMs',
        'ANTIBODY_unrecognized_feature_bc_frac',
        'ANTIBODY_feature_reads_in_cells',
        'ANTIBODY_multi_filtered_bcs_median_counts',
    ],
    'Custom': [
        'Custom_multi_transcriptome_conf_mapped_reads_frac',
        'Custom_frac_feature_reads_usable',
        'Custom_feature_reads_usable_per_cell',
        'Custom_unrecognized_feature_bc_frac',
        'Custom_feature_reads_in_cells',
        'Custom_multi_filtered_bcs_median_counts'
    ]
}


def add_data(websummary_json, alarm_list, input_data):
    """Adds data to global dictionary"""
    if input_data is None:
        return

    if ALARMS in input_data:
        alarm_list.extend(input_data[ALARMS])
        del input_data[ALARMS]
    websummary_json.update(input_data)
    return


def hero_metrics(metadata, sample_data, species_list):
    if sample_data is None or sample_data.summary is None:
        return None

    # For FB-only web summaries, only use cell counts and total antibody reads
    if 'filtered_bcs_transcriptome_union' not in sample_data.summary:
        data = {}
        for key in ANTIBODY_CELL_CALLING_METRIC_KEYS:
            metrics = metadata.gen_metric_list(
                sample_data.summary, [key], species_list)
            for metric in metrics:
                # remove the prefix, so the key name matches with the input json expects
                new_key = metric.key.replace('ANTIBODY_', '')
                data[new_key] = metric.gen_metric_dict()
        return data

    data = {}
    for key in ['filtered_bcs_transcriptome_union', 'multi_transcriptome_total_raw_reads_per_filtered_bc']:
        metrics = metadata.gen_metric_list(
            sample_data.summary, [key], species_list)
        for metric in metrics:
            data[metric.key] = metric.gen_metric_dict()

    is_barnyard = len(species_list) > 1
    if not is_barnyard:
        median_unique_genes = 'filtered_bcs_median_unique_genes_detected'
        metrics = metadata.gen_metric_list(
            sample_data.summary, [median_unique_genes], species_list)
        for metric in metrics:
            data[median_unique_genes] = metric.gen_metric_dict()

    alarm_keys = ['filtered_bcs_transcriptome_union']
    alarms = metadata.gen_metric_list(
        sample_data.summary, alarm_keys, species_list)
    new_alarms = [metric.alarm_dict for metric in alarms if metric.alarm_dict]
    if new_alarms:
        data[ALARMS] = new_alarms

    return data


def pipeline_info_table(sample_data, sample_properties, pipeline):
    """Generates a table of general pipeline information.
    """
    assert isinstance(sample_properties, wsp.SampleProperties)
    if sample_data is None or sample_data.summary is None:
        return None

    rows = [
        ['Sample ID', sample_properties.sample_id],
        ['Sample Description', sample_properties.sample_desc],
        ['Chemistry', sample_data.summary.get('chemistry_description')],
    ]

    if sample_data.summary.get('spatial_slide_info', None) is not None:
        rows.append(['Slide Serial Number', sample_data.summary['spatial_slide_info']])

    if isinstance(sample_properties, wsp.ExtendedCountSampleProperties):
        if sample_properties.reference_path:
            rows.append(['Reference Path', sample_properties.reference_path])
        # This was meant to be enabled in 3.1 but due to a bug was not included./
        #if sample_properties.barcode_whitelist:
        #    rows.append(
        #        ['Barcode Whitelist', sample_properties.barcode_whitelist])

    # Find references in the summary
    if isinstance(sample_properties, wsp.AggrCountSampleProperties) and \
            not sample_properties.genomes:
        rows.append([cr_constants.REFERENCE_TYPE, 'Not applicable for aggr with feature barcoding-only samples'])
    elif pipeline in [shared_constants.PIPELINE_AGGR, shared_constants.PIPELINE_REANALYZE]:
        genomes = sample_properties.genomes
        if genomes is not None:
            rows.append([cr_constants.REFERENCE_TYPE,
                         cr_reference.get_ref_name_from_genomes(genomes)])
    else:
        reference_metric_prefixes = [cr_constants.REFERENCE_METRIC_PREFIX,
                                     vdj_constants.REFERENCE_METRIC_PREFIX]
        # Find all references in the summary
        for prefix in reference_metric_prefixes:
            ref_type_key = '%s%s' % (prefix, cr_constants.REFERENCE_TYPE_KEY)
            if ref_type_key in sample_data.summary:
                ref_type = sample_data.summary[ref_type_key]

                ref_version_key = '%s%s' % (
                    prefix, cr_constants.REFERENCE_VERSION_KEY)
                if ref_version_key in sample_data.summary:
                    ref_version = '-%s' % sample_data.summary.get(
                        ref_version_key)
                else:
                    ref_version = ''

                ref_name_key = '%s%s' % (
                    prefix, cr_constants.REFERENCE_GENOMES_KEY)
                if ref_name_key in sample_data.summary:
                    ref_name = sample_data.summary.get(ref_name_key)
                    if isinstance(ref_name, list):
                        ref_name = cr_reference.get_ref_name_from_genomes(
                            ref_name)

                    rows.append([ref_type, '%s%s' % (ref_name, ref_version)])

    # add pipeline version
    rows.append(['Pipeline Version', sample_properties.version])

    pipeline_info = {
        "header": ["Sample"],
        "rows": rows,
    }
    return {"pipeline_info_table": pipeline_info}


def create_table_with_alarms(table_key, title, metric_keys, alarm_keys,
                             metadata, sample_data, species_list):
    """
    Sequencing info for GEX
    """
    if sample_data is None or sample_data.summary is None:
        return None

    data_dict = {}

    metrics = metadata.gen_metric_list(
        sample_data.summary, metric_keys, species_list)
    if metrics:
        data_dict["help"] = {'title': title,
                             'data': metadata.gen_metric_helptext(metric_keys)}
        data_dict["table"] = {
            "rows": [[metric.name, metric.value_string] for metric in metrics]}

    if not data_dict:
        return None

    result = {table_key: data_dict}

    # Alerts.
    if alarm_keys:
        alarms = metadata.gen_metric_list(
            sample_data.summary, alarm_keys, species_list)
        # If a metric is from a barnyard and the cumulative version of the metric should be tested
        # we do not check the non-cumulative metrics
        alarms = [x for x in alarms if not (x.is_barnyard and
                                            x.parent_metric_info.include_cumulative and
                                            not x.is_cumulative)]
        new_alarms = [
            metric.alarm_dict for metric in alarms if metric.alarm_dict]

        if new_alarms:
            result[ALARMS] = new_alarms

    return result


def sequencing_table(metadata, sample_data, species_list):
    """ Sequencing info for GEX """
    return create_table_with_alarms("sequencing", "Sequencing",
                                    SEQUENCING_METRIC_KEYS, SEQUENCING_ALARM_KEYS,
                                    metadata, sample_data, species_list)


def feature_barcode_sequencing_table(metadata, sample_data, species_list, feature_barcode):
    metric_keys = ['{}_{}'.format(feature_barcode, i)
                   for i in SEQUENCING_METRIC_KEYS]
    alarm_keys = ['{}_{}'.format(feature_barcode, i)
                  for i in SEQUENCING_ALARM_KEYS]

    return create_table_with_alarms("{}_sequencing".format(feature_barcode.upper()),
                                    "{} Sequencing".format(
                                        FB_DISPLAY_NAME[feature_barcode]),
                                    metric_keys, alarm_keys,
                                    metadata, sample_data, species_list)


def feature_barcode_application_table(metadata, sample_data, species_list, feature_barcode):
    """ Feature barcoding application metric """
    return create_table_with_alarms("{}_application".format(feature_barcode.upper()),
                                    "{} Application".format(
                                        FB_DISPLAY_NAME[feature_barcode]),
                                    FB_APP_METRIC_KEYS[feature_barcode], None,
                                    metadata, sample_data, species_list)


def mapping_table(metadata, sample_data, species_list):
    """ Mapping info table """
    metric_keys = [
        'genome_mapped_reads_frac',
        'genome_conf_mapped_reads_frac',
        'intergenic_conf_mapped_reads_frac',
        'intronic_conf_mapped_reads_frac',
        'exonic_conf_mapped_reads_frac',
        'transcriptome_conf_mapped_reads_frac',
        'on_target_conf_mapped_reads_frac',
        'antisense_reads_frac',
    ]

    alarm_keys = [
        'transcriptome_conf_mapped_reads_frac',
        'antisense_reads_frac',
    ]

    return create_table_with_alarms("mapping", "Mapping",
                                    metric_keys, alarm_keys,
                                    metadata, sample_data, species_list)

def spot_calling_table(metadata, sample_data, species_list, zoom_images, metric_keys, alarm_keys):
    """cell calling data (table and plot). """
    #TODO: Barnyard not currently in spatial
    is_barnyard = len(species_list) > 1
    if is_barnyard:
        metric_keys.insert(0, 'filtered_bcs_transcriptome_union')

    data_dict = create_table_with_alarms("cells", "Spots",
                                         metric_keys, alarm_keys,
                                         metadata, sample_data, species_list)

    if data_dict is None:
        return None

    # add image
    data_dict["cells"]['zoom_images'] = zoom_images
    return data_dict

def cell_calling_table(metadata, sample_data, sample_properties, species_list, metric_keys, alarm_keys):
    """cell calling data (table and plot). """

    is_barnyard = len(species_list) > 1
    if is_barnyard:
        metric_keys.insert(0, 'filtered_bcs_transcriptome_union')

    table_dict = create_table_with_alarms("cells", "Cells",
                                         metric_keys, alarm_keys,
                                         metadata, sample_data, species_list)

    if table_dict is None:
        return None

    # the data we are interested in is in data_dict["cells"]
    data_dict = table_dict["cells"]
    to_return = {}
    # Be sure to bubble up alarms
    if ALARMS in table_dict:
        to_return[ALARMS] = table_dict[ALARMS]


    chart = {
        'config': pltly.PLOT_CONFIG,
        'layout': {
            'title': 'Barcode Rank Plot',
            'xaxis': {
                'title': 'Barcodes',
                'type': 'log',
                "showline": True,
                "zeroline": False,
                "fixedrange": False,
            },
            'yaxis': {
                'title': 'UMI counts',
                'type': 'log',
                "showline": True,
                "zeroline": False,
                "fixedrange": False,
            },
            'hovermode': 'closest',
        },
        'data': [],
    }

    knee_plot = cr_webshim.plot_barcode_rank(
        chart, sample_properties, sample_data)
    if knee_plot:
        data_dict['barcode_knee_plot'] = knee_plot
        data_dict['help']['data'] = data_dict['help']['data'] + RANK_PLOT_HELP
    to_return["cells"] = data_dict
    return to_return


def batch_correction_table(metadata, sample_data, species_list):
    metric_keys = [
        "batch_effect_score_before_correction",
        "batch_effect_score_after_correction",
    ]

    return create_table_with_alarms("batch_correction", "Chemistry Batch Correction",
                                    metric_keys, None,
                                    metadata, sample_data, species_list)


def aggregation_table(metadata, sample_data, sample_properties):
    """
    Report normalization metrics in aggr. The trick here is to use the batch prefix as the
    a species/genome prefix, and define these metric as species_specific.

    TODO: the above trick doesn't generate good web summary if it's a barnyard aggr sample.
    """
    if not isinstance(sample_properties, wsp.AggrCountSampleProperties):
        return None
    metric_keys = [
        "pre_normalization_total_reads",
        "post_normalization_total_reads",
        "pre_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc",
        "post_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc",
    ] + AGGREGATION_METRIC_KEYS

    alarm_keys = [
        'lowest_frac_reads_kept',
    ]
    batches = sample_properties.agg_batches
    return create_table_with_alarms("aggregation", "Aggregation",
                                    metric_keys, alarm_keys,
                                    metadata, sample_data, batches)


def feature_barcode_aggregation_table(metadata, sample_data, sample_properties, feature_barcode):
    if not isinstance(sample_properties, wsp.AggrCountSampleProperties):
        return None
    metric_keys = ['{}_{}'.format(feature_barcode, i)
                   for i in AGGREGATION_METRIC_KEYS]
    batches = sample_properties.agg_batches
    return create_table_with_alarms("{}_aggregation".format(feature_barcode),
                                    "{} Aggregation".format(
                                        FB_DISPLAY_NAME[feature_barcode]),
                                    metric_keys, None,
                                    metadata, sample_data, batches)
