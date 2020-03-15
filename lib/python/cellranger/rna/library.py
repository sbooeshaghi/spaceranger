#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
import json
import re

GENE_EXPRESSION_LIBRARY_TYPE = 'Gene Expression'
CRISPR_LIBRARY_TYPE = 'CRISPR Guide Capture'
ANTIBODY_LIBRARY_TYPE = 'Antibody Capture'
FEATURETEST_LIBRARY_TYPE = 'FEATURETEST'
CUSTOM_LIBRARY_TYPE = 'Custom'
VDJ_LIBRARY_TYPE = 'VDJ'
DEFAULT_LIBRARY_TYPE = GENE_EXPRESSION_LIBRARY_TYPE

MULTI_REFS_PREFIX = 'multi'

CUSTOM_METRIC_PREFIX = CUSTOM_LIBRARY_TYPE
DISPLAY_PREFIX_CUSTOM = 'Custom'
CRISPR_METRIC_PREFIX = 'CRISPR'
DISPLAY_PREFIX_CRISPR = 'CRISPR:'
ANTIBODY_METRIC_PREFIX = 'ANTIBODY'
DISPLAY_PREFIX_ANTIBODY = 'Antibody:'
FEATURETEST_METRIC_PREFIX = 'FEATURETEST'

LIBRARY_TYPE = 'library_type'

RECOGNIZED_FEATURE_TYPES = [GENE_EXPRESSION_LIBRARY_TYPE,
                            CRISPR_LIBRARY_TYPE,
                            ANTIBODY_LIBRARY_TYPE,
                            FEATURETEST_LIBRARY_TYPE,
                           ]
FEATURE_LIBRARY_TYPES = [CRISPR_LIBRARY_TYPE,
                         ANTIBODY_LIBRARY_TYPE,
                         FEATURETEST_LIBRARY_TYPE]

metric_prefix_map = {
    GENE_EXPRESSION_LIBRARY_TYPE:'',
    CRISPR_LIBRARY_TYPE: CRISPR_METRIC_PREFIX,
    ANTIBODY_LIBRARY_TYPE: ANTIBODY_METRIC_PREFIX,
    CUSTOM_LIBRARY_TYPE: CUSTOM_METRIC_PREFIX,
    FEATURETEST_LIBRARY_TYPE: FEATURETEST_METRIC_PREFIX
}

report_prefix_map = {
    GENE_EXPRESSION_LIBRARY_TYPE: '',
    CRISPR_LIBRARY_TYPE: DISPLAY_PREFIX_CRISPR,
    ANTIBODY_LIBRARY_TYPE: DISPLAY_PREFIX_ANTIBODY,
    CUSTOM_LIBRARY_TYPE: DISPLAY_PREFIX_CUSTOM
}

# 'target_set_name' should be a key in library_info
TARGET_SET_KEY = 'target_set_name'
DEFAULT_TARGET_SET = ''

def _get_prefix(lib_type, sep, prefix_map):
    if lib_type == GENE_EXPRESSION_LIBRARY_TYPE:
        return ''
    else:
        value = prefix_map.get(lib_type, lib_type)
        return value + sep

def get_library_type_metric_prefix(lib_type):
    """Get the metric prefix for a given library type.
    """
    return _get_prefix(lib_type, '_', metric_prefix_map)


def get_library_type_report_prefix(lib_type):
    """Gets the prefix to be used in displayed reports """
    return _get_prefix(lib_type, ' ', report_prefix_map)

def add_multi_prefix(metric):
    """Appends the prefix for cumulative metrics onto a metric name """
    return "{}_{}".format(MULTI_REFS_PREFIX, metric)

def add_species_prefix(species, metric):
    """Append the species/genome name to the front of a metric name """
    return '{}_{}'.format(species, metric)

def get_bam_library_info(bam):
    """Get the library info from a BAM's comment lines.
    Args:
      bam (pysam.AlignmentFile): BAM file
    Returns:
      list of dicts
    """
    comments = bam.header['CO']
    libraries = []
    for comment in comments:
        m = re.match(r'^library_info:(.+)$', comment)
        if m:
            libraries.append(json.loads(m.group(1)))
    return libraries

def has_genomes(library_type):
    """Do genomes make sense for a library type"""
    return library_type == GENE_EXPRESSION_LIBRARY_TYPE

def sorted_library_types(library_info):
    """Sorted list of unique library types in library_info"""
    return sorted(set(lib[LIBRARY_TYPE] for lib in library_info))

def has_target_set(library):
    return TARGET_SET_KEY in library and library[TARGET_SET_KEY] != DEFAULT_TARGET_SET
