#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#
""" Constants for feature-barcoding analysis """
from cellranger.rna.library import CRISPR_METRIC_PREFIX, ANTIBODY_METRIC_PREFIX


PREFIX_FROM_FEATURE_TYPE = {
    CRISPR_METRIC_PREFIX: CRISPR_METRIC_PREFIX,
    ANTIBODY_METRIC_PREFIX: ANTIBODY_METRIC_PREFIX,
}
