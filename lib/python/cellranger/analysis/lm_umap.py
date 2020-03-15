#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.analysis.io as analysis_io
import cellranger.analysis.constants as analysis_constants
import cellranger.io as cr_io

import collections
import numpy as np
import os
import umap
import martian

UMAP = collections.namedtuple('UMAP', ['transformed_umap_matrix',
                                       'name',  # Human readable form
                                       'key',   # Machine queryable form, must be unique
                                       ])


def run_umap(transformed_pca_matrix, name='UMAP', key='UMAP', input_pcs=None, n_neighbors=None,
             min_dist=None, metric=None, umap_dims=None, random_state=None):

    if umap_dims is None:
        umap_dims = analysis_constants.UMAP_N_COMPONENTS

    if n_neighbors is None:
        n_neighbors = analysis_constants.UMAP_DEFAULT_N_NEIGHBORS

    if min_dist is None:
        min_dist = analysis_constants.UMAP_MIN_DIST

    if metric is None:
        metric = analysis_constants.UMAP_DEFAULT_METRIC
    else:
        metric = metric.encode("utf-8")

    if random_state is None:
        random_state = analysis_constants.RANDOM_STATE

    if input_pcs is not None:
        transformed_pca_matrix = transformed_pca_matrix[:, :input_pcs]

    try:
        umap_reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist,
                                 n_components=umap_dims, metric=metric,
                                 random_state=np.random.RandomState(random_state))

        transformed_umap_matrix = umap_reducer.fit_transform(transformed_pca_matrix)
    except ValueError:
        martian.log_info('Failed to run UMAP with default setting on this data. Trying random initialization.')

        umap_reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist,
                                 n_components=umap_dims, metric=metric,
                                 random_state=np.random.RandomState(random_state),
                                 init='random')

        transformed_umap_matrix = umap_reducer.fit_transform(transformed_pca_matrix)
    except:
        martian.log_info('Failed to run UMAP on this data. Filling in 0s as UMAP embeddings.')
        transformed_umap_matrix = np.zeros((transformed_pca_matrix.shape[0], umap_dims), dtype=np.float32)

    return UMAP(transformed_umap_matrix, name=name, key=key)


def save_umap_csv(umap_obj, barcodes, base_dir):
    """Save a UMAP object to CSV"""
    # Preserve backward compatibility with pre-3.0 CSV files
    #   where the CSV directory was named "2_components" and the HDF5 dataset was named "_2"
    key = umap_obj.key + '_components'

    umap_dir = os.path.join(base_dir, key)
    cr_io.makedirs(umap_dir, allow_existing=True)

    matrix_fn = os.path.join(umap_dir, 'projection.csv')
    n_umap_components = umap_obj.transformed_umap_matrix.shape[1]
    matrix_header = ['Barcode'] + ['UMAP-%d' % (i+1) for i in range(n_umap_components)]
    analysis_io.save_matrix_csv(matrix_fn, umap_obj.transformed_umap_matrix, matrix_header, barcodes)


def save_umap_h5(umap_obj, f):
    """Save a UMAP object to HDF5"""
    group = f.create_group(f.root, analysis_constants.ANALYSIS_H5_UMAP_GROUP)
    analysis_io.save_h5(f, group, umap_obj.key, umap_obj)
