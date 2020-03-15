#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
import itertools
import numpy as np
import tables
import cellranger.analysis.pca as cr_pca
import cellranger.analysis.lm_umap as cr_umap
import cellranger.analysis.io as analysis_io
import cellranger.h5_constants as h5_constants
import cellranger.analysis.constants as analysis_constants
import cellranger.matrix as cr_matrix
import cellranger.io as cr_io
import cellranger.rna.library as rna_library

__MRO__ = """
stage RUN_UMAP(
    in  h5    matrix_h5,
    in  h5    pca_h5,
    in  bool  skip,
    in  int   random_seed,
    in  int   n_neighbors,
    in  int   input_pcs,
    in  int   max_dims,
    in  float min_dist,
    in  string metric,
    in  bool  is_antibody_only,
    out h5    umap_h5,
    out path  umap_csv,
    src py    "stages/analyzer/run_umap",
) split using (
    in  int   umap_dims,
    in  string feature_type,
)
"""

def split(args):
    if args.skip:
        return {'chunks': []}

    feature_ref = cr_matrix.CountMatrix.load_feature_ref_from_h5_file(args.matrix_h5)
    feature_types = sorted(list(set(fd.feature_type for fd in feature_ref.feature_defs)))

    chunks = []
    matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.matrix_h5)
    mem_gb = max(matrix_mem_gb, h5_constants.MIN_MEM_GB)

    min_umap_dims = analysis_constants.TSNE_N_COMPONENTS
    max_umap_dims = args.max_dims if args.max_dims is not None else min_umap_dims
    for umap_dims, feature_type in itertools.product(range(min_umap_dims, max_umap_dims+1),
                                                     feature_types):
        chunks.append({
            'umap_dims': umap_dims,
            'feature_type': feature_type,
            '__mem_gb': mem_gb,
            '__threads': 4,
        })
    return {'chunks': chunks, 'join': {'__mem_gb' : 1}}

def lower_no_space(s):
    return s.replace(" ", "_").lower()

def get_umap_name(feature_type, n_components):
    return '%s_%d-d' % (lower_no_space(feature_type), n_components)

def get_umap_key(feature_type, n_components):
    if feature_type == rna_library.DEFAULT_LIBRARY_TYPE:
        # Preserve backward HDF5 compatibility with pre-3.0 HDF5 files
        #   where the CSV directory was named "2_components" and the HDF5 dataset was named "_2"
        return str(n_components)
    else:
        return '%s_%d' % (lower_no_space(feature_type), n_components)

def main(args, outs):
    if args.skip:
        return

    umap_dims = args.umap_dims
    if args.is_antibody_only:
        matrix = cr_matrix.CountMatrix.load_h5_file(args.matrix_h5)
        matrix = matrix.select_features_by_type(rna_library.ANTIBODY_LIBRARY_TYPE)
        matrix.m.data = np.log2(1 + matrix.m.data)
        umap_input = matrix.m.transpose().todense()

        name = get_umap_name(rna_library.ANTIBODY_LIBRARY_TYPE, args.umap_dims)
        key = str(args.umap_dims)

        umap = cr_umap.run_umap(umap_input, name=name, key=key, input_pcs=args.input_pcs, n_neighbors=args.n_neighbors,
            min_dist=args.min_dist, metric=args.metric, umap_dims=umap_dims, random_state=args.random_seed)

        filters = tables.Filters(complevel = h5_constants.H5_COMPRESSION_LEVEL)
        with tables.open_file(outs.umap_h5, 'w', filters = filters) as f:
            cr_umap.save_umap_h5(umap, f)

        cr_umap.save_umap_csv(umap, matrix.bcs, outs.umap_csv)
        return

    if args.feature_type == rna_library.GENE_EXPRESSION_LIBRARY_TYPE:
        # Use PCA for gene expression
        pca = cr_pca.load_pca_from_h5(args.pca_h5)
        umap_input = pca.transformed_pca_matrix
    else:
        # Use feature space for other feature types
        # Assumes other feature types are much lower dimension than gene expression
        matrix = cr_matrix.CountMatrix.load_h5_file(args.matrix_h5)
        matrix = matrix.select_features_by_type(args.feature_type)
        matrix.m.data = np.log2(1 + matrix.m.data)
        umap_input = matrix.m.transpose().todense()

    name = get_umap_name(args.feature_type, args.umap_dims)
    key = get_umap_key(args.feature_type, args.umap_dims)

    umap = cr_umap.run_umap(umap_input, name=name, key=key, input_pcs=args.input_pcs, n_neighbors=args.n_neighbors,
            min_dist=args.min_dist, metric=args.metric, umap_dims=umap_dims, random_state=args.random_seed)

    filters = tables.Filters(complevel=h5_constants.H5_COMPRESSION_LEVEL)
    with tables.open_file(outs.umap_h5, 'w', filters=filters) as f:
        cr_umap.save_umap_h5(umap, f)

    matrix_bcs = cr_matrix.CountMatrix.load_bcs_from_h5_file(args.matrix_h5)
    cr_umap.save_umap_csv(umap, matrix_bcs, outs.umap_csv)

def join(args, outs, chunk_defs, chunk_outs):
    if args.skip:
        return

    chunk_h5s = [chunk_out.umap_h5 for chunk_out in chunk_outs]
    chunk_csv_dirs = [chunk_out.umap_csv for chunk_out in chunk_outs]
    analysis_io.combine_h5_files(chunk_h5s, outs.umap_h5, [analysis_constants.ANALYSIS_H5_UMAP_GROUP])
    for csv_dir in chunk_csv_dirs:
        cr_io.copytree(csv_dir, outs.umap_csv, allow_existing=True)
