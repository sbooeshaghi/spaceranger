#!/usr/bin/env python
#
# Copyright (c) 2018 10x Genomics, Inc. All rights reserved.
#
import cPickle
from collections import defaultdict
from itertools import izip
import sys

import martian
import numpy as np
import tables

import cellranger.utils as cr_util
import cellranger.analysis.pca as cr_pca
import cellranger.h5_constants as h5_constants
import cellranger.analysis.constants as analysis_constants
import cellranger.analysis.batch_correction as cr_batch_correction

__MRO__  = """
stage CORRECT_CHEMISTRY_BATCH(
    in  pickle dimred_matrix,
    in  pickle matrix_barcode_feature_info,
    in  map[]  library_info,
    in  int    cbc_knn,
    in  float  cbc_alpha,
    in  float  cbc_sigma,
    in  bool   cbc_realign_panorama,
    in  bool   skip,
    out float  batch_score_before_correction,
    out float  batch_score_after_correction,
    out h5     aligned_pca_h5,
    out path   aligned_pca_csv,
    src py     "stages/analyzer/correct_chemistry_batch",
) split (
    in  int    batch_id,
    in  map    batch_to_bc_indices,
    in  pickle ordered_dimred_matrix,
    in  pickle idx_to_batch_id,
    in  bool   need_reorder_barcode,
    in  pickle barcode_reorder_index,
    out binary batch_nearest_neighbor,
)
"""

NUM_ENTRIES_PER_MEM_GB = 2000000
JOIN_MEM_GB = 64

def option(arg, default):
    return arg if arg is not None else default


def split(args):
    if args.skip:
        return {'chunks': []}

    gg_id_to_batch_id, batch_id_to_name = {}, {}

    for lib in args.library_info:
        gg_id_to_batch_id[lib['gem_group']] = lib['batch_id']
        batch_id_to_name[lib['batch_id']] = lib['batch_name']

    # load the barcodes
    with open(args.matrix_barcode_feature_info) as fp:
        bc_feature_info = cPickle.load(fp)
        bcs = bc_feature_info.get('barcodes')

    batch_ids = np.array([gg_id_to_batch_id[cr_util.split_barcode_seq(bc)[1]] for bc in bcs])

    with open(args.dimred_matrix) as fp:
        dimred_matrix = cPickle.load(fp)

    # re-order matrix such that barcodes from same batch are grouped together
    new_bc_indices = None
    batch_to_bc_indices = []
    idx_to_batch_id = np.full(dimred_matrix.shape[0], 0, dtype=np.int8)

    base = 0
    for b_id in range(len(batch_id_to_name)):
        batch_bc_indices = np.where(batch_ids == b_id)[0]
        if batch_bc_indices.shape[0] == 0:
            continue

        new_bc_indices = batch_bc_indices if new_bc_indices is None else np.append(new_bc_indices, batch_bc_indices)
        batch_to_bc_indices.append((base, base + batch_bc_indices.shape[0]))
        idx_to_batch_id[base : base + batch_bc_indices.shape[0]] = b_id
        base += len(batch_bc_indices)

    # 1. check if needs re-order; 2. if needs re-order, store the original order
    need_reorder_barcode = (not np.all(np.diff(new_bc_indices)>=0))
    if need_reorder_barcode:
        dimred_matrix = dimred_matrix[new_bc_indices]
        barcode_reorder_index = np.argsort(new_bc_indices)

        barcode_reorder_index_file = martian.make_path('barcode_reorder_index.pickle')
        with open(barcode_reorder_index_file, 'wb') as fp:
            cPickle.dump(barcode_reorder_index, fp, cPickle.HIGHEST_PROTOCOL)

        ordered_dimred_matrix_file = martian.make_path('ordered_dimred_matrix.pickle')
        with open(ordered_dimred_matrix_file, 'wb') as fp:
            cPickle.dump(dimred_matrix, fp, cPickle.HIGHEST_PROTOCOL)
    else:
        barcode_reorder_index_file, ordered_dimred_matrix_file = None, None

    idx_to_batch_id_file = martian.make_path('idx_to_batch_id.pickle')
    with open(idx_to_batch_id_file, 'wb') as fp:
        cPickle.dump(idx_to_batch_id, fp, cPickle.HIGHEST_PROTOCOL)

    nitem = dimred_matrix.shape[0]
    nbatch = len(batch_to_bc_indices)
    cbc_knn = option(args.cbc_knn, analysis_constants.CBC_KNN)
    matrix_mem_gb = sys.getsizeof(dimred_matrix) / 1e9 # float(nitem * ndim) / NUM_ENTRIES_PER_MEM_GB
    # 72 for size of tuple, 32 * 2 for size of 2 np.int64's, and 40% for inefficient dictionaries
    nn_mem_gb = 1.4 * nbatch * nitem * cbc_knn * (72 + 2 * 32) / 1e9
    # presuming all in one batch, dimred_matrix, cur_matrix, ref_matrix
    main_mem_gb = max(int(3.0 * matrix_mem_gb + nn_mem_gb + 1.0), h5_constants.MIN_MEM_GB)

    chunks = []
    for batch_id in xrange(len(batch_to_bc_indices)):
        chunks.append({
            '__mem_gb': main_mem_gb,
            'batch_id': batch_id,
            'batch_to_bc_indices': batch_to_bc_indices,
            'ordered_dimred_matrix': ordered_dimred_matrix_file,
            'idx_to_batch_id': idx_to_batch_id_file,
            'need_reorder_barcode': need_reorder_barcode,
            'barcode_reorder_index': barcode_reorder_index_file,
        })

    return {
        'chunks': chunks,
        'join': {'__mem_gb': JOIN_MEM_GB}
    }

def main(args, outs):
    if args.skip:
        return

    dimred_matrix_file = args.dimred_matrix if args.ordered_dimred_matrix is None else args.ordered_dimred_matrix
    with open(dimred_matrix_file) as fp:
        dimred_matrix = cPickle.load(fp)

    cbc_knn = option(args.cbc_knn, analysis_constants.CBC_KNN)

    batch_to_bc_indices = args.batch_to_bc_indices
    batch_start_idx  = batch_to_bc_indices[args.batch_id][0]
    batch_end_idx = batch_to_bc_indices[args.batch_id][1]
    cur_matrix = dimred_matrix[batch_start_idx:batch_end_idx,:]

    # nearest neighbor pair: stores the nearest neighbors from match_i to match_j
    # key = (batch_i, batch_j), values = set((idx_i, idx_j), ...), the index here is the global index
    batch_nearest_neighbor = defaultdict(set)

    from_idx, to_idx = None, None
    # Batch balanced KNN
    for batch in xrange(len(args.batch_to_bc_indices)):
        if batch == args.batch_id:
            continue

        ref_matrix = dimred_matrix[batch_to_bc_indices[batch][0]:batch_to_bc_indices[batch][1],]
        nn_idx_right = cr_batch_correction.find_knn(cur_matrix, ref_matrix, cbc_knn)

        # convert index (in cur_matrix and ref_matrix) to global index (in dimred_matrix)
        nn_idx_left = np.repeat(np.arange(cur_matrix.shape[0]) + batch_start_idx, cbc_knn)
        nn_idx_right += batch_to_bc_indices[batch][0]

        from_idx = nn_idx_left if from_idx is None else np.concatenate([from_idx, nn_idx_left])
        to_idx = nn_idx_right if to_idx is None else np.concatenate([to_idx, nn_idx_right])

        for i, j in izip(from_idx, to_idx):
            batch_nearest_neighbor[(args.batch_id, batch)].add((i,j))

    outs.batch_nearest_neighbor = martian.make_path('batch_nearest_neighbor.binary')
    with open(outs.batch_nearest_neighbor, 'wb') as fp:
        cr_batch_correction.serialize_batch_nearest_neighbor(fp, batch_nearest_neighbor)

    return


def join(args, outs, chunk_defs, chunk_outs):
    #TODO: clean this up
    #pylint: disable=too-many-locals
    if args.skip:
        return

    chunk_def = chunk_defs[0]
    batch_to_bc_indices = chunk_def.batch_to_bc_indices

    dimred_matrix_file = args.dimred_matrix if chunk_def.ordered_dimred_matrix is None else chunk_def.ordered_dimred_matrix
    with open(dimred_matrix_file) as fp:
        dimred_matrix = cPickle.load(fp)

    with open(chunk_def.idx_to_batch_id) as fp:
        idx_to_batch_id = cPickle.load(fp)

    # batch score before correction
    outs.batch_score_before_correction = cr_batch_correction.batch_effect_score(dimred_matrix, idx_to_batch_id)

    nn_pairs = {}
    for chunk_out in chunk_outs:
        with open(chunk_out.batch_nearest_neighbor, 'rb') as fp:
            batch_nearest_neighbor = cr_batch_correction.deserialize_batch_nearest_neighbor(fp)
        for k,v in batch_nearest_neighbor.iteritems():
            nn_pairs[k] = v

    mutual_nn = {} # mnn between batches
    overlap_percentage = {} # percentage of matching cells (max of batch i and batch j)

    for i in xrange(len(batch_to_bc_indices)):
        batch_i_size = batch_to_bc_indices[i][1] - batch_to_bc_indices[i][0]
        for j in xrange(len(batch_to_bc_indices)):
            if i >= j:
                continue

            if (i,j) not in nn_pairs or (j,i) not in nn_pairs:
                continue

            nn_ij = nn_pairs[(i,j)]
            nn_ji = set([(y,x) for x,y in nn_pairs[(j,i)]])
            mutual_nn[(i,j)] = nn_ij & nn_ji

            batch_j_size = batch_to_bc_indices[j][1] - batch_to_bc_indices[j][0]
            overlap_percentage[(i,j)] = max(
                float(len(set([idx for idx, _ in mutual_nn[(i,j)]]))) / batch_i_size,
                float(len(set([idx for _, idx in mutual_nn[(i,j)]]))) / batch_j_size
            )

    cbc_alpha = option(args.cbc_alpha, analysis_constants.CBC_ALPHA)
    cbc_realign_panorama = option(args.cbc_realign_panorama, analysis_constants.CBC_REALIGN_PANORAMA)
    cbc_sigma = option(args.cbc_sigma, analysis_constants.CBC_SIGMA)

    align_orders = [k for k,v in sorted(overlap_percentage.items(), key=lambda x: x[1], reverse=True) if v > cbc_alpha]

    ## panorama stitch ##
    aligned_dimred_matrix = dimred_matrix
    panoramas = [] # a list of stitched panoramas
    batch_id_to_alignment_count = defaultdict(int)

    for (i,j) in align_orders:
        panorama_idx_i, panorama_idx_j = None, None
        for idx, panorama in enumerate(panoramas):
            if i in panorama:
                assert panorama_idx_i is None, "batch is found in multiple panorama"
                panorama_idx_i = idx

            if j in panorama:
                assert panorama_idx_j is None, "batch is found in multiple panorama"
                panorama_idx_j = idx

        if panorama_idx_i is None:
            panoramas.append(set([i]))
            panorama_idx_i = len(panoramas) - 1

        if panorama_idx_j is None:
            panoramas.append(set([j]))
            panorama_idx_j = len(panoramas) - 1

        # re-align within panorama is enabled
        if cbc_realign_panorama is True:
            batch_id_to_alignment_count[i] += 1
            batch_id_to_alignment_count[j] += 1
            if batch_id_to_alignment_count[i] > 3 and batch_id_to_alignment_count[j] > 3:
                continue
        # re-align within panorama is disable
        else:
            if panorama_idx_i == panorama_idx_j:
                continue

        # use the panorama with the larger shape as reference
        panorama_i_size = sum([batch_to_bc_indices[b][1] - batch_to_bc_indices[b][0] for b in panoramas[panorama_idx_i]])
        panorama_j_size = sum([batch_to_bc_indices[b][1] - batch_to_bc_indices[b][0] for b in panoramas[panorama_idx_j]])
        if panorama_i_size < panorama_j_size:
            panorama_idx_i, panorama_idx_j = panorama_idx_j, panorama_idx_i

        # merge panorama j (current) to panorama i (reference)
        # ordered batch id <-> ordered index range
        batches_in_panorama_j = sorted(list(panoramas[panorama_idx_j]))
        cur_submatrix_idx = np.concatenate([np.arange(batch_to_bc_indices[b][0], batch_to_bc_indices[b][1])
                                            for b in batches_in_panorama_j])

        matches = []
        for ref in panoramas[panorama_idx_i]:
            for cur in panoramas[panorama_idx_j]:
                if ref < cur and (ref, cur) in mutual_nn:
                    matches.extend([(c,r) for r,c in mutual_nn[(ref, cur)]])
                if ref > cur and (cur, ref) in mutual_nn:
                    matches.extend(mutual_nn[(cur, ref)])

        mnn_cur_idx = [a for a, _ in matches]
        mnn_ref_idx = [b for _, b in matches]

        corr_vector = cr_batch_correction.correction_vector(aligned_dimred_matrix, cur_submatrix_idx,
                                                            mnn_cur_idx, mnn_ref_idx, cbc_sigma)
        assert(corr_vector.shape[0] == cur_submatrix_idx.shape[0])
        assert(corr_vector.shape[1] == aligned_dimred_matrix.shape[1])

        base = 0
        for b in batches_in_panorama_j:
            batch_size = batch_to_bc_indices[b][1] - batch_to_bc_indices[b][0]
            aligned_dimred_matrix[batch_to_bc_indices[b][0]:batch_to_bc_indices[b][1], :] += corr_vector[base:(base + batch_size), :]
            base += batch_size

        # merge panoramas and delete panorama j
        if panorama_idx_i != panorama_idx_j:
            panoramas[panorama_idx_i].update(panoramas[panorama_idx_j])
            panoramas.pop(panorama_idx_j)

    # batch score after correction
    outs.batch_score_after_correction = cr_batch_correction.batch_effect_score(aligned_dimred_matrix, idx_to_batch_id)

    if chunk_def.need_reorder_barcode:
        with open(chunk_def.barcode_reorder_index) as fp:
            barcode_reorder_index = cPickle.load(fp)
            aligned_dimred_matrix = aligned_dimred_matrix[barcode_reorder_index]

    # TODO: currently save the aligned matrix (dimred) into pca h5 format, for downstream analysis
    aligned_pca = cr_pca.PCA(aligned_dimred_matrix, np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0))
    pca_map = {aligned_dimred_matrix.shape[1] : aligned_pca}
    filters = tables.Filters(complevel = h5_constants.H5_COMPRESSION_LEVEL)
    with tables.open_file(outs.aligned_pca_h5, 'w', filters = filters) as h5_file:
        cr_pca.save_pca_h5(pca_map, h5_file)

    # load the barcodes and feature info
    with open(args.matrix_barcode_feature_info) as fp:
        bc_feature_info = cPickle.load(fp)
        bcs = bc_feature_info.get('barcodes')
        features = bc_feature_info.get('features')

    cr_pca.save_pca_csv_with_bc_feature(pca_map, bcs, features, outs.aligned_pca_csv)
