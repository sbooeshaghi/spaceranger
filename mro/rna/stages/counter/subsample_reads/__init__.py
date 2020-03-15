#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
import cPickle
import json
import cellranger.constants as cr_constants
from cellranger.molecule_counter import MoleculeCounter
import cellranger.subsample as cr_ss
import tenkit.safe_json as tk_safe_json

__MRO__ = """
stage SUBSAMPLE_READS(
    in  h5     molecule_info,
    in  csv    filtered_barcodes,
    in  csv[]  target_features,
    in  string target_mode,
    out json   summary,
    out pickle merged_metrics,
    src py     "../rna/stages/counter/subsample_reads",
) split (
    in  int    chunk_start,
    in  int    chunk_len,
    in  map[]  subsample_info,
    out pickle metrics,
) using (
    mem_gb   = 2,
    volatile = strict,
)
"""


def split(args):
    # construct a range of subsampling depths and write out metadata
    is_targeted = args.target_mode is not None
    subsamplings = cr_ss.construct_all_subsamplings(
        args.molecule_info, args.filtered_barcodes, is_targeted)

    if len(subsamplings) == 0:
        return {"chunks": []}

    # Split the molecule info h5 into equi-RAM chunks
    chunks = []
    tgt_chunk_len = cr_constants.NUM_MOLECULE_INFO_ENTRIES_PER_CHUNK
    mc = MoleculeCounter.open(args.molecule_info, "r")
    for chunk_start, chunk_len in mc.get_chunks(tgt_chunk_len, preserve_boundaries=True):
        chunks.append(
            {
                "chunk_start": chunk_start,
                "chunk_len": chunk_len,
                "subsample_info": subsamplings,
                # The estimate_mem_gb only count the memory usage for the MoleculeCounter object, which is
                # under-estimated the actual memory usage.
                # Based on memory profiling with test case fuzzer_114, actual memory usage is ~4x more
                # than estimate_mem_gb (without cap), here set scale = 6.
                "__mem_gb": MoleculeCounter.estimate_mem_gb(chunk_len, scale=6),
            }
        )
    mc.close()
    return {"chunks": chunks, "join": {"__mem_gb": 6}}


def main(args, outs):
    if args.target_mode is None:
        feature_indices = None
    else:
        feature_indices = set()
        for fn in args.target_features:
            with open(fn, 'r') as f:
                feature_indices.update(int(line.strip()) for line in f)
        feature_indices = sorted(feature_indices)
    data = cr_ss.run_subsampling(
        args.molecule_info,
        args.subsample_info,
        args.filtered_barcodes,
        feature_indices,
        args.chunk_start,
        args.chunk_len)
    with open(outs.metrics, "w") as f:
        cPickle.dump(data, f, protocol=cPickle.HIGHEST_PROTOCOL)


def join(args, outs, chunk_defs, chunk_outs):
    # Merge tallies
    metrics = [chunk.metrics for chunk in chunk_outs]
    data = cr_ss.join_metrics(metrics)

    subsample_info = chunk_defs[0].subsample_info if len(chunk_defs) > 0 else []
    summary = cr_ss.calculate_subsampling_metrics(
        data, args.molecule_info, args.filtered_barcodes, subsample_info, args.target_mode
    )

    with open(outs.summary, "w") as f:
        json.dump(tk_safe_json.json_sanitize(summary), f, indent=4, sort_keys=True)

    with open(outs.merged_metrics, "w") as f:
        cPickle.dump([data, subsample_info], f, protocol=cPickle.HIGHEST_PROTOCOL)
