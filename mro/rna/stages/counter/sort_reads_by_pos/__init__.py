#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import os
import resource
import tenkit.bam as tk_bam
import cellranger.io as cr_io

__MRO__ = """
stage SORT_BY_POS(
    in  bam[]   inputs,
    in  int     num_threads,
    in  int     mem_gb,
    out bam     output,
    out bam.bai index,
    src py      "stages/counter/sort_reads_by_pos",
) split using (
    in  bam     chunk_input,
)
"""

def split(args):
    chunks = []
    for chunk_input in args.inputs:
        chunks.append({
            'chunk_input': chunk_input,
            '__mem_gb': 2,
        })
    join_def = {
        '__threads': args.num_threads,
        '__mem_gb': 4,
    }
    return {'chunks': chunks, 'join': join_def}

def main(args, outs):
    args.coerce_strings()
    bam_prefix = os.path.splitext(outs.output)[0]
    tk_bam.sort(str(args.chunk_input), str(bam_prefix))

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    input_bams = [str(chunk.output) for chunk in chunk_outs]
    tk_bam.merge(outs.output, input_bams, args.__threads)
    outs.index = outs.output + '.bai'
    tk_bam.index(outs.output)
