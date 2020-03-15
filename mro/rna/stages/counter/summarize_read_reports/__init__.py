#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

# This stage only exists to force earlier VDR of prior stages.

import itertools
import cellranger.chemistry as cr_chem
import cellranger.utils as cr_utils
import cellranger.io as cr_io
import martian

__MRO__ = """
stage SUMMARIZE_READ_REPORTS(
    in  json     extract_reads_summary,
    in  json     barcode_counts,
    in  json     feature_counts,
    in  int[]    gem_groups,
    in  string[] library_types,
    in  string[] library_ids,
    in  string[] read_groups,
    in  map      align,
    in  string[] bam_comments,
    in  fastq[]  read1s,
    in  fastq[]  read2s,
    in  fastq[]  tags,
    in  bool     retain_fastqs,
    in  map      chemistry_def,
    out json     summary,
    out json     barcode_counts,
    out json     feature_counts,
    out int[]    gem_groups,
    out string[] library_types,
    out string[] library_ids,
    out string[] read_groups,
    out map      align,
    out string[] bam_comments,
    out fastq[]  read1s,
    out fastq[]  read2s,
    out fastq[]  tags,
    src py       "../rna/stages/counter/summarize_read_reports",
) split (
    in  fastq    chunk_read1,
    in  fastq    chunk_read2,
    in  fastq    chunk_tags,
) using (
    volatile = strict,
)
"""

def split(args):
    paired_end = cr_chem.is_paired_end(args.chemistry_def)
    if paired_end:
        assert len(args.read1s) == len(args.read2s)
    assert len(args.tags) == len(args.read1s)

    chunks = []

    for read1, read2, tags in itertools.izip_longest(args.read1s, args.read2s, args.tags):
        # Conditionally retain input fastqs.
        # Always retain tag fastq.
        chunks.append({
            'chunk_read1': read1 if args.retain_fastqs else None,
            'chunk_read2': read2 if paired_end and args.retain_fastqs else None,
            'chunk_tags': tags,
        })

    return {'chunks': chunks}

def main(args, outs):
    if args.chunk_read1 is not None:
        # Ensure same extension
        _, extension = cr_utils.splitexts(args.chunk_read1)
        read1 = martian.make_path('read1' + extension)
        cr_io.copy(args.chunk_read1, read1)
        outs.read1s = [read1]
    else:
        outs.read1s = []
    if args.chunk_read2 is not None:
        _, extension = cr_utils.splitexts(args.chunk_read2)
        read2 = martian.make_path('read2' + extension)
        cr_io.copy(args.chunk_read2, read2)
        outs.read2s = [read2]
    else:
        outs.read2s = []
    if args.chunk_tags is not None:
        _, extension = cr_utils.splitexts(args.chunk_tags)
        tags = martian.make_path('tags' + extension)
        cr_io.copy(args.chunk_tags, tags)
        outs.tags = [tags]
    else:
        outs.tags = []

def join(args, outs, chunk_defs, chunk_outs):
    cr_io.copy(args.extract_reads_summary, outs.summary)
    cr_io.copy(args.barcode_counts, outs.barcode_counts)
    cr_io.copy(args.feature_counts, outs.feature_counts)

    outs.gem_groups = args.gem_groups
    outs.library_types = args.library_types
    outs.library_ids = args.library_ids
    outs.read_groups = args.read_groups
    outs.align = args.align
    outs.bam_comments = args.bam_comments

    outs.read1s = sum([co.read1s for co in chunk_outs], [])
    outs.read2s = sum([co.read2s for co in chunk_outs], [])
    outs.tags = sum([co.tags for co in chunk_outs], [])
