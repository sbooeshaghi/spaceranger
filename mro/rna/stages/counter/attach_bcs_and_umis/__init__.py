#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
from collections import defaultdict
import itertools
import json
import martian
import os
import sys
import tenkit.log_subprocess as tk_subproc
import tenkit.safe_json as tk_safe_json
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
import cellranger.rna.library as rna_library
import cellranger.reference as cr_reference
import cellranger.utils as cr_utils

__MRO__ = '''
stage ATTACH_BCS_AND_UMIS(
    in  bam[]    genome_inputs,
    in  fastq[]  tags,
    in  path     reference_path,
    in  csv      feature_reference,
    in  int[]    gem_groups,
    in  string[] library_types,
    in  string[] library_ids,
    in  map      chemistry_def,
    in  map      annotation_params,
    in  string   barcode_whitelist,
    in  json     barcode_counts,
    in  json     feature_counts,
    in  float    barcode_confidence_threshold,
    in  int      umi_min_qual_threshold,
    in  string[] bam_comments,
    in  bool     rescue_multimappers,
    in  bool     correct_barcodes,
    in  bool     skip_metrics,
    in  map      skip_translate,
    in  bool     is_antibody_only,
    in  bool     is_pd,
    in  map[]    library_info,
    out bam[]    output,
    out int[]    num_alignments,
    out json     summary,
    out csv      barcodes_detected,
    src py       "../rna/stages/counter/attach_bcs_and_umis",
) split (
    in  bam      chunk_genome_input,
    in  fastq    chunk_tags,
    in  int      gem_group,
    in  string   library_type,
    in  string   library_id,
    in  json     library_info_json,
    in  json     bam_comments_json,
    in  bool     do_skip_translate,
    out bincode chunked_reporter,
) using (
    # No index file is generated for the bam.
    mem_gb   = 2,
    volatile = strict,
)
'''

def split(args):
    # Write BAM comments to json file
    bam_comment_fn = martian.make_path('bam_comments.json')
    with open(bam_comment_fn, 'w') as f:
        json.dump(args.bam_comments, f)

    # Write library info to a file
    libraries_fn = martian.make_path('libraries.json')
    with open(libraries_fn, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(args.library_info), f, indent=4, sort_keys=True)

    chunks = []
    for chunk_genome_input, tags, gem_group, library_type, library_id, in itertools.izip_longest(
            args.genome_inputs, args.tags, args.gem_groups, args.library_types, args.library_ids):

        gem_group_str = str(gem_group)
        if args.skip_translate is not None and \
           gem_group_str in args.skip_translate and \
           library_type in args.skip_translate[gem_group_str]:
            this_skip_translate = args.skip_translate[gem_group_str][library_type]
        elif not args.skip_translate and args.is_antibody_only:
            this_skip_translate = True
            # 3'v3 libraries go through the barcode translation, since the skip_translate flag is correctly
            # set to `False` after the CHECK_BARCODES_COMPATIBILITY stage is run on the metasample.
            # That stage cannot be run for single libraries (including antibody-only libraries), so we do the translation here instead.
            barcode_translate_path = os.path.join(cr_constants.BARCODE_WHITELIST_TRANSLATE_PATH, args.barcode_whitelist + '.txt.gz')
            if os.path.exists(barcode_translate_path): # this is only the case for 3'v3 FB libraries
                this_skip_translate = False
        else:
            this_skip_translate = True

        chunks.append({
            'chunk_genome_input': chunk_genome_input,
            'chunk_tags': tags,
            'gem_group': gem_group,
            'library_type': library_type,
            'library_id': library_id,
            'library_info_json': libraries_fn,
            'bam_comments_json': bam_comment_fn,
            'do_skip_translate': this_skip_translate,
            '__mem_gb': 4,
        })
    join_def = {
        '__mem_gb': 12,
    }
    return {'chunks': chunks, 'join': join_def}

def main(args, outs):
    gene_index_tab = martian.make_path('gene_index.tab')
    convert_pickle_to_rust_index(cr_utils.get_reference_genes_index(args.reference_path), gene_index_tab)

    if args.barcode_whitelist is None:
        barcode_whitelist = 'null'
    elif not os.path.exists(args.barcode_whitelist):
        barcode_whitelist = cr_utils.get_barcode_whitelist_path(args.barcode_whitelist)
    else:
        barcode_whitelist = args.barcode_whitelist

    output = martian.make_path('output.bam')
    outs.output = [output]
    chunked_reporter = martian.make_path('chunked_reporter.bincode')
    chunk_metadata = martian.make_path('chunk_metadata.json')

    cmd = [
        'annotate_reads', 'main',
        args.chunk_genome_input,
        args.chunk_tags,
        output,
        chunked_reporter,
        args.reference_path,
        gene_index_tab,
        args.barcode_counts,
        barcode_whitelist,
        str(args.gem_group),
        chunk_metadata,
        cr_chem.get_strandedness(args.chemistry_def),
        args.feature_counts,
        args.library_type or rna_library.DEFAULT_LIBRARY_TYPE,
        args.library_id,
        args.library_info_json,
        '--bam-comments', args.bam_comments_json,
    ]

    if cr_chem.get_endedness(args.chemistry_def) == cr_constants.FIVE_PRIME:
        cmd.append('--fiveprime')
    else:
        # Redundant with the above flag as of now, but might want to
        # add a martian arg to control this of we want this only in spatial
        # and not in count
        cmd.append('--tso-tag')

    if args.do_skip_translate:
        cmd.append('--skip-translate')
    if args.feature_reference is not None:
        cmd.extend(['--feature-ref', args.feature_reference])
    if args.is_pd:
        cmd.append("--is-pd")

    print >> sys.stderr, 'Running', ' '.join(map(lambda x: "'%s'" % x, cmd))
    tk_subproc.check_call(cmd, cwd=os.getcwd())
    with open(chunk_metadata) as f:
        metadata = json.load(f)
    outs.num_alignments = [metadata['num_alignments']]

def join(args, outs, chunk_defs, chunk_outs):
    outs.output = sum([co.output for co in chunk_outs], [])
    outs.num_alignments = sum([co.num_alignments for co in chunk_outs], [])

    # Write chunk info to a temporary file for the rust code to consume
    chunk_metrics = []
    for chunk_def, chunk_out in zip(chunk_defs, chunk_outs):
        chunk_metrics.append({
            'metrics': chunk_out.chunked_reporter,
            'library_type': chunk_def.library_type,
        })
    with open('chunk_metrics.json', 'w') as f:
        json.dump(chunk_metrics, f)

    cmd = [
        'annotate_reads', 'join',
        'chunk_metrics.json',
        outs.summary,
        outs.barcodes_detected,
    ]
    if args.is_pd:
        cmd.append("--is-pd")
    print >> sys.stderr, 'Running', ' '.join(cmd)
    tk_subproc.check_call(cmd, cwd=os.getcwd())

def convert_pickle_to_rust_index(pickle_path, out_path):
    # NOTE: we could possibly avoid this by using serde-pickle
    # NOTE2: Nope, can't deserialize custom classes apparently.
    gene_index = cr_reference.GeneIndex.load_pickle(pickle_path)
    with open(out_path, 'w') as writer:
        writer.write('\t'.join(["transcript_id", "gene_id", "gene_name", "transcript_length"]))

        # NOTE: It is critical that these genes be in the same order that they are used in
        # FeatureReference construction (i.e., gene_index.genes)
        gene2transcripts = defaultdict(list)
        for (transcript_id, transcript) in gene_index.transcripts.iteritems():
            gene2transcripts[transcript.gene.id].append((transcript_id, transcript))

        for gene in gene_index.genes:
            for (transcript_id, transcript) in gene2transcripts[gene.id]:
                writer.write('\n' + '\t'.join([transcript_id, transcript.gene.id, transcript.gene.name, str(transcript.length)]))
