#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# 1. Demultiplex reads using the I1 reads, if present. Initially we will detect
#    common sample indicies by looking at the reads. In  the future we will
#    accept a sample sheet
# 2. Put the output FASTQ files in a canonical location
#
# FASTQ file demultiplexer/interleaver for Illumina files.
# Takes a list of fastq files with filenames of the form:
# <prefix>_S0_L001_R1_001.fastq
# where S denotes something, L denotes lane,
# R1/R2/I1/I2 denotes Read 1, Read 2, Index Read 1, Index Read 2.
# The final 001 denotes the file number and is used to split large outputs into
# multiple files.
# If you supply multiple files that differ only by their file number, they will
# be demultiplexed
# in order and the sequences concatenated, dropping the file number index.
#
# All input fastq files must have the same <prefix> string.
#
# The tool will read an index file to determine which are the 'common' barcodes.
# Reads matching the common barcodes will be put into files labelled with the
# barcode sequence. The remaining reads will be put labelled with barcode 'X'.
#
import os
import sys
import itertools
import collections
import json
import subprocess
import numpy
import glob
import gzip
import csv
import collections
import martian
import shutil
from future_builtins import zip
import tenkit.cache as tk_cache
import tenkit.dict_utils as tk_dict
import tenkit.seq as tk_seq
import tenkit.stats as tk_stats
from tenkit.fasta import IlmnFastqFile
from tenkit.constants import DEMULTIPLEX_INVALID_SAMPLE_INDEX, \
                             MAX_INDICES_TO_DEMUX, \
                             SAMPLE_INDEX_MAX_HAMMING_DISTANCE

import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt


__MRO__ = """
stage DEMULTIPLEX(
    in  path     raw_fastq_path,
    in  float    sample_index_error_rate,
    in  bool     interleave,
    in  bool     rc_i2_read,
    in  bool     split_by_tile,
    out path     demultiplexed_fastq_path,
    out json     demultiplex_summary,
    out string[] common_i7_bcs,
    out string[] common_i5_bcs,
    src py       "stages/bcl_processor/demultiplex",
) split (
    in  bool   demultiplex,
    in  string common_i7_bcs,
    in  string common_i5_bcs,
    in  string input_files,
    in  string read_types,
    in  int    chunk_number,
    in  string tile_folder,
)
"""


def join(args, outs, chunk_defs, chunk_outs):
    os.makedirs(outs.demultiplexed_fastq_path)
    
    # Move output file to final location
    for chunk_out in chunk_outs:
        for f in os.listdir(chunk_out.demultiplexed_fastq_path):
            in_file = os.path.join(chunk_out.demultiplexed_fastq_path, f)
            # if this is a tile
            if os.path.isdir(in_file):
                target_dir = os.path.join(outs.demultiplexed_fastq_path, os.path.basename(in_file))
                if not os.path.exists(target_dir):
                    os.makedirs(target_dir)
                for g in os.listdir(in_file):
                    tile_file = os.path.join(in_file, g)
                    shutil.move(tile_file, target_dir)
            else: 
                shutil.move(in_file, outs.demultiplexed_fastq_path)

    # Combine result data
    r = {'num_reads':0, 'num_clusters': 0, 'invalid_count':0, 'sample_index_counts':{}}
    for chunk_out in chunk_outs:
        # We count each end of a paired-end read separately in the summary file.
        summary_counts = json.load(open(chunk_out.demultiplex_summary))
        num_clusters = sum(summary_counts.values())
        num_reads = 2 * num_clusters
        invalid_reads = summary_counts[DEMULTIPLEX_INVALID_SAMPLE_INDEX]
        del summary_counts[DEMULTIPLEX_INVALID_SAMPLE_INDEX]
        summary_counts = {k:2*v for (k,v) in summary_counts.iteritems()}
        r['num_clusters'] += num_clusters
        r['num_reads'] += num_reads
        r['invalid_count'] += invalid_reads
        r['sample_index_counts'] = tk_dict.add_dicts(r['sample_index_counts'], summary_counts, depth=1)
    r['invalid_frac'] = tk_stats.robust_divide(r['invalid_count'], r['num_clusters'])

    json.dump(r, open(outs.demultiplex_summary, 'w'))
    outs.common_i7_bcs = chunk_defs[0].common_i7_bcs
    # Only in case of a dual-index configurations
    if hasattr(chunk_defs[0], 'common_i5_bcs'):
        outs.common_i5_bcs = chunk_defs[0].common_i5_bcs


def main(args, outs):
    if args.demultiplex:
        main_demultiplex_go(args, outs)
    else:
        main_demultiplex(args, outs)

class FastqRow:
    def __init__(self, header, seq, qual):
        self.header = header
        self.seq = seq
        self.qual = qual

    def write(self, f):
        f.write(self.header + "\n")
        f.write(self.seq + "\n")
        f.write("+\n")
        f.write(self.qual + "\n")

class FastqParser:
    def __init__(self, infile, rc=False):
        self.file = infile
        self.rc = rc

    def read_fastq(self):
        if self.file[-2:] == "gz":
            proc = martian.Popen(["gunzip", "--stdout", self.file], stdout=subprocess.PIPE)
            reader = proc.stdout
        else:
            reader = file(self.file, "r")

        while True:
            header = reader.next().strip()
            seq = reader.next().strip()
            reader.next() # incr line
            qual = reader.next().strip()

            if self.rc:
                seq = tk_seq.get_rev_comp(seq)
                qual = qual[::-1]

            yield FastqRow(header, seq, qual)

        reader.close()

class FindCommonBarcodes:
    ''' (int): Sample size constant for the class '''
    SAMPLE_SIZE = int(2e6) # dual-index: we probably need to make this 4e6

    def __init__(self, sample_size=SAMPLE_SIZE):
        self.sample_size = sample_size

    def hamming_distance_between_barcodes(self, bc1, bc2):
        # assert len(bc1) == len(bc2)
        return sum([bc1[i] != bc2[i] for i in range(len(bc1))])
    
    def hamming_distance_to_barcode_list(self, ref_bc, bc_list):
        return map(lambda bc: self.hamming_distance_between_barcodes(ref_bc, bc), bc_list)

    def find_given_hamming_distance_to_barcode_list(self, ref_bc, bc_list, hd=1):
        hd_list = self.hamming_distance_to_barcode_list(ref_bc, bc_list)
        return [index for index, value in enumerate(hd_list) if value == hd]

    def get_index_counts(self, fastq_groups):
        ''' fastq_groups (array(tuples(str))): array with pairs of SIs files '''
        index_counts = {}

        for grp in fastq_groups:
            n = 0
            if (len(grp) == 2):
                i7_fq = grp[0]
                i5_fq = grp[1]
                for i7, i5 in zip(i7_fq.read_fastq(), i5_fq.read_fastq()):
                    if n > self.sample_size:
                        break
                    try:
                        index_counts[(i7.seq, i5.seq)] += 1
                    except KeyError:
                        index_counts[(i7.seq, i5.seq)] = 1
                    n += 1
            elif (len(grp) == 1):
                i7_fq = grp[0]
                universal_i5 = '*'
                for i7 in i7_fq.read_fastq():
                    if n > self.sample_size:
                        break
                    try:
                       index_counts[(i7.seq, universal_i5)] += 1
                    except KeyError:
                        index_counts[(i7.seq, universal_i5)] = 1
                    n += 1

        return sorted(index_counts.items(), key=lambda el: (-el[1], el[0][0]))

    # Look at a bunch of index reads and choose the commonly occuring barcodes.
    # return (common_barcodes, rare_barcodes)
    def pick_common_indexes(self, fastq_groups):
        index_counts = self.get_index_counts(fastq_groups)

        total_counts = sum(v for (k,v) in index_counts)

        group_counter = {}
        for ic in index_counts:
            if ic[0][0] in group_counter:
                group_counter[ic[0][0]] += ic[1]
            else:
                group_counter[ic[0][0]] = ic[1]
        group_counts = sorted(group_counter.items(), key=lambda el: -el[1])

        # dual-index: you need to do this per group
        c = 0
        prop_i7 = 0.90
        i90 = 0
        for i, grp in enumerate(group_counts):
            c += grp[1]
            if c > prop_i7 * total_counts:
                i90 = i
                break

        # number of barcodes that account for 90% of reads
        martian.log_info("CI index is {:d}".format(i90))
        martian.log_info("CI barcode is {}".format(group_counts[i90]))

        # median # of observations of barcodes accounting for the 90%
        num_obs_good_bcs = numpy.median(
            [ grp_count for (bc, grp_count) in group_counts[:(i90+1)] ]
            )
        martian.log_info("Median counts of good barcodes in {:d} reads: {:f}".format(
                total_counts, num_obs_good_bcs)
            )
        # Establish a threshold of .4% of median of high-occurrence SIs
        # or a minimum of 20
        min_obs_bc = max(num_obs_good_bcs / 250, 20)
        martian.log_info(
            "Minimum observed good barcodes counts is {:f}".format(min_obs_bc))

        # only demultiplex a reasonable number of sample indices
        if len(group_counts) > MAX_INDICES_TO_DEMUX:
            min_obs_bc = max(min_obs_bc, float(group_counts[MAX_INDICES_TO_DEMUX][1]))
        martian.log_info(
            "Minimum observed good barcodes counts (after reasonable adjustment) is {:f}".
            format(min_obs_bc))

        # Find barcodes at 1-HD and log them and their counts
        # Start with i7s
        #hd_dict = {}
        #ic_list = map(lambda grp: grp[0], group_counts)
        #for i7 in ic_list[:
        #    hd_list = self.find_given_hamming_distance_to_barcode_list(i7, ic_list, hd=1)
        #    hd_dict[i7] = [ic_list[i] for i in hd_list]
        #hamming_dist_i7_csv_path = martian.make_path('hamming_distance_i7.csv')
        #with open(hamming_dist_i7_csv_path, 'w') as f:
        #    header = ['I7_bc', 'I7_hd_corrections']
        #    csv.writer(f).writerow(header)
        #    for k1 in hd_dict.keys():
        #        csv.writer(f).writerow([k1, str(hd_dict[k1])])
        # Continue with i5s
        i7i5_store = {}
        i7i5_index_pos = {}
        for pos, grp in enumerate(index_counts):
            (i7, i5) = grp[0]
            if i7 not in i7i5_store:
                i7i5_store[i7] = []
                i7i5_index_pos[i7] = []
            i7i5_store[i7].append(i5)
            i7i5_index_pos[i7].append(pos)
        #hd_dict = {}
        #for i7 in ic_list:
        #    ic_list2 = i7i5_store[i7]
        #    for i5 in ic_list2:
        #        hd_list = self.find_given_hamming_distance_to_barcode_list(i5, ic_list2, hd=1)
        #        try:
        #            hd_dict[i5].extend([ic_list2[i] for i in hd_list])
        ##        except KeyError:
        ##            hd_dict[i5] = [ic_list2[i] for i in hd_list]
        #hamming_dist_i5_csv_path = martian.make_path('hamming_distance_i5.csv')
        #with open(hamming_dist_i5_csv_path, 'w') as f:
        #    header = ['I5_bc', 'I5_hd_corrections']
        #    csv.writer(f).writerow(header)
        #    for k2 in hd_dict.keys():
        #        csv.writer(f).writerow([k2, str(hd_dict[k2])])

        major_counts = {}
        group_totals = {}
        i5_counts = {}
        i5_numkeys = {}
        for grp in index_counts:
            (k1, k2) = grp[0]
            cnt = grp[1]
            if k1 not in major_counts:
                major_counts[k1] = cnt
            elif major_counts[k1] < cnt:
                major_counts[k1] = cnt
            if k1 not in group_totals:
                group_totals[k1] = 0
            group_totals[k1] += cnt
            if k2 not in i5_counts:
                i5_counts[k2] = 0
                i5_numkeys[k2] = 0
            i5_counts[k2] += cnt
            i5_numkeys[k2] += 1  
                
        cardinality_csv_path = martian.make_path('dual_index_cardinality_i7.csv')
        with open(cardinality_csv_path, 'w') as csvfile:
            header = ['I7_bc', 'I5_keys', 'total_counts']
            csv.writer(csvfile).writerow(header)
            for k1 in group_totals:
                i5_key_counts = len(i7i5_store[k1])
                csv.writer(csvfile).writerow([k1, i5_key_counts, group_totals[k1]])
        
        cardinality_csv_path = martian.make_path('dual_index_cardinality_i5.csv')
        with open(cardinality_csv_path, 'w') as csvfile:
            header = ['I5_bc', 'I7_keys', 'total_counts']
            csv.writer(csvfile).writerow(header)
            for i5 in i5_counts.keys():
                csv.writer(csvfile).writerow([i5, i5_numkeys[i5], i5_counts[i5]])

        index_counts_csv_path = martian.make_path('dual_index_counts.csv')
        group_counts_path = martian.make_path('group_counts.csv')
        with open(index_counts_csv_path, 'w') as f1, open(group_counts_path, 'w') as f2:
            header = ['I7_bc', 'I5_bc', 'counts', 'proportion']
            csv.writer(f1).writerow(header)
            header = ['I7_bc', 'I5_key_counts', 'I5_major_count', 'I5_total_count', 'is_larger_90pct']
            csv.writer(f2).writerow(header)
            for i7, tot in group_totals.items():
                major_count = major_counts[i7]
                include_count = major_count > 0.90 * tot
                i5_key_counts = 0
                for i, pos_i5 in enumerate(i7i5_index_pos[i7]):
                    grp = index_counts[pos_i5]
                    i5 = grp[0][1]
                    cnt_i7_i5 = grp[1]
                    proportion_i7_i5 = '{:0.5f}'.format(float(cnt_i7_i5) / tot)
                    csv.writer(f1).writerow([i7, i5, str(cnt_i7_i5), proportion_i7_i5])
                csv.writer(f2).writerow([i7, i5_key_counts, major_count, tot, include_count])

        # Histogram of counts for the i7s
        ci_hist_png_path = martian.make_path('dual_index_histogram_ci_i7.png')
        plt.figure()
        d = map(lambda t: t[1], group_counts)
        lim_thres = group_counts[i90+1][1]
        max_h = max([h.get_height() for h in sns.distplot(d, kde=False, color='blue', hist_kws={'log':True}, kde_kws={'cumulative': True}).patches])
        plt.plot([min_obs_bc, min_obs_bc], [0, max_h])
        plt.plot([lim_thres, lim_thres], [0, max_h])
        plt.title('Distribution of i7 counts')
        plt.xlabel("value")
        plt.ylabel("Frequency")
        plt.savefig(ci_hist_png_path)

        ci_hist_png_path = martian.make_path('dual_index_histogram_ci_i5.png')
        plt.figure()
        c = 0
        d = sorted(i5_counts.values(), reverse=True)
        lim_thres = 0
        prop_i5 = 0.90
        for i, i5_d in enumerate(d):
            c += i5_d
            if c > prop_i5 * total_counts:
                lim_thres = d[i-1]
                break
        martian.log_info("lim_thres for i5 index is {}".format(lim_thres))
        max_h = max([h.get_height() for h in sns.distplot(d, kde=False, color='blue', hist_kws={'log':True}, kde_kws={'cumulative': True}).patches])
        plt.plot([lim_thres, lim_thres], [0, max_h])
        plt.title('Distribution of i5 counts')
        plt.xlabel("value")
        plt.ylabel("Frequency")
        plt.savefig(ci_hist_png_path)

        good_bcs = []
        noise_bcs = []
        for i, grp in enumerate(index_counts):
            (i7, i5) = grp[0]
            cnt = grp[1]
            
            # To keep the old behaviour, this should be based in the counts of i7,
            # not the combination (i7, i5)
            total_cnt_i7 = group_totals[i7]
            if total_cnt_i7 > min_obs_bc:
                if cnt > 0.90 * total_cnt_i7:
                    # We have found directly the maximum index
                    good_bcs.append("_".join([i7, i5]))
                    martian.log_info("count of {:d}({:s}) is {:d}".format(i, "_".join((i7, i5)), cnt))
                else:
                    # This is not the maximum index
                    num_i5s_i7 = len(i7i5_store[i7])
                    # If I am not maximum and not the first one, 
                    # how am I doing in terms of HAMMING_DISTANCE
                    # SAMPLE_INDEX_MAX_HAMMING_DISTANCE
                    if self.hamming_distance_between_barcodes(i7i5_store[i7][0], i5) <= \
                        SAMPLE_INDEX_MAX_HAMMING_DISTANCE:
                        # Let's keep these ones even if they are away from
                        # the major i5. 
                        # Correct or no correct? Possibly correct.
                        # Check # of secondary keys in `dual_index_cardinality.csv`
                        good_bcs.append("_".join([i7, i5]))
                        martian.log_info("count of {:d}({:s}) is {:d}".format(i, "_".join((i7, i5)), cnt))
                    else:
                        # Worth checking for error correction(s)?
                        noise_bcs.append("_".join([i7, i5]))
            elif total_cnt_i7 <= min_obs_bc:
                noise_bcs.append("_".join([i7, i5]))

        return (good_bcs, noise_bcs)


# Demultiplex a series of FASTQ iterators.
# The index iterator must be the first iterator
# The filename Vector{String} must all be in the same order of read type.
# Interleave map tells which output to write each of the seq_interator entries to.
def process_fastq_chunk(seq_iters, filenames, no_match_filenames, file_cache,
    _interleave_map, summary_counts, max_reads = -1):

    #out_streams = { k:[ gzip.open(x, open_file_mode) for x in v ] for (k,v) in filenames.items() }
    #no_match_out_streams = [ gzip.open(x, open_file_mode) for x in no_match_filenames ]
    valid_bcs = set(filenames.keys())

    if _interleave_map is None:
        interleave_map = range(len(seq_iters))
    else:
        interleave_map = _interleave_map

    read_iterators = zip(*seq_iters)
    n = 0

    for read_set in read_iterators:
        # Log the counts for each sample index
        bc_seq = read_set[0].seq
        if bc_seq in valid_bcs:
            summary_counts[bc_seq] += 1
        else:
            summary_counts[DEMULTIPLEX_INVALID_SAMPLE_INDEX] += 1

        #target_streams = out_streams.get(bc_seq, no_match_out_streams)
        tfn = filenames.get(bc_seq, no_match_filenames)
        target_streams = [file_cache.get(x) for x in tfn]

        for i in range(len(read_set)):
            target_index = interleave_map[i]
            read_set[i].write(target_streams[target_index])

        n += 1
        if (n%10**5) == 0:
            martian.log_info("Reads processed %i" % n)

        if max_reads > 0 and n >= max_reads:
            break

# Demultiplex a series of FASTQ iterators.
# The index iterator must be the first iterator
# The filename Vector{String} must all be in the same order of read type.
# Interleave map tells which output to write each of the seq_interator entries to.
def process_fastq_chunk_no_demult(seq_iters, filenames, file_cache,
    _interleave_map, summary_counts, max_reads = -1):

    if _interleave_map is None:
        interleave_map = range(len(seq_iters))
    else:
        interleave_map = _interleave_map

    read_iterators = zip(*seq_iters)
    n = 0

    for read_set in read_iterators:
        # Log the counts for each sample index
        summary_counts[DEMULTIPLEX_INVALID_SAMPLE_INDEX] += 1

        target_streams = [file_cache.get(x) for x in filenames]

        for i in range(len(read_set)):
            target_index = interleave_map[i]
            read_set[i].write(target_streams[target_index])

        n += 1
        if (n%10**5) == 0:
            martian.log_info("Reads processed %i" % n)

        if max_reads > 0 and n >= max_reads:
            break




def groupby(f, items):
    groups = collections.defaultdict(list)
    for i in items:
        groups[f(i)].append(i)
    return groups


def _tile_for_fastq_file(args, ilmnFastqFile):
    relpath = os.path.relpath(ilmnFastqFile.filename, args.raw_fastq_path)
    return relpath.split(os.path.sep)[0]

def split(args):
    """
    Split side effect is to sample and then estimate common sample
    indices (mislabeled as 'bcs' in the code) for subsequent processing,
    and pass the `common_bcs` argument to chunks
    """
    # Code supports non-interleaved mode, but we're not currently passing that argument
    #do_interleave = True

    if args.split_by_tile:
        file_glob = os.path.join(args.raw_fastq_path, "Tile*", "Project_*", "*", "*.fastq*")
    else:
        file_glob = os.path.join(args.raw_fastq_path, "Project_*", "*", "*.fastq*")

    files = glob.glob(file_glob)
    if len(files) == 0:
        martian.throw("No FASTQ files were found for this run. Perhaps there was an error in bcl2fastq, or the input data is bad?")

    file_info = [ IlmnFastqFile(x) for x in files ]

    # Some check for consistency of inputs
    #prefixes = set([x.prefix for x in file_info])

    # May need to revisit handling of multiple lanes in the future!
    # if not args.collapse_lanes and len(prefixes) > 1:
    #    martian.log_info("Observed multiple prefixes: %s" % prefixes)
    #    return 1

    if args.split_by_tile:
        file_groups = groupby(lambda x: (_tile_for_fastq_file(args, x), x.s, x.lane, x.group), file_info).items()
        # order by tile/lane/group
        file_groups.sort(key = lambda(k, files): (k[0], k[2], k[3]))
    else:
        file_groups = groupby(lambda x: (x.s, x.lane, x.group), file_info).items()
        # Order the demultiplex by the group filename
        file_groups.sort(key = lambda (k,files): files[0].group)

    num_files_per_group = [len(v) for (k,v) in file_groups]

    if len(set(num_files_per_group)) > 1:
        martian.throw("You are missing or have extra fastq file! Check your input files")

    read_sets = [tuple(sorted(f.read for f in grp_files)) for (grp, grp_files) in file_groups]

    if len(set(read_sets)) > 1:
        martian.throw("You don't have the same set of reads for all read groups! Check your input files!")

    # The list of read_types we are getting, eg. ["R1", "I1", "I2", "R2"]
    read_types = read_sets[0]
    martian.log_info("Read types are: {}".format(read_types))

    # this is specified upstream by the COMPUTE_DEMUX_PARAMS stage,
    # default behavior is to match the index read that has the same
    # length of the default SI length, which is 8 bases; this may
    # not be the case in Neo?
    #
    # dual-index: may need to add a si2_read_type as determined
    # by upstream stage
    index_read = args.si_read_type.split("_")
    martian.log_info("Index read found: {:s}".format(str(index_read)))
    demultiplex =  any(elem in read_types for elem in index_read)
    martian.log_info("Demultiplex is: {}".format(demultiplex))

    #if not (index_read in read_types):
    if not demultiplex:
        martian.log_info("Supplied read types: {:s}".format(str(read_types)))
        martian.log_info("Copying reads with no demultiplexing")
        good_bcs = []

    else:
        # Set up everything we need for demultiplexing
        # Figure out which barcodes are well-represented in the index file
        # We will only demultiplex these ones.
        sort_read_types = []
        sort_read_types.extend(sorted([x for x in read_types if x not in index_read]))
        read_types = index_read + sort_read_types

        # pick out SI files
        # dual-index: will need to pick out i5 files
        # need to sort too as I1 is i7 and I2 is i5
        index_files_for_calibration = [
            sorted([f for f in grp if f.read in index_read], key=lambda f: f.read)
                                        for (k,grp) in file_groups ]

        # pick out the index read files, and then count over them
        # with the pick_common_indexes method
        # Increase number of reads to 4e6
        bcFind = FindCommonBarcodes()
        if args.rc_i2_read:
            # Check with Agora: I think this RC only happens there
            # It does as we sequence ATAC only in NovaSeq
            # WHAT IF agora would ship dual_index? 
            # Check with Vijay, it should not happen as Agora reads CB in I1
            bcFastqs = [map(lambda f: FastqParser(f.filename, rc=(f.read == "I2")), grp) 
                        for grp in index_files_for_calibration]
        else:
            bcFastqs = [map(lambda f: FastqParser(f.filename), grp) 
                        for grp in index_files_for_calibration]
            
        (good_bcs, bad_bcs) = bcFind.pick_common_indexes(bcFastqs)

        martian.log_info("Got {:d} common barcodes".format(len(good_bcs))) # Only if a list()
        martian.log_info("Read types: {:s}".format(str(read_types)))
        # by this point, 'good_bcs' has the sample indices that

        # Log the structure(s) created
        # Save the good stuff for inspection
        goodstuff_json_path = martian.make_path('good_bc_counter.json')
        with open(goodstuff_json_path, 'w') as f:
            json.dump(good_bcs, f)
        # Save noise for inspection
        garbage_json_path = martian.make_path('noise_bc_counter.json')
        with open(garbage_json_path, 'w') as f:
            json.dump(bad_bcs, f)

        # Split the barcodes into the two groups. i5s will be empty if the entire flowcell is single-index
        # py3-compat: force copy list out of map obj
        index_lol = list(map(lambda k: k.split('_'), good_bcs)) 
        i7_good_bcs = [pair[0] for pair in index_lol]
        i5_good_bcs = [pair[1] for pair in index_lol if len(pair) > 1]  # Will be empty if single-index

    chunk_defs = []
    chunk_number = 0
    chunk_len = 1

    for chunk_start in range(0, len(file_groups), chunk_len):
        grps = file_groups[chunk_start:(chunk_start+chunk_len)]

        # each demultiplex chunk will get the list of highest-occurring SIs
        # (good_bcs) determined above, and will limit demultiplexing to
        # those reads
        chunk = {'demultiplex': demultiplex, 'common_i7_bcs': i7_good_bcs, 'common_i5_bcs': i5_good_bcs, 
                'read_types': read_types, 'chunk_number': chunk_number, '__mem_gb': 2}

        chunk['input_files'] = [f.filename for (grp, file_list) in grps for f in file_list]
        if args.split_by_tile:
            tiles = [_tile_for_fastq_file(args, f) for (grp, file_list) in grps for f in file_list]
            tile_set = set(tiles)
            if len(tile_set) > 1:
                martian.throw("File list spans multiple tiles")
            chunk['tile_folder'] = tiles[0]
        else:
            chunk['tile_folder'] = None
        chunk_defs.append(chunk)
        chunk_number += 1

    return {'chunks': chunk_defs}

def main_demultiplex_go(args, outs):
    data = {
        # dual-index: common sample indices may need to change into a dict with
        # i7 to None (all OK) or set of legal i5s
        #'common_sample_indices': args.common_bcs,
        'primary_sample_indices': args.common_i7_bcs,
        'dual_sample_indices': args.common_i5_bcs,
        'file_groups': [],
    }
    file_info = [IlmnFastqFile(x) for x in args.input_files]
    file_groups = groupby(lambda x: (x.s, x.lane, x.group), file_info).items()
    for (_, lane, _), input_files in file_groups:
        # read_type should be all reads from sequencer
        files = {read_type: [f for f in input_files if f.read == read_type][0].filename for read_type in args.read_types}
        data['file_groups'].append({
            'lane': lane,
            'files': files,
        })

    input_json_path = martian.make_path('godemux_input.json')
    with open(input_json_path, 'w') as f:
        json.dump(data, f)

    output_dir = outs.demultiplexed_fastq_path
    if args.split_by_tile:
        output_dir = os.path.join(output_dir, args.tile_folder)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # send common_sis, and file groups to godemux
    subproc_args = ['godemux', input_json_path, output_dir,
                    outs.demultiplex_summary, '--demult-read', args.si_read_type,
                    '--chunk', str(args.chunk_number)]
    if args.rc_i2_read:
        subproc_args += ['--rci2read']
    martian.check_call(subproc_args)

# DEPRECATED
# This code is only here for the case where demultiplex = False
def main_demultiplex(args, outs):
    martian.log_info("INFO: {:s} METHOD IS DEPRECATED!".format(sys._getframe(0).f_code.co_name))
    martian.log_info("You should check that `args.demultiplex` is `True`. Currently `{}`".format(args.demultiplex))

    do_interleave = True
    file_info = [ IlmnFastqFile(x) for x in args.input_files ]
    file_groups = groupby(lambda x: (x.s, x.lane, x.group), file_info).items()

    demultiplex = args.demultiplex
    read_types = args.read_types
    good_bcs = args.common_bcs

    # For no interleaving:
    interleave_map = range(len(args.read_types))
    output_reads = args.read_types

    if not ("R1" in read_types) or not ("R2" in read_types):
        martian.throw("You requested interleaving, but you don't have R1 and R2 read types")

    r1_slot = read_types.index("R1")
    r2_slot = read_types.index("R2")
    interleave_map[r2_slot] = r1_slot
    output_reads = [ read_types[idx] for idx in numpy.unique(interleave_map) ]

    # Create output path
    os.mkdir(outs.demultiplexed_fastq_path)
    output_path = outs.demultiplexed_fastq_path

    # counts of each valid barcode and non-matching barcodes
    summary_counts = { bc:0 for bc in good_bcs }
    summary_counts[DEMULTIPLEX_INVALID_SAMPLE_INDEX] = 0

    with tk_cache.FileHandleCache(open_func=gzip.open) as file_cache:
        # Iterate over the file groups
        for (k, input_files) in file_groups:
            # original path:
            # <path>/<prefix>_S0_L001_R1_001.fastq
            # new path:
            # <outpath>/read-<read_id>_si-xxxxx_lane-<lane>_chunk-<chunk>.fastq
            # input_files should have constant prefix, S, and L
            # sort input_files to match the read_types
            read_to_file_dict = { x.read:x for x in input_files }
            input_files = [ read_to_file_dict[rt] for rt in read_types ]
            output_files = [ read_to_file_dict[rt] for rt in output_reads ]

            def output_file(path, in_file, barcode):
                if do_interleave and in_file.read[0] == "R":
                    read = "RA"
                else:
                    read = in_file.read

                # Chunk over lanes to get some parallelism to speed up alignment
                f = "read-%s_si-%s_lane-%03d-chunk-%03d.fastq.gz" % (read, barcode, in_file.lane, args.chunk_number)
                return os.path.join(path, f)

            if args.rc_i2_read:
                # For NextSeq we need to RC the I2 read
                input_iters = [ FastqParser(f.filename, rc=(f.read == "I2")).read_fastq() for f in input_files ]
            else:
                input_iters = [ FastqParser(f.filename).read_fastq() for f in input_files ]

            martian.log_info("Demultiplexing from: %s" % input_files[0].filename)

            if demultiplex:
                bc_files = { bc: [output_file(output_path, f, bc) for f in output_files] for bc in good_bcs }
                err_files = [ output_file(output_path, f, "X") for f in output_files ]
                process_fastq_chunk(input_iters, bc_files, err_files, file_cache, interleave_map, summary_counts)

            else:
                out_files = [ output_file(output_path, f, 'X') for f in output_files ]
                process_fastq_chunk_no_demult(input_iters, out_files, file_cache, interleave_map, summary_counts)

        output_files = file_cache.have_opened

    # Write out the summary counts to JSON
    with open(outs.demultiplex_summary, "w") as f:
        json.dump(summary_counts, f)
