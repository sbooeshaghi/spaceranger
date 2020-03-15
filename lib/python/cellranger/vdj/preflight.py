#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import os
import socket
from cellranger.preflight import PreflightException
import cellranger.vdj.constants as vdj_constants
import cellranger.vdj.reference as vdj_ref
from tenkit.seq import get_rev_comp


def check_refdata(reference_path, denovo):
    hostname = socket.gethostname()

    if reference_path is None and not denovo:
        raise PreflightException(
            "Must specify --reference unless --denovo is specified.")

    if reference_path is None:
        return

    print "Checking reference_path (%s) on %s..." % (reference_path, hostname)

    required_files = [
        vdj_constants.REFERENCE_FASTA_PATH,
    ]

    for filename in required_files:
        p = os.path.join(reference_path, filename)

        if not os.path.isfile(p):
            raise PreflightException(
                "Your reference does not contain the expected files, including %s, or they are not readable. Please check your reference folder on %s." % (p, hostname))


def check_inner_enrichment_primers(primers_file, reference_path):
    """
    Check that the path is valid, contains only expected characters (ACGT) and targets C-regions
    """

    # 1. Need not specify inner enrichment primers for standard human and mouse VDJ
    if primers_file is None:
        if reference_path is None:
            # If no reference is specified (in denovo mode), make sure primers are specified
            raise PreflightException("You need to specify inner enrichment primers (using --inner-enrichment-primers flag) when a reference is not specified.")
        else:
            # Make sure that we find at least one internally designed primer that targets
            # at least one C-Region in the input reference
            for feat in vdj_ref.get_vdj_feature_iter(reference_path):
                if feat.region_type=='C-REGION':
                    for primer in vdj_constants.VDJ_KNOWN_INNER_PRIMERS:
                        primer_rc = get_rev_comp(primer)
                        if primer_rc in feat.sequence:
                            return

        raise PreflightException("Inner enrichment primers are required for species other than human or mouse for which primers are not provided by 10x Genomics. None of the constant regions in the reference (%s) is targeted by the known primers." % reference_path)

    hostname = socket.gethostname()
    print "Checking enrichment primers (%s) on %s..." % (
        primers_file, hostname)

    # 2. If specified, make sure that the path exists
    if not os.path.isfile(primers_file):
        raise PreflightException(
            "The file specifying inner enrichment primers (%s), does not exists or is not readable. Please check your path on %s." % (primers_file, hostname))

    # 3. Make sure that the file is a newline separated list of ACGT sequences
    inner_primers = []
    with open(primers_file) as f:
        for i, line in enumerate(f.readlines()):
            seq = line.strip()
            if len(seq) == 0:
                raise PreflightException(
                    "Line number %s in the inner enrichment primers file (%s) is empty. You should specify a newline separated list of primers." % (i+1, primers_file))
            for j, base in enumerate(seq):
                if base not in {'A', 'C', 'G', 'T'}:
                    raise PreflightException("Inner enrichment primers file (%s) contain non ACGT characters, which are not supported (Found %s in line %s, character %s). You should specify a newline separated list of primers." % (
                        primers_file, base, i+1, j+1))
            inner_primers.append(seq)

    if not inner_primers: # Empty file
        raise PreflightException("Inner enrichment primers file (%s) contains zero entries. You should specify at least one primer" % primers_file)

    if reference_path:
        # 4. Make sure that every primer targets at least 1 constant region in the reference.
        invalid = []
        for primer in inner_primers:
            # The inner primers are the reverse primers
            primer_rc = get_rev_comp(primer)
            found = False
            for feat in vdj_ref.get_vdj_feature_iter(reference_path):
                if feat.region_type == 'C-REGION':
                    if primer_rc in feat.sequence:
                        found = True
                        break
            if not found:
                invalid.append(primer)

        if invalid:
            invalid_str = ', '.join(invalid)
            raise PreflightException(
                "None of the C-REGIONs in the reference %s is targeted by the following inner enrichment primer(s): %s." % (reference_path, invalid_str))



def check_chain(chain):
    if chain not in vdj_constants.CHAIN_TYPE_SPECS:
        raise PreflightException(
            "Must specify --chain as one of: " + ", ".join(vdj_constants.CHAIN_TYPE_SPECS) + ".")
