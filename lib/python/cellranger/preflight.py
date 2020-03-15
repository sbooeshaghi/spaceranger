#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
import os
import json
import socket
import subprocess
import tempfile
import tenkit.log_subprocess as tk_subproc
import tenkit.preflight as tk_preflight
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
import cellranger.rna.library as rna_library
from cellranger.feature_ref import FeatureDefException
from cellranger import csv_utils
import cellranger.rna.feature_ref as rna_feature_ref
from cellranger.spatial.utils import parse_slide_sample_area_id, call_gprreader, read_from_json

ALLOWED_LIBRARY_TYPES = [
    rna_library.GENE_EXPRESSION_LIBRARY_TYPE,
    rna_library.CRISPR_LIBRARY_TYPE,
    rna_library.ANTIBODY_LIBRARY_TYPE,
    rna_library.FEATURETEST_LIBRARY_TYPE,
    rna_library.CUSTOM_LIBRARY_TYPE
]

class PreflightException(Exception):
    def __init__(self, msg):
        super(PreflightException, self).__init__(msg)
        self.msg = msg
    def __str__(self):
        return self.msg

def check(result):
    """ Translate (ok,msg) style from tenkit into an exception """
    ok, msg = result
    if not ok:
        raise PreflightException(msg)

def is_int(s):
    try:
        int(s)
    except ValueError:
        return False
    return True

def check_gex_or_ab_present(sample_def):
    # Skip this check if library type is all VDJ.
    if all(x.get(rna_library.LIBRARY_TYPE) == 'VDJ' for x in sample_def):
        return
    # At least one "Gene Expression" or "Antibody Capture" library is required.
    # Treat an empty library_type as GENE_EXPRESSION. This is done to maintain backward compatibility
    # with historical, v2 samples that didn't have a library_type
    if not any(x.get(rna_library.LIBRARY_TYPE) == rna_library.GENE_EXPRESSION_LIBRARY_TYPE or x.get(rna_library.LIBRARY_TYPE) == rna_library.ANTIBODY_LIBRARY_TYPE or x.get(rna_library.LIBRARY_TYPE) is None for x in sample_def):
        raise PreflightException("You must specify >= 1 input library with either library_type == '%s' or library_type == '%s' to run 'cellranger count'" % (rna_library.GENE_EXPRESSION_LIBRARY_TYPE, rna_library.ANTIBODY_LIBRARY_TYPE))

def check_sample_def(sample_defs, feature_ref=None, pipeline=None):
    hostname = socket.gethostname()

    check(tk_preflight.check_gem_groups(sample_defs))
    if pipeline != cr_constants.PIPELINE_VDJ:
        check_gex_or_ab_present(sample_defs)
    # Check uniqueness of sample_def entries
    sd_entries = sorted([
        (sd.get("read_path"),
         sd.get("sample_names"),
         sd.get("sample_indices"),
         sd.get("lanes")) for sd in sample_defs
    ])

    for i in range(len(sd_entries) - 1):
        if sd_entries[i] == sd_entries[i + 1]:
            msg = "Duplicated entry in the input FASTQ data. Please use a unique combination of fastq path and sample name."
            msg += "\nPath: %s" % sd_entries[i][0]
            msg += "\nNote in demux mode, a unique combination fastq path, sample indices, and lanes is required."
            raise PreflightException(msg)

    print "Checking FASTQ folder..."
    for sample_def in sample_defs:
        read_path = sample_def["read_path"]
        if read_path.strip() == "":
            raise PreflightException("Empty fastq path specifed. Please specify an absolute path.")
        if not read_path.startswith('/'):
            raise PreflightException("Specified FASTQ folder must be an absolute path: %s" % read_path)
        if not os.path.exists(read_path):
            raise PreflightException("On machine: %s, specified FASTQ folder does not exist: %s" % (hostname, read_path))
        if not os.access(read_path, os.X_OK):
            raise PreflightException("On machine: %s, cellranger does not have permission to open FASTQ folder: %s" % (hostname, read_path))
        if not os.listdir(read_path):
            raise PreflightException("Specified FASTQ folder is empty: " + read_path)

        lanes = sample_def["lanes"]
        if lanes is not None:
            for lane in lanes:
                if not is_int(lane):
                    raise PreflightException("Lanes must be a comma-separated list of numbers.")

        check(tk_preflight.check_sample_indices(sample_def))

        if pipeline == cr_constants.PIPELINE_COUNT:
            allowed_public = [x for x in ALLOWED_LIBRARY_TYPES]
            allowed_public.remove(rna_library.FEATURETEST_LIBRARY_TYPE)
            options = ", ".join(("'%s'" % x for x in allowed_public))
            library_type = sample_def.get(rna_library.LIBRARY_TYPE, None)

            # Check for empty library_type
            if library_type == '':
                msg = ("library_type field may not be an empty string."
                    "\nThe 'library_type' field in the libraries csv"
                    " must be one of %s") % options
                raise PreflightException(msg)

            # Check for a valid library_type
            if not (library_type is None or library_type in ALLOWED_LIBRARY_TYPES):
                msg = ("Unknown library_type: '%s'."
                    "\nThe 'library_type' field in the libraries csv"
                    " must be one of %s") % \
                    (library_type, options)
                raise PreflightException(msg)

            # Check that the library_type exists in the feature_ref
            if feature_ref is not None and \
            library_type is not None and \
            library_type != rna_library.GENE_EXPRESSION_LIBRARY_TYPE:

                if not any(x.feature_type == library_type for x in feature_ref.feature_defs):
                    msg = "You declared a library with library_type = '%s', but there are no features declared with that feature_type in the feature reference." % library_type
                    msg += "\nCheck that the 'library_type' field in the libraries csv matches at least 1 entry in the 'feature_type' field in the feature reference csv"
                    raise PreflightException(msg)

        elif pipeline == cr_constants.PIPELINE_VDJ:
            # library type can be missing, or VDJ
            library_type = sample_def.get(rna_library.LIBRARY_TYPE, None)
            if library_type is not None and not (library_type == rna_library.VDJ_LIBRARY_TYPE):
                msg = "You declared a library with library_type = '%s'. For the vdj pipeline, the library_type field in sample_def must be missing or '%s'" % (library_type, rna_library.VDJ_LIBRARY_TYPE)
                raise PreflightException(msg)


def check_refdata(reference_path):
    hostname = socket.gethostname()

    if reference_path is None:
        raise PreflightException("Must specify a transcriptome reference path.")

    print "Checking reference_path (%s) on %s..." % (reference_path, hostname)

    required_files = [
        cr_constants.REFERENCE_METADATA_FILE,
        cr_constants.REFERENCE_FASTA_PATH,
        cr_constants.REFERENCE_GENES_GTF_PATH,
        cr_constants.REFERENCE_GENES_INDEX_PATH,
    ]
    for filename in required_files:
        p = os.path.join(reference_path, filename)
        if not os.path.isfile(p):
            raise PreflightException("Your reference does not contain the expected files, or they are not readable. Please check your reference folder on %s." % hostname)

    for filename in cr_constants.STAR_REQUIRED_FILES:
        p = os.path.join(reference_path, cr_constants.REFERENCE_STAR_PATH, filename)
        if not os.path.exists(p):
            raise PreflightException("Your reference doesn't appear to be indexed. Please run the mkreference tool")

def expand_libraries_csv(csv_path):
    required_cols = set(cr_constants.LIBRARIES_CSV_FIELDS)
    reader = csv_utils.load_csv_filter_comments(csv_path, "libraries", required_cols)

    libraries = []

    for row in reader:
        print row.keys()

        for key in row.iterkeys():
            if key is None or row[key] is None:
                msg = "Invalid libraries CSV file: incorrrect number of columns on line number (after excluding comment lines) %d" % reader.line_num
                raise PreflightException(msg)
            row[key] = row[key].strip()

        if row['sample'].strip() == "":
            raise PreflightException('Empty sample field in libraries csv. Please specify an non-empty sample value for each library.')

        library = {
            "fastq_mode": "ILMN_BCL2FASTQ",
            "gem_group": None,
            "lanes": None,
            "read_path": row["fastqs"],
            "sample_names": [row["sample"]],
            rna_library.LIBRARY_TYPE: row[rna_library.LIBRARY_TYPE],
            "sample_indices": ["any"],
        }

        libraries.append(library)

    return libraries


def check_chemistry(name, custom_def, allowed_chems):
    check(cr_chem.check_chemistry_defs())
    check(cr_chem.check_chemistry_arg(name, allowed_chems))

    if name == cr_chem.CUSTOM_CHEMISTRY_NAME:
        check(cr_chem.check_chemistry_def(custom_def))

def check_read_lengths_vs_chemistry(name, allowed_chems, r1_length, r2_length):
    if name == cr_chem.CUSTOM_CHEMISTRY_NAME:
        return None
    else:
        chem = cr_chem.get_chemistry(name)

    sets = [("R1", "--r1-length", r1_length), ("R2", "--r2-length", r2_length)]

    for (read_type, flag, user_value) in sets:
        if user_value is not None:
            read_def = None
            if cr_chem.get_rna_read_def(chem).read_type == read_type:
                read_def = cr_chem.get_rna_read_def(chem)
            elif cr_chem.get_rna_read2_def(chem).read_type == read_type:
                read_def = cr_chem.get_rna_read2_def(chem)

            min_rna_len = cr_constants.MINIMUM_TRIMMED_READ_LENGTH_PREFLIGHT
            # Check that alignable sequence after bc/umi removal and hard-trimming is long enough
            if read_def is not None and read_def.offset + min_rna_len > user_value:
                msg = "You selected %s %d. In the selected chemistry is '%s', the RNA read on %s starts at position %d, leaving %d bp of alignable sequence after trimming.  At least %d bp of alignable sequence are required." % \
                (flag, user_value, name, read_type, read_def.offset, max(0, user_value - read_def.offset), cr_constants.MINIMUM_TRIMMED_READ_LENGTH_PREFLIGHT)
                msg += "\n"
                msg += "Please check your %s setting" % flag
                return msg
    return None

def check_feature_ref(transcriptome_ref_path, feature_ref_path):

    if not os.path.isfile(feature_ref_path):
        raise PreflightException("Could not find the feature reference file %s" % feature_ref_path)

    if not os.access(feature_ref_path, os.R_OK):
        raise PreflightException("feature reference is not readable, please check file permissions: %s" % feature_ref_path)

    try:
        feature_ref = rna_feature_ref.from_transcriptome_and_csv(
            transcriptome_ref_path,
            feature_ref_path)

        rna_feature_ref.FeatureExtractor(feature_ref)
        return feature_ref

    except FeatureDefException as e:
        raise PreflightException(str(e))

def check_environment():
    check(tk_preflight.check_open_fh())

def record_package_versions():
    for package in cr_constants.PACKAGE_VERSION_CMDS:
        name = package['name']
        cmd = package['cmd']

        version = tk_subproc.check_output(cmd, shell=True)
        print '%s: %s' % (name, version)

def check_read_length(x):
    if x < 1:
        raise PreflightException("Specified read length must be greater than or equal to 1.")

#############################################################################################
def check_feature_preflights(sample_def, feature_reference):
    # If any non "Gene Expression" libraries are present then the feature-ref is required.
    if any((x.get("library_type") != None and x.get("library_type") != rna_library.GENE_EXPRESSION_LIBRARY_TYPE) for x in sample_def):
        if feature_reference is None:
            raise PreflightException("You must specify --feature-ref when using Cell Ranger with feature barcoding libraries.")

def check_sample_info(sample_def, reference_path, feature_reference=None):
    if feature_reference is not None:
        print "Checking feature definition file..."
        feature_ref = check_feature_ref(reference_path, feature_reference)
    else:
        feature_ref = None

    print "Checking sample info..."
    check_sample_def(sample_def, feature_ref, pipeline=cr_constants.PIPELINE_COUNT)

def check_spatial_arguments(tissue_image_path, loupe_alignment_file, gpr_file, slide_serial_capture_area):
    """ Check for the existence of GPR files if given to the pipeline """

    if not os.path.exists(tissue_image_path):
        raise PreflightException("The image file %s does not exist" %tissue_image_path)

    if loupe_alignment_file is not None:
        if not os.path.exists(loupe_alignment_file):
            raise PreflightException("The loupe alignment file %s does not exist" %loupe_alignment_file)
        else:
            try:
                with open(loupe_alignment_file) as file_handle:
                    if True not in (spot.get('tissue') for spot in json.load(file_handle)['oligo']):
                        raise PreflightException("No spots were selected as 'in-tissue' in the Loupe alignment file.")
            except IOError: # file is unreadable
                raise PreflightException("Cannot read Loupe alignment file: %s" %loupe_alignment_file)
            except (ValueError, KeyError): # file is not a json file or is missing required keys
                raise PreflightException("Loupe alignment file is invalid: %s" %loupe_alignment_file)


    if gpr_file is not None:
        if not os.path.exists(gpr_file):
            raise PreflightException("The slide file %s does not exist" %gpr_file)

    _check_slide_capture_area_ids(slide_serial_capture_area,loupe_alignment_file,gpr_file)

def _check_slide_capture_area_ids(slide_serial_capture_area, loupe_alignment_file, gpr_file):
    """ Check that slide_serial_capture_area from commandline arguments matches loupe_alignment_file
    and/or gpr_file if provided."""

    cmdline_slide_id, cmdline_area_id = None, None
    alignment_file_slide_id, alignment_file_area_id = None, None
    gpr_slide_id = None

    if slide_serial_capture_area:
        cmdline_slide_id, cmdline_area_id = parse_slide_sample_area_id(slide_serial_capture_area)

    if loupe_alignment_file is not None:

        with open(loupe_alignment_file) as file_handle:
            alignment_file_data = json.load(file_handle)
            if 'serialNumber' in alignment_file_data:
                alignment_file_slide_id = alignment_file_data['serialNumber']
                if not alignment_file_slide_id: # make sure empty strings are interpretated as None
                    alignment_file_slide_id = None
            if 'area' in alignment_file_data:
                alignment_file_area_id = alignment_file_data['area']
                if not alignment_file_area_id: # make sure empty strings are interpretated as None
                    alignment_file_area_id = None

        if cmdline_slide_id != alignment_file_slide_id or cmdline_area_id != alignment_file_area_id:
            raise PreflightException(
                "You provided {}, but the loupe alignment file {}".format(
                    "--unknown-slide" if not cmdline_slide_id and not cmdline_area_id else \
                    "--slide={} --area={}".format(cmdline_slide_id, cmdline_area_id),
                    "says {}-{}".format(
                        alignment_file_slide_id,
                        alignment_file_area_id
                        ) if alignment_file_slide_id is not None or alignment_file_area_id is not None else \
                    "is empty or malformed"
                ))

    if gpr_file is not None:

        # calling gpr reader
        out_path = tempfile.mkdtemp()
        # set to B1 we just want to read the slide file and get the slide id, we don't care about the area right now
        area_id = 'B1'
        call_gprreader("read", gpr_file, area_id, out_path)
        _, basename_ext = os.path.split(gpr_file)
        basename, _ = os.path.splitext(basename_ext)
        gpr_json = os.path.join(out_path, basename + '_' + area_id + '.json')
        gpr_data = read_from_json(gpr_json)
        os.remove(gpr_json)
        os.rmdir(out_path)
        gpr_slide_id = gpr_data['headers']['barcode']
        if not gpr_slide_id: # make sure empty strings are interpretated as None
            gpr_slide_id = None

        if cmdline_slide_id != gpr_slide_id:
            raise PreflightException(
                "You provided {}, but the slide file {}".format(
                    "--unknown-slide" if not cmdline_slide_id and not cmdline_area_id else \
                    "--slide={}".format(cmdline_slide_id),
                    "says {}".format(gpr_slide_id) if gpr_slide_id is not None else \
                    "is empty or malformed"
                ))

def check_image_open(tissue_image_path, loupe_alignment_file=None):
    """ Check that tiffer can open an image and it's largest dimension is > 2000px """

    tissue_image_checksum = None
    if loupe_alignment_file is not None:
        with open(loupe_alignment_file, 'r') as loupefile:
            loupe_data = json.load(loupefile)
        if 'checksum' in loupe_data:
            tissue_image_checksum = loupe_data['checksum']

    if tissue_image_checksum is not None:
        call = ["tiffer",
                "verify",
                tissue_image_path,
                "--size",
                "2000",
                "--checksum",
                tissue_image_checksum]
    else:
        call = ["tiffer",
                "verify",
                tissue_image_path,
                "--size",
                "2000"]

    unicode_call = [arg.encode('utf-8') for arg in call]
    try:
        tk_subproc.check_output(unicode_call, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError, e:
        raise PreflightException("Tissue image incorrectly formatted:\n%s" % e.output)

def check_common_preflights(reference_path,
                          chemistry_name, custom_chemistry_def, allowed_chems,
                          r1_length, r2_length, check_executables,
                          recovered_cells, force_cells):

    print "Checking reference..."
    check_refdata(reference_path)

    print "Checking chemistry..."
    check_chemistry(chemistry_name, custom_chemistry_def, allowed_chems)

    if r1_length is not None:
        print "Checking read 1 length..."
        check_read_length(r1_length)
    if r2_length is not None:
        print "Checking read 2 length..."
        check_read_length(r2_length)

    # Open file handles limit - per CELLRANGER-824, only check this on the execution machine.
    # We can tell if we're on the execution machine by looking at args.check_executables
    if check_executables:
        print "Checking system environment..."
        check_environment()

    print "Checking optional arguments..."
    if recovered_cells is not None and force_cells is not None:
        raise PreflightException("Cannot specify both --force-cells and --expect-cells in the same run.")
