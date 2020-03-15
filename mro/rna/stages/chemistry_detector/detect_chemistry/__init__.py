#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
# Given a set of FASTQs, infer 3' vs 5' chemistry.
# If a V(D)J reference is also given, try to distinguish 5' unbiased gene exprssion vs VDJ-enriched.

import martian
import json
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
import cellranger.rna.library as rna_library
import cellranger.fastq as cr_fastq
import cellranger.preflight as cr_preflight
import cellranger.io as cr_io
from . import chemistry_helpers as helpers

__MRO__ = """
stage DETECT_CHEMISTRY(
    in  string sample_id,
    in  map[]  sample_def,
    in  path   reference_path,
    in  path   vdj_reference_path,
    in  string chemistry_name_spec,
    in  string[] allowed_chems,
    in  int    r1_length,
    in  int    r2_length,
    out json   summary,
    out string chemistry_type,
    out txt    report,
    out bool   is_antibody_only,
    src py     "stages/chemistry_detector/detect_chemistry",
)
"""

def main(args, outs):
    #TODO: cleanup
    #pylint: disable=too-many-branches,too-many-statements

    outs.is_antibody_only = False

    # Check chemistry restrictions
    if args.allowed_chems is not None and \
       args.chemistry_name_spec not in args.allowed_chems:
        martian.exit("The chemistry name '%s' is not allowed for this pipeline. The allowed values are: %s" % (args.chemistry_name_spec, ', '.join(args.allowed_chems)))

    ## Run preflight checks
    try:
        helpers.run_preflight_checks(args)
    except cr_preflight.PreflightException as e:
        martian.exit(e.msg)

    ## Find the input fastqs
    # 'count' requires library_type to be set. 'vdj' doesn't require a library_type, but only supports VDJ libraries, so let any sample_def entries
    # that don't have library_type set into the detection loop.
    detect_library_types = [rna_library.GENE_EXPRESSION_LIBRARY_TYPE, rna_library.VDJ_LIBRARY_TYPE, None]
    gex_or_vdj_defs = [x for x in args.sample_def if x.get(rna_library.LIBRARY_TYPE) in detect_library_types]
    antibody_defs   = [x for x in args.sample_def if x.get(
        rna_library.LIBRARY_TYPE) == rna_library.ANTIBODY_LIBRARY_TYPE]

    chemistry_name = args.chemistry_name_spec
    report = ''
    metrics = {}

    if len(gex_or_vdj_defs) == 0 and len(antibody_defs) >= 1:
        outs.is_antibody_only = True
        martian.update_progress('Only Antibody Capture library detected...')

        if args.chemistry_name_spec in cr_chem.AUTO_CHEMISTRY_NAMES:
            report += 'Only Antibody Capture library detected.'

            for sd_idx, sd in enumerate(args.sample_def):
                fq_spec = cr_fastq.FastqSpec.from_sample_def(sd)

                # Infer chemistry for each sample index/name (aka fastq group)
                for group, group_spec in fq_spec.get_group_spec_iter():
                    try:
                        chemistry_name = cr_chem.infer_sc3p_chemistry(group_spec)
                    except cr_chem.NoInputFastqsException:
                        # It's okay for a single sample index/name to be absent
                        continue

            if chemistry_name == 'SC3Pv3':
                outs.chemistry_type = chemistry_name
                report += "\nThe chemistry version or sequencing configuration is likely %s" % cr_chem.get_chemistry_description_from_name(chemistry_name)
            else:
                outs.chemistry_type = "SC-FB"
                report += "\nThe chemistry version is either Single Cell 3' v2 or Single Cell 5' v1"

            # Write report file
            martian.log_info(report)
            with open(outs.report, 'w') as f:
                f.write(report + "\n")
            return

        else:
            ## If chem explicitly specified, just check it and finish
            if args.chemistry_name_spec not in cr_chem.AUTO_CHEMISTRY_NAMES or \
               args.chemistry_name_spec == cr_chem.CUSTOM_CHEMISTRY_NAME:
                ok, msg = cr_chem.check_chemistry_arg(args.chemistry_name_spec)
                if not ok:
                    martian.exit(msg)

            # Check that there is a reasonable whitelist hit rate for explicitly set chemistries
            if args.chemistry_name_spec != cr_chem.CUSTOM_CHEMISTRY_NAME:
                for sd_idx, sd in enumerate(args.sample_def):
                    fq_spec = cr_fastq.FastqSpec.from_sample_def(sd)

                    # Check that chemistry correct rate is reasonable.
                    for group, group_spec in fq_spec.get_group_spec_iter():
                        res = cr_chem.check_whitelist_match(args.chemistry_name_spec, group_spec)
                        if res is not None:
                            martian.exit(res)

            # Write empty json
            with open(outs.summary, 'w') as f:
                 json.dump({}, f)

            outs.chemistry_type = args.chemistry_name_spec
            outs.report = None
            return

    ## If chem explicitly specified, just check it and finish
    if args.chemistry_name_spec not in cr_chem.AUTO_CHEMISTRY_NAMES or \
       args.chemistry_name_spec == cr_chem.CUSTOM_CHEMISTRY_NAME:
        ok, msg = cr_chem.check_chemistry_arg(args.chemistry_name_spec)
        if not ok:
            martian.exit(msg)

        # Check that there is a reasonable whitelist hit rate for explicitly set chemistries
        if args.chemistry_name_spec != cr_chem.CUSTOM_CHEMISTRY_NAME:
            for sd_idx, sd in enumerate(args.sample_def):
                fq_spec = cr_fastq.FastqSpec.from_sample_def(sd)

                # Check that chemistry correct rate is reasonable.
                for group, group_spec in fq_spec.get_group_spec_iter():
                    res = cr_chem.check_whitelist_match(args.chemistry_name_spec, group_spec)
                    if res is not None:
                        martian.exit(res)

        # Write empty json
        with open(outs.summary, 'w') as f:
            json.dump({}, f)

        outs.chemistry_type = args.chemistry_name_spec
        outs.report = None
        return

    chunks = helpers.find_fastqs(gex_or_vdj_defs)

    if args.chemistry_name_spec == 'auto':
        (txome_idx, vdj_idx) = helpers.prepare_transcriptome_indexes(args.reference_path, args.vdj_reference_path)

        auto_chemistries = {}
        for (idx, sd) in enumerate(gex_or_vdj_defs):
            chunks = helpers.find_fastqs([sd])

            sd_report = "\nDetect Report -- %s (%s):\n" % (sd["read_path"], sd.get(rna_library.LIBRARY_TYPE))
            chemistry_name, _report, metrics = helpers.infer_sc3p_or_sc5p(chunks, txome_idx, vdj_idx)
            sd_report += _report
            report += sd_report
            auto_chemistries[idx] = chemistry_name
            if not chemistry_name:
                err_msg = ("Unable to detect the chemistry for the following dataset. "
                           "Please validate it and/or specify the chemistry via the --chemistry argument.\n"
                           + sd_report)
                martian.exit(err_msg)


        if len(set(auto_chemistries.itervalues())) > 1:
            c = ', '.join(map(str, set(auto_chemistries.itervalues())))
            s = '\n'.join("  Sample def %d: %s" % (idx, chem) for (idx, chem) in sorted(auto_chemistries.iteritems()))

            any_failed = any(c is None for c in auto_chemistries.itervalues())

            if not any_failed:
                martian.exit("Detected conflicting chemistry types (%s). Please run these data separately. %s" % (c, s))
            else:
                martian.exit("Detected conflicting chemistry types (%s). Please run these data separately and/or specify the chemistry via the --chemistry argument. %s" % (c, s))

        else:
            chemistry_name = auto_chemistries[0]


    # Further refinement:
    #   - Detect the sequencing configuration for SC5P (SC5P-PE vs SC5P-R2)
    #   - Detect the sequencing configuration for SCVDJ (SCVDJ vs SCVDJ-R2)
    #
    # The chemistry/seq-config must be consistent across all sample defs
    if chemistry_name in cr_chem.AUTO_CHEMISTRY_NAMES:
        # Map (sample_def_idx, fastq_group_name) => chemistry_name
        group_chem = {}
        group_exception = {}

        for sd_idx, sd in enumerate(args.sample_def):
            fq_spec = cr_fastq.FastqSpec.from_sample_def(sd)

            # Infer chemistry for each sample index/name (aka fastq group)
            for group, group_spec in fq_spec.get_group_spec_iter():
                try:
                    group_chem[(sd_idx, group)] = cr_chem.infer_chemistry(chemistry_name, group_spec)

                except cr_chem.NoInputFastqsException:
                    # It's okay for a single sample index/name to be absent
                    continue

                except cr_chem.NoChemistryFoundException as e:
                    # It's okay for a single sample index to be unclassifiable
                    group_chem[(sd_idx, group)] = None
                    group_exception[(sd_idx, group)] = e
                    continue

        if len(group_chem) == 0:
            # Could not find any FASTQs
            martian.exit(cr_constants.NO_INPUT_FASTQS_MESSAGE)

        martian.log_info("Detected chemistries:")
        for (i, g) in group_chem.iteritems():
            martian.log_info("%s: %s" % (str(i), str(g)))

        found_chemistries = filter(lambda x: x is not None, group_chem.itervalues())

        # Check for zero chemistry types
        if len(found_chemistries) == 0:
            s = ', '.join("Sample def %d/%s: %s" % (i,g,e) for ((i,g),e) in sorted(group_exception.iteritems()))
            martian.exit("Unable to auto-detect chemistry. %s" % s)

        # Check for multiple chemistry types
        if len(set(found_chemistries)) > 1:
            detected_chemistries = map(str, sorted(list(set(group_chem.itervalues()))))
            c = ', '.join(detected_chemistries)
            s = ', '.join("Sample def %d/%s: %s" % (i,g,v) for ((i,g),v) in sorted(group_chem.iteritems()))

            any_failed = any(c is None for c in group_chem.itervalues())

            if set(detected_chemistries) == set(["SC5P-PE", "SC5P-R2"]):
                msg = "'cellranger count' doesn't support a mixture of 5' paired end (SC5P-PE) and 5' R2 (SC5P-R2) read types. "
                msg += "To process this combination of data, you will need to use 5' single-end mode. Specify '--chemistry SC5P-R2' on the 'cellranger count' command line."
                martian.exit(msg)

            if not any_failed:
                martian.exit("Detected conflicting chemistry types (%s). Please run these data separately. %s" % (c, s))
            else:
                martian.exit("Detected conflicting chemistry types (%s). Please run these data separately and/or specify the chemistry via the --chemistry argument. %s" % (c, s))

        chemistry_name = found_chemistries[0]

        report += "\nThe chemistry version or sequencing configuration is likely %s" % cr_chem.get_chemistry_description_from_name(chemistry_name)

    outs.chemistry_type = chemistry_name

    # Write report file
    martian.log_info(report)

    with open(outs.report, 'w') as f:
        f.write(report + "\n")

    # Write summary JSON
    metrics['chemistry'] = chemistry_name
    with open(outs.summary, 'w') as f:
        json.dump(metrics, f)

    # Check the read-length arguments to make sure they're compatible with the selected chemistry.
    msg = cr_preflight.check_read_lengths_vs_chemistry(chemistry_name,
                                                 args.allowed_chems,
                                                 args.r1_length,
                                                 args.r2_length)
    if msg is not None:
        martian.exit(msg)

def join(args, outs, chunk_defs, chunk_outs):
    cr_io.copy(chunk_outs[0].summary, outs.summary)
    if chunk_outs[0].report is not None:
        cr_io.copy(chunk_outs[0].report, outs.report)
    outs.chemistry_type = chunk_outs[0].chemistry_type
