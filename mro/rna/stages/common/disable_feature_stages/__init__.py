# Looks at the sample def and determines which feature-counter calls can be disabled
import cellranger.rna.library as rna_library

__MRO__ = """
stage DISABLE_FEATURE_STAGES(
    in  map[]  sample_def,
    out bool  disable_crispr,
    out bool  disable_antibody,
    out bool  disable_internalqc,
    src py     "stages/common/disable_feature_stages",
)
"""

def main(args, outs):
    sample_def = args.sample_def
    library_types = [x.get('library_type') for x in sample_def
                                if x.get('library_type') is not None]

    found_crispr = rna_library.CRISPR_LIBRARY_TYPE in library_types
    found_antibody = rna_library.ANTIBODY_LIBRARY_TYPE in library_types
    found_internalqc = rna_library.FEATURETEST_LIBRARY_TYPE in library_types

    outs.disable_crispr = not(found_crispr)
    outs.disable_antibody = not(found_antibody)
    outs.disable_internalqc = not(found_internalqc)
