"""
Parses a target feature file -- this file should be a line by line file with gene names or IDs.
There will be one of these per sample, in the sample_def.

To make this compatible with Feature Barcoding, this stage will add in any features to the
target set that were defined in the optional input feature_reference file.

"""

import martian
import cellranger.io as cr_io
import cellranger.utils as cr_utils
import cellranger.reference as cr_reference
import cellranger.sample_def as cr_sample_def
import cellranger.rna.feature_ref as rna_feature_ref

__MRO__ = """
stage PARSE_TARGET_FEATURES(
    in  map[] sample_def,
    in  path  reference_path,
    in  csv   feature_reference,
    in  bool  is_antibody_only,
    out csv[] target_features,
    out csv[] off_target_features,
    out csv   target_gene_ids,
    out bool  disable_targeted,
    src py    "stages/common/parse_target_features",
) using (
    mem_gb = 2,
)
"""

def main(args, outs):
    outs.target_features = []
    outs.off_target_features = []

    ref_idx = cr_utils.get_reference_genes_index(args.reference_path)
    ref_gene_index = cr_reference.GeneIndex.load_pickle(ref_idx)
    gene_name_map = {gene.name: i for i, gene in enumerate(ref_gene_index.genes)}

    # keep track of all ints in order to find off-target genes
    all_gene_ints = {ref_gene_index.gene_id_to_int(x) for x in ref_gene_index.get_gene_ids()}
    # we keep track of the union of target features for this analysis
    # this union will be what gets sliced out of the count matrix
    union_target_gene_ids = set()

    # if a feature reference is used, bring in all of those IDs and
    # concatenate them to the target set
    if args.feature_reference is not None:
        csv_feature_defs, _ = rna_feature_ref.parse_feature_def_file(
            args.feature_reference,
            index_offset=len(ref_gene_index.genes))
        feature_ref_ids = {x.id for x in csv_feature_defs}
        feature_ref_indices = [x.index for x in csv_feature_defs]
        union_target_gene_ids.update(feature_ref_ids)
        all_gene_ints.update(feature_ref_indices)
    else:
        feature_ref_ids = set()
        feature_ref_indices = []

    # FIXME(joey.arthur) handle this more generically --  pipeline should know
    #     that gene expression features are the only ones currently being targeted
    if args.is_antibody_only is True:
        # stick with just the feature reference and stop here
        # this is necessary because target_gene_ids will be
        # used later to subset the feature-barcode matrix
        cr_io.write_target_features_csv(outs.target_gene_ids, feature_ref_ids,
                                        header='Gene')
        outs.disable_targeted = True
        return

    # do we have any targeted analyses here?
    # If not, then just write every gene ID to a file so
    # that SC_RNA_ANALYZER runs on the full gene set
    outs.disable_targeted = all(cr_sample_def.get_target_set(sd) is None
                                for sd in args.sample_def)
    if outs.disable_targeted:
        all_gene_ids = set(ref_gene_index.get_gene_ids()) | feature_ref_ids
        # this is used to hook in to the existing gene subsetting functionality of the Matrix object
        cr_io.write_target_features_csv(
            outs.target_gene_ids, all_gene_ids, header='Gene')
        return

    for i, sample_def in enumerate(args.sample_def):
        target_genes_filename = cr_sample_def.get_target_set(sample_def)
        if target_genes_filename is None:  # no targeting for this library
            target_features = all_gene_ints
            off_target_features = []
            # TODO: should it be possible to perform analysis only on the target set
            #   when a non-targeted library is combined with a targeted one?
            union_target_gene_ids.update(set(ref_gene_index.get_gene_ids()))
        else:  # targeted library
            target_features = feature_ref_indices[:]
            with open(target_genes_filename, 'r') as fh:
                for gene_id in fh:
                    gene_id = gene_id.rstrip()
                    # we support gene_id, gene_id without suffix, and gene_name
                    gene_ints = [ref_gene_index.gene_id_to_int(gene_id),
                                 ref_gene_index.gene_id_to_int(gene_id.split('.')[0]),
                                 ref_gene_index.gene_id_to_int(gene_name_map.get(gene_id, None))]
                    valid_gene_ints = [x for x in gene_ints if x is not None]
                    if len(valid_gene_ints) == 0:
                        martian.exit('Gene {} not seen in reference'.format(gene_id))
                    j = valid_gene_ints[0]
                    target_features.append(j)
                    union_target_gene_ids.add(ref_gene_index.genes[j].id)
                target_features = sorted(target_features)

            off_target_features = all_gene_ints - set(target_features)
            assert len(off_target_features) != len(all_gene_ints)

        on_target = martian.make_path('{}.on_target.csv'.format(i))
        off_target = martian.make_path('{}.off_target.csv'.format(i))
        cr_io.write_target_features_csv(on_target, target_features)
        cr_io.write_target_features_csv(off_target, off_target_features)
        outs.target_features.append(on_target)
        outs.off_target_features.append(off_target)

    # this is used to hook in to the existing gene subsetting functionality of the Matrix object
    cr_io.write_target_features_csv(outs.target_gene_ids, union_target_gene_ids, header='Gene')
