# These define information about the sample required to generate a web summary
#pylint: disable=too-few-public-methods,missing-docstring
import martian
from cellranger import utils

class SampleProperties(object):
    def __init__(self, sample_id, sample_desc, version_from_git=False):
        self.sample_id = sample_id
        self.sample_desc = sample_desc
        if not version_from_git:
            self.version = martian.get_pipelines_version()
        else:
            self.version = utils.get_version()

class CountSampleProperties(SampleProperties):
    """ Various versions of this class are passed around for Count, Aggr, Reanalyze, Spatial
    web summaries, etc."""
    def __init__(self, sample_id, sample_desc, genomes, version_from_git=False):
        super(CountSampleProperties, self).__init__(sample_id, sample_desc,
                                                    version_from_git=version_from_git)
        self.genomes = genomes


class ExtendedCountSampleProperties(CountSampleProperties):
    """ Properties for a count run"""
    def __init__(self, sample_id, sample_desc, genomes, barcode_whitelist,
                 reference_path, version_from_git=False):
        super(ExtendedCountSampleProperties, self).__init__(sample_id, sample_desc, genomes,
                                                            version_from_git=version_from_git)
        self.barcode_whitelist = barcode_whitelist
        self.reference_path = reference_path


class AggrCountSampleProperties(CountSampleProperties):
    """ Properties from an Aggr Run"""
    def __init__(self, sample_id, sample_desc, genomes, agg_batches, version_from_git=False):
        super(AggrCountSampleProperties, self).__init__(sample_id, sample_desc, genomes,
                                                        version_from_git=version_from_git)
        self.agg_batches = agg_batches


class VdjSampleProperties(SampleProperties):
    def __init__(self, sample_id, sample_desc, chain_type, version_from_git=False):
        super(VdjSampleProperties, self).__init__(sample_id, sample_desc,
                                                  version_from_git=version_from_git)
        self.chain_type = chain_type


class SampleDataPaths(object):
    def __init__(self, summary_path=None, barcode_summary_path=None, analysis_path=None,
                 vdj_clonotype_summary_path=None,
                 vdj_barcode_support_path=None,
                 filtered_barcodes_path=None,
                 vdj_cell_barcodes_path=None):
        assert filtered_barcodes_path is None or vdj_cell_barcodes_path is None
        self.summary_path = summary_path
        self.barcode_summary_path = barcode_summary_path
        self.analysis_path = analysis_path
        self.vdj_clonotype_summary_path = vdj_clonotype_summary_path
        self.vdj_barcode_support_path = vdj_barcode_support_path
        self.filtered_barcodes_path = filtered_barcodes_path
        self.vdj_cell_barcodes_path = vdj_cell_barcodes_path
