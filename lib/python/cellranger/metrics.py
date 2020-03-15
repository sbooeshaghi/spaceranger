class Metrics(object):
    def update(self, other):
        for k,v in other.__dict__.iteritems():
            if v is not None:
                setattr(self, k, v)

class BarcodeFilterResults(Metrics):
    def __init__(self, default_value = None):
        self.filtered_bcs = default_value
        self.filtered_bcs_lb = default_value
        self.filtered_bcs_ub = default_value
        self.max_filtered_bcs = default_value
        self.filtered_bcs_var = default_value
        self.filtered_bcs_cv = default_value

    @staticmethod
    def init_with_constant_call(n_bcs):
        res = BarcodeFilterResults()
        res.filtered_bcs = n_bcs
        res.filtered_bcs_lb = n_bcs
        res.filtered_bcs_ub = n_bcs
        res.max_filtered_bcs = 0
        res.filtered_bcs_var = 0
        res.filtered_bcs_cv = 0
        return res

    def to_dict_with_prefix(self, i):
        return {'gem_group_%d_%s' % (i, key): value for key, value in self.__dict__.iteritems()}
