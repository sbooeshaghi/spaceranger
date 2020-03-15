#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
import os
from collections import namedtuple, OrderedDict
import h5py
import cellranger.io as cr_io
import cellranger.h5_constants as h5_constants

FeatureDef = namedtuple(
    'FeatureDef', ['index', 'id', 'name', 'feature_type', 'tags'])

# Required HDF5 datasets
REQUIRED_DATASETS = ['id', 'name', 'feature_type']

# These feature tag keys are reserved for internal use.
GENOME_FEATURE_TAG = 'genome'
DEFAULT_FEATURE_TAGS = [GENOME_FEATURE_TAG]
RESERVED_TAGS = DEFAULT_FEATURE_TAGS

class FeatureDefException(Exception):
    pass


class TargetFeatureSets(object):
    """
    Each target set should be a list of integers with a dictionary key representing the name.
    This name will also exist in the library_info.
    """
    def __init__(self, target_features):
        self.target_feature_sets = dict()
        if target_features is not None:
            for name, s in target_features.iteritems():
                assert all(isinstance(x, int) for x in s)
                self.target_feature_sets[name] = set(s)

    def __eq__(self, other):
        if not isinstance(other, TargetFeatureSets):
            return False
        return self.target_feature_sets == other.target_feature_sets

    @property
    def names(self):
        return self.target_feature_sets.keys()

    def iteritems(self):
        for key, val in self.target_feature_sets.iteritems():
            yield key, val

    def get(self, key):
        return self.target_feature_sets.get(key)

    def update(self, other):
        self.target_feature_sets.update(other)

    def get_all_target_features(self):
        return set.union(*self.target_feature_sets.values())


class FeatureReference(object):
    '''Store a list of features (genes, antibodies, etc).'''

    def __init__(self, feature_defs, all_tag_keys, target_features=None):
        '''Create a FeatureReference.

        Args:
            feature_defs (list of FeatureDef): All feature definitions.
            all_tag_keys (list of str): All optional tag keys.
            target_features (dictionary of list of int): Optional target set(s). Each target set
            should be a list of integers with a dictionary key representing the name. This name
            will also exist in the library_info.
        '''
        self.feature_defs = feature_defs
        self.all_tag_keys = all_tag_keys

        if target_features is not None:
            self.target_features = TargetFeatureSets(target_features)
        else:
            self.target_features = None

        # Assert uniqueness of feature IDs
        id_map = {}
        for fdef in self.feature_defs:
            if fdef.id in id_map:
                this_fd_str = 'ID: %s; name: %s; type: %s' % (
                    fdef.id, fdef.name, fdef.feature_type)
                seen_fd_str = 'ID: %s; name: %s; type: %s' % (
                    id_map[fdef.id].id, id_map[fdef.id].name, id_map[fdef.id].feature_type)
                raise FeatureDefException('Found two feature definitions with the same ID: '
                                          '(%s) and (%s). All feature IDs must be distinct.'
                                          % (this_fd_str, seen_fd_str))
            id_map[fdef.id] = fdef

        self.id_map = id_map

    def get_count_of_feature_type(self, feature_type):
        ''' Count how many features in the matrix are of a given type. (e.g. "Gene Expression")'''
        total = 0
        for feature in self.feature_defs:
            if feature.feature_type == feature_type:
                total += 1
        return total

    def __eq__(self, other):
        return self.feature_defs == other.feature_defs and \
               self.all_tag_keys == other.all_tag_keys and \
               (self.target_features is None) == (other.target_features is None) and \
               self.target_features == other.target_features

    def __ne__(self, other):
        return not self == other

    @staticmethod
    def addtags(feature_ref, new_tags, new_labels=None):
        """Add new tags and corresponding labels to existing feature_ref
           If new labels are None, empty strings are supplied by default

        Args:
            feature_ref: a FeatureReference instance
            new_tags: a list of new tags
            new_labels: per feature list of label values corresponding to the new tags
        """
        assert len(new_tags) > 0
        for tag in new_tags:
            if tag in feature_ref.all_tag_keys:
                raise ValueError("tag {} is already present in feature_ref")

        if new_labels is not None:
            if len(feature_ref.feature_defs) != len(new_labels):
                raise ValueError("number of labels does not match number of features "
                                 "in feature_ref")
            for labels in new_labels:
                assert len(labels) == len(new_tags)
            use_labels = new_labels
        else:
            # initialize to empty
            use_labels = [[""] * len(new_tags) for _ in xrange(len(feature_ref.feature_defs))]
        assert len(feature_ref.feature_defs) == len(use_labels)

        augmented_features = []
        for fd, newvals in zip(feature_ref.feature_defs, use_labels):
            A = {a: b for a, b in zip(new_tags, newvals)}
            A.update(fd.tags)
            augmented_features.append(FeatureDef(index=fd.index,
                                                 id=fd.id,
                                                 name=fd.name,
                                                 feature_type=fd.feature_type,
                                                 tags=A))

        return FeatureReference(feature_defs=augmented_features,
                                all_tag_keys=feature_ref.all_tag_keys + new_tags,
                                target_features=feature_ref.target_features)

    @staticmethod
    def join(feature_ref1, feature_ref2):
        """Concatenate two feature references, requires unique ids and identical tags"""
        assert feature_ref1.all_tag_keys == feature_ref2.all_tag_keys
        feature_defs1 = feature_ref1.feature_defs
        feature_defs2 = feature_ref2.feature_defs

        if feature_ref1.target_features is None:
            combined_target_features = feature_ref2.target_features
        elif feature_ref2.target_features is None:
            combined_target_features = feature_ref1.target_features
        else:
            combined_target_features = feature_ref1.target_features
            # if feature_ref2 has the same keys, they will be over-written
            combined_target_features.update(feature_ref2.target_features)
        return FeatureReference(feature_defs=feature_defs1 + feature_defs2,
                                all_tag_keys=feature_ref1.all_tag_keys,
                                target_features=combined_target_features)

    @classmethod
    def empty(cls):
        return cls(feature_defs=[], all_tag_keys=[], target_features=None)

    def get_num_features(self):
        return len(self.feature_defs)

    def get_feature_ids_by_type(self, feature_type):
        """ Return a list of feature ids of a particular feature type (e.g. "Gene Expression")"""
        return [f.id for f in self.feature_defs if f.feature_type==feature_type]

    def get_indices_for_type(self, feature_type):
        indices = []
        for feature in self.feature_defs:
            if feature.feature_type == feature_type:
                indices.append(feature.index)
        return indices

    def get_genomes(self, feature_type=None):
        """Get sorted list of genomes. Empty string is for reverse compatibility.
        A specific feature type can optionally be specified.
        """
        genomes = set(f.tags.get(GENOME_FEATURE_TAG, '') for f in self.feature_defs
                      if (feature_type is None or f.feature_type == feature_type))
        return sorted(genomes)

    def get_target_features(self):
        """Gets all on-target features in this FeatureReference, collapsed across samples"""
        return set.union(*self.target_features.target_feature_sets.values())

    def get_off_target_features(self):
        """Gets all off-target features in this FeatureReference, collapsed across samples"""
        all_features = set(range(self.get_num_features()))
        return all_features - self.get_target_features()

    def to_hdf5(self, group):
        """Write to an HDF5 group."""

        # Write required datasets
        for col in REQUIRED_DATASETS:
            data = [getattr(f, col) for f in self.feature_defs]
            cr_io.create_hdf5_string_dataset(group, col, data, compression=True)

        # Write tag datasets
        for col in self.all_tag_keys:
            # Serialize missing data as empty unicode string
            data = [f.tags.get(col, '') for f in self.feature_defs]
            cr_io.create_hdf5_string_dataset(group, col, data, compression=True)

        # Write target_features as a new sub-group
        if self.target_features is not None:
            target_set_group = group.create_group(h5_constants.H5_TARGET_SET_ATTR)
            for key, val in self.target_features.iteritems():
                cr_io.create_hdf5_string_dataset(target_set_group, key, map(str, sorted(val)), compression=True)

        # Record names of all tag columns
        cr_io.create_hdf5_string_dataset(
            group, '_all_tag_keys', self.all_tag_keys)

    @classmethod
    def from_hdf5(cls, group):
        """Load from an HDF5 group.
        Args:
            group (h5py.Dataset): Group to load from.
        Returns:
            feature_ref (FeatureReference): New object.
        """
        # FIXME: ordering may not be guaranteed in python3
        def load_str(node):
            if node.shape is None:
                return []
            if node.dtype.char == 'S':
                return cr_io.read_hdf5_string_dataset(node)
            else:
                return node[:]

        def format_path(path):
            # strip off leading slash and convert to str (could be unicode)
            return str(path[1:])

        def h5py_dataset_iterator(group, prefix=''):
            for key in group:
                item = group[key]
                path = '/'.join([prefix, key])
                if isinstance(item, h5py.Dataset):
                    yield format_path(path), load_str(item)
                elif isinstance(item, h5py.Group):
                    # TODO: in python 3.3 this can be a yield from statement
                    for subpath, subitem in h5py_dataset_iterator(item, path):
                        yield subpath, subitem

        data = OrderedDict(h5py_dataset_iterator(group))

        # Load Tag Keys
        all_tag_keys = data['_all_tag_keys']

        # Load Target Sets, if they exist
        target_features = dict()
        for key, val in data.iteritems():
            if h5_constants.H5_TARGET_SET_ATTR in key:
                key = os.path.basename(key)
                target_features[key] = map(int, val)
        if len(target_features) == 0:
            target_features = None

        # Build FeatureDefs
        feature_defs = []
        num_features = len(data[REQUIRED_DATASETS[0]])

        for i in xrange(num_features):
            fdef_dict = {col: data[col][i] for col in REQUIRED_DATASETS}
            tags = {k: data[k][i] for k in all_tag_keys if len(data[k][i]) > 0}
            fdef_dict['index'] = i
            fdef_dict['tags'] = tags
            feature_defs.append(FeatureDef(**fdef_dict))

        return cls(feature_defs, all_tag_keys=all_tag_keys, target_features=target_features)
