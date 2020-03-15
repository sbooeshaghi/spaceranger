""" Python versions of the React Components available in the web summary.
The Web Summary data is serialized as JSON and for simple objects we just pass dictionaries around,
 but also use typed classes for some aspects of the Components. """

import json
import numpy as np

from websummary import summarize
import cellranger.websummary.plotly_tools as pltly

class React10X(object): #pylint: disable=too-few-public-methods
    """A base class used to indicate data is related to React Components used in the websummary"""
    def __init__(self):
        pass

    def to_dict_for_json(self):
        """ We currently serialize our Web Summary data with json.dumps
         which only understands a few built-in types.  To support serialization of the web summary
         data, types can implement this method to provide custom serialization of their values.

         Caution: Edits to this returned value will be reflected in the original instance
         """
        return vars(self)


class ReactComponentEncoder(json.JSONEncoder):
    # pylint: disable=too-few-public-methods,arguments-differ
    """An Encoder class to enable json.dumps to interpret React10X objects."""
    # pylint: disable=method-hidden
    def default(self, obj):
        if isinstance(obj, React10X):
            return obj.to_dict_for_json()
        return json.JSONEncoder.default(self, obj)


class ClusteringData(React10X): #pylint: disable=too-few-public-methods
    """POD class for data on a clustering"""
    def __init__(self, key, clustering, data):
        super(ClusteringData, self).__init__()
        self.key = key
        self.name = clustering.description
        self.data = data


class Clusterings(React10X): #pylint: disable=too-few-public-methods
    """Stores similar data for a collection of different possible clusterings such
    as differential expression tables or plots."""
    def __init__(self):
        super(Clusterings, self).__init__()
        self.clusterings = []

    def add_clustering(self, data):
        assert isinstance(data, ClusteringData)
        self.clusterings.append(data)


class SharedCoordinatePlotCollection(Clusterings): #pylint: disable=too-few-public-methods
    """ Class used to store a collection of plots that all share
    the same x,y coordinates and layout/config.  Used for serialization to the WebSummary
    where this type is interpreted and expanded in the JavaScript code. """
    def __init__(self, x, y, config, layout, plt_type, marker):
        """

        :param x: Vector of doubles
        :param y: Vector of doubles
        :param config: dict with plotly config elements
        :param layout: dict with plotly layout elements
        :param plt_type: either `scatter` or `scattergl`
        :param marker: The marker to use in the scatter plot JSON, e.g. {"opacity": 0.9, "size": 4}
        """
        super(SharedCoordinatePlotCollection, self).__init__()
        self.x = x
        self.y = y
        self.config = config
        self.layout = layout
        self.type = plt_type
        self.marker = marker

    def add_clustering(self, key, clustering): #pylint: disable=arguments-differ
        new_clustering = ClusteringData(key, clustering,
                                    self._clustering_to_clusters(clustering))
        super(SharedCoordinatePlotCollection, self).add_clustering(new_clustering)

    def _clustering_to_clusters(self, clustering):
        clustering_labels = clustering.clusters
        num_cells = len(clustering_labels)
        data = []
        for i in range(max(clustering_labels)):
            index = i + 1
            name = "Cluster {}".format(index)
            indices = np.where(clustering_labels == index)[0]
            prop = len(indices) * 1.0 / num_cells
            new_cluster = {
                "name": name,
                "indices": list(indices),
                "type": self.type,
                "mode": "markers",
                "marker": self.marker,
                "text": "{}: {:.1%}".format(name, prop)
            }
            data.append(new_cluster)
        return data

    def to_dict_for_json(self):
        """ We drop the unusesd variables"""
        del self.type
        del self.marker
        return vars(self)


class ClusteringSelector(React10X):
    """Mimics the data requirements and model of the React Component in ClusteringSelector.js"""
    def __init__(self, plot_help_txt, table_help_txt):
        """
        Initialize the data for a ClusteringSelector Component, which displays a top row of plots
        one each on the left and right side, and a bottom row with a table.  Both the plots and the
        table have a header above them.  The plot on the right side can either be a single plot
        that does not change as different clusterings are selected (as in Cell Ranger) or a list of
        plots equal in size to the number of clusters used (as in Space Ranger).

        :param plot_help_txt: dict of data needed for the React component `DynamicHelptext` that
        will be displayed above the plots
        :param table_help_txt: dict of data needed for the React Component `HeaderWithHelp` that
        will be displayed above the table.
        """
        super(ClusteringSelector, self).__init__()
        self._left_plots = None
        self._right_plots = None
        self._tables = None
        self.table_help_txt = table_help_txt
        self.plot_help_txt = plot_help_txt

    def _validate_same_clusters(self):
        """The plots/tables should all be arrays with the same
         cluster represented at each position in them.  """
        to_check = [x for x in [self._left_plots, self._right_plots,
                                self._tables] if x and isinstance(x, Clusterings)]
        if len(to_check) > 1:
            baseline = to_check[0].clusterings
            for alt in to_check[1:]:
                #pylint: disable=deprecated-lambda
                difs = map(lambda x, y: x.key != y.key or x.name != y.name,
                           baseline,
                           alt.clusterings)
                if any(difs):
                    raise ValueError("The clusterings for the plots and tables of a "
                                     "clustering selector must be the same and in order.")

    @property
    def right_plots(self):
        return self._right_plots

    @right_plots.setter
    def right_plots(self, value):
        """The right plots must be a SharedCoordinatePlotCollection """
        assert isinstance(value, SharedCoordinatePlotCollection)
        self._right_plots = value
        self._validate_same_clusters()

    @property
    def left_plots(self):
        return self._left_plots

    @left_plots.setter
    def left_plots(self, value):
        """
        :param value: Either a single plot, or a list of plots.
        :return:
        """
        assert isinstance(value, (dict, SharedCoordinatePlotCollection))
        self._left_plots = value
        self._validate_same_clusters()

    @property
    def tables(self):
        return self.tables

    @tables.setter
    def tables(self, value):
        assert isinstance(value, Clusterings)
        self._tables = value
        self._validate_same_clusters()

    def to_dict_for_json(self):
        l_plots = self._left_plots
        if isinstance(l_plots, SharedCoordinatePlotCollection):
            l_plots = l_plots.to_dict_for_json()
        return {
            "left_plots" : l_plots,
            "right_plots" : self._right_plots,
            "plot_help" : self.plot_help_txt,
            "table_help" : self.table_help_txt,
            "tables" : self._tables
        }



def _change_all_plots(d, plt_func):
    """
    Recursive function to update all elements in a web summary
    :param d: the data to iterate/recurse through
    :param plt_func: A function that takes a dictionary or object representing plot data and
    modifies it
    """
    #pylint: disable=invalid-name
    if isinstance(d, dict) and all(key in d for key in ["data", "layout", "config"]):
        plt_func(d)
    elif isinstance(d, dict):
        for _, v in d.iteritems():
            _change_all_plots(v, plt_func)
    elif isinstance(d, list):
        for v in d:
            _change_all_plots(v, plt_func)
    elif isinstance(d, SharedCoordinatePlotCollection):
        plt_func(d)
    elif isinstance(d, ClusteringSelector):
        _change_all_plots(d.left_plots, plt_func)
        _change_all_plots(d.right_plots, plt_func)


class WebSummaryData(React10X):
    """ Class to store the web summary data as an object instead of as pure JSON.
        This class is a builder and will be invalid after any call to to_dict_for_json
        as this method may modify or change the data to fit conventions of the
        web summary"""
    _SUMMARY_TAB = 'summary_tab'
    _ANALYSIS_TAB = 'analysis_tab'
    # Only used by space ranger
    _ALIGNMENT_TAB = 'alignment_tab'
    _ALARMS = 'alarms'
    _SPATIAL_TAB = 'spatial_tab'
    _CLUSTERING_SELECTOR = 'clustering_selector'
    _VDJ_TAB = 'vdj_analysis_tab'
    _DIAGNOSTICS = 'diagnostics'

    def __init__(self, sample_properties, command, pipeline):
        super(WebSummaryData, self).__init__()
        self._data = {
            WebSummaryData._ALARMS: {WebSummaryData._ALARMS: []},
            'sample': {
                        'id': sample_properties.sample_id,
                        'description': sample_properties.sample_desc,
                        'command': command,
                        'subcommand': pipeline,
                    },
            WebSummaryData._SUMMARY_TAB: {},
            WebSummaryData._ANALYSIS_TAB: {},
            WebSummaryData._ALIGNMENT_TAB: {},
            WebSummaryData._SPATIAL_TAB: {},
            WebSummaryData._VDJ_TAB: {},
            WebSummaryData._DIAGNOSTICS: {},
        }
        self._clustering_selector = None

    def _clear_if_empty(self, value):
        if not self._data[value]:
            del self._data[value]

    def _set_all_plots_font(self, font_name=pltly.DEFAULT_WEB_FONT):
        format_fn = lambda x: pltly.add_font_to_layout(x, font=font_name)
        _change_all_plots(self._data, format_fn)

    # pylint: disable=invalid-name
    def to_dict_for_json(self):
        """ Should only be called once at the end. """
        self._check_valid()
        self._format_clustering_selector()
        to_filter = [WebSummaryData._ALIGNMENT_TAB,
                     WebSummaryData._ANALYSIS_TAB,
                     WebSummaryData._SPATIAL_TAB,
                     WebSummaryData._VDJ_TAB,
                     WebSummaryData._DIAGNOSTICS]
        for k in to_filter:
            self._clear_if_empty(k)
        self._set_all_plots_font()
        to_return = self._data
        # Invalidate future calls by deleting data leading to errors
        del self._data
        return to_return

    def _check_valid(self):
        """ Make sure we are not calling class after formatting and altering it."""
        if not hasattr(self, "_data"):
            raise KeyError("The web summary was accessed after the data was returned and "
                           "the instance invalidated.")


    def _format_clustering_selector(self):
        """If a clustering selector is present, prepare it for serialization to JSON
        and invalidate the reference. """
        if self._clustering_selector:
            self.analysis_tab[WebSummaryData._CLUSTERING_SELECTOR] = self._clustering_selector
            self._clustering_selector = None

    def _get_safe(self, key):
        self._check_valid()
        return self._data.get(key)

    def _set_tab(self, key, value):
        self._check_valid()
        self._data[key] = value

    @property
    def summary_tab(self):
        return self._get_safe(WebSummaryData._SUMMARY_TAB)

    @property
    def alarms(self):
        self._check_valid()
        return self._data[WebSummaryData._ALARMS][WebSummaryData._ALARMS]

    @property
    def analysis_tab(self):
        return self._get_safe(WebSummaryData._ANALYSIS_TAB)

    @property
    def diagnostics(self):
        return self._get_safe(WebSummaryData._DIAGNOSTICS)

    @diagnostics.setter
    def diagnostics(self, value):
        if value is not None:
            self._set_tab(WebSummaryData._DIAGNOSTICS, value)

    @property
    def clustering_selector(self):
        self._check_valid()
        return self._clustering_selector

    @clustering_selector.setter
    def clustering_selector(self, value):
        if value is None: # We can continually set this to None
            self._clustering_selector = value
        else:
            self._check_valid()
            assert isinstance(value, ClusteringSelector)
            self._clustering_selector = value

    @property
    def alignment_tab(self):
        return self._get_safe(WebSummaryData._ALIGNMENT_TAB)

    @alignment_tab.setter
    def alignment_tab(self, value):
        self._set_tab(WebSummaryData._ALIGNMENT_TAB, value)

    @property
    def vdj_tab(self):
        return self._get_safe(WebSummaryData._VDJ_TAB)

    @vdj_tab.setter
    def vdj_tab(self, value):
        self._set_tab(WebSummaryData._VDJ_TAB, value)

    @property
    def spatial_tab(self):
        self._check_valid()
        return self._data[WebSummaryData._SPATIAL_TAB]

def write_html_file(filename, websummary_data):
    """
    Write out a web summary to an HTML file
    :param filename: string of output filename
    :param websummary_data: an instance of WebSummaryData
    :return:
    """
    contents = """<div data-key="summary" data-component="CellRangerSummary">"""
    with open(filename, 'w') as outfile:
        summarize.generate_html_summary({'summary' : websummary_data},
                                        contents, None, outfile, cls=ReactComponentEncoder)
