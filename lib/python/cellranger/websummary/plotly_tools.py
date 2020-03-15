"""Shared ploting code for the web summary."""
import copy

BUTTON_RESET_SCALED_2D = "resetScale2d"
BUTTON_TO_IMAGE = "toImage"
TO_IMAGE_BUTTON_OPTIONS = "toImageButtonOptions"
MODE_BAR_BUTTONS = "modeBarButtons"
CONFIG = "config"
SHOW_AXIS_DRAG_HANDLES = "showAxisDragHandles"

PLOT_CONFIG = {
    "staticPlot": False,
    "displayModeBar": True,
    MODE_BAR_BUTTONS: [[BUTTON_TO_IMAGE, BUTTON_RESET_SCALED_2D]],
    SHOW_AXIS_DRAG_HANDLES: True,
    TO_IMAGE_BUTTON_OPTIONS: {"width": 730, "height": 650},      # ensure that this is square
    "scrollZoom": False,
}

# Make a spatial config that is similar but with no zoom or reset axis button
SPATIAL_PLOT_CONFIG = copy.deepcopy(PLOT_CONFIG)
SPATIAL_PLOT_CONFIG[MODE_BAR_BUTTONS] = [[]]
del SPATIAL_PLOT_CONFIG[TO_IMAGE_BUTTON_OPTIONS]
SPATIAL_PLOT_CONFIG[SHOW_AXIS_DRAG_HANDLES] = False

DEFAULT_WEB_FONT = "DIN Next LT Pro"

def add_font_to_layout(chart, font=DEFAULT_WEB_FONT):
    """ Updates a charts `layout` dictionary to use a specific font."""
    fnt = {"family": font}
    if isinstance(chart, dict):
        layout = chart['layout']
    elif hasattr(chart, "layout"):
        layout = chart.layout
    else:
        raise Exception("Must pass an object or dictionary to update layout font.")
    font_dict = layout.get('font', {})
    font_dict.update(fnt)
    layout['font'] = font_dict
