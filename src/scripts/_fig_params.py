"""
Figure parameters for kinisi to help make nice plots.

Copyright (c) Andrew R. McCluskey and Benjamin J. Morgan

Distributed under the terms of the MIT License

@author: Andrew R. McCluskey
"""

from matplotlib import rcParams

# blue, orange, green, pink, dark green, grey
colors = ["#0099c8", "#D55E00", "#029E73", "#CC78BC", "#006646", "#949494"]
# colors = ["#0099c8", "#003366", "#99be00", "#006646", "#ff7d00", "#821482"]
FONTSIZE = 8
NEARLY_BLACK = "#161616"
LIGHT_GREY = "#F5F5F5"
WHITE = "#ffffff"
NBINS = 30
ALPHABET = 'abcdefghijklmnopqrstuvwxyz'
CREDIBLE_INTERVALS = [[30.9, 69.1, 0.6], [15.9, 84.1, 0.5], [6.7, 93.3, 0.4], 
                      [2.3, 97.7, 0.3], [0.6, 99.4, 0.2], [0.15, 99.85, 0.1]]

MASTER_FORMATTING = {
    "axes.formatter.limits": (-3, 3),
    "xtick.major.pad": 1,
    "ytick.major.pad": 1,
    "ytick.color": NEARLY_BLACK,
    "xtick.color": NEARLY_BLACK,
    "axes.labelcolor": NEARLY_BLACK,
    "axes.spines.bottom": True,
    "axes.spines.left": True,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 1.0,
    "axes.axisbelow": True,
    "legend.frameon": False,
    'axes.edgecolor': NEARLY_BLACK,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "mathtext.fontset": "custom",
    "font.size": FONTSIZE,
    "font.family": "sans-serif",
    "savefig.bbox": "tight",
    "axes.facecolor": WHITE,
    "axes.labelpad": 2.0,
    "axes.labelsize": FONTSIZE,
    "axes.titlepad": 8,
    "axes.titlesize": FONTSIZE,
    "axes.grid": False,
    "grid.color": WHITE,
    "lines.markersize": 1.0,
    "xtick.major.size": 2.0,
    "xtick.major.width": 1.0,
    "ytick.major.size": 2.0,
    "ytick.major.width": 1.0,
    "lines.scale_dashes": False,
    "xtick.labelsize": FONTSIZE,
    "ytick.labelsize": FONTSIZE,
    "legend.fontsize": FONTSIZE,
    "lines.linewidth": 1,
}
for k, v in MASTER_FORMATTING.items():
    rcParams[k] = v

# 510.0pt. \textwidth
# 246.0pt. \columnwidth

COLUMN_WIDTH = 246.0 / 72.27  # inches
PAGE_WIDTH = 510.0 / 72.27  # inches
GOLDEN_RATIO = (1 + 5**0.5) / 2  # Approx. 1.618

def mid_points(axis):
    """
    Find the mid point in an axis, in axes space.
    """
    x = (axis.get_window_extent().x1 - axis.get_window_extent().x0) / 2 + axis.get_window_extent().x0
    y = (axis.get_window_extent().y1 - axis.get_window_extent().y0) / 2 + axis.get_window_extent().y0
    return x, y