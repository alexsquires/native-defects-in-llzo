from matplotlib import rc, rcParams, cycler
from collections import OrderedDict

# ---------------------------------------------------
# Color sets
# ---------------------------------------------------
# Standard tableau 20 set
tableau = OrderedDict([
    ('blue',        '#1F77B4'),
    ('orange',      '#FF7F0E'),
    ('green',       '#2CA02C'),
    ('red',         '#D62728'),
    ('purple',      '#9467BD'),
    ('brown',       '#8C564B'),
    ('pink',        '#E377C2'),
    ('grey',        '#7F7F7F'),
    ('yellow',      '#BCBD22'),
    ('turquoise',   '#17BECF'),
    ('lightblue',   '#AEC7E8'),
    ('lightorange', '#FFBB78'),
    ('lightgreen',  '#98DF8A'),
    ('lightred',    '#FF9896'),
    ('lightpurple', '#C5B0D5'),
    ('lightbrown',  '#C49C94'),
    ('lightpink',   '#F7B6D2'),
    ('lightgrey',   '#C7C7C7'),
    ('lightyellow', '#DBDB8D'),
    ('lightturquoise', '#9EDAE5')
])

# Updated Tableau 10 set
# https://www.tableau.com/about/blog/2016/7/colors-upgrade-tableau-10-56782
tableau10_16 = OrderedDict([
    ('blue',   '#4e79a7'),
    ('orange', '#f28e2b'),
    ('green',  '#59a14f'),
    ('red',    '#e15759'),
    ('yellow', '#edc948'),
    ('purple', '#b07aa1'),
    ('cyan',   '#76b7b2'),
    ('pink',   '#ff9da7'),
    ('brown',  '#9c755f'),
    ('grey',   '#bab0ac')
])

fontsize=16

nearly_black = '#161616'
light_grey = '#EEEEEE'
lighter_grey = '#F5F5F5'

formatting = { 'axes.formatter.limits': (-3,3),
               'axes.edgecolor': 'grey',
               'axes.labelcolor': nearly_black,
               'axes.facecolor': 'white',
               'axes.labelpad': 10.0,
               'axes.labelsize': fontsize,
               'axes.titlepad': 30,
               'xtick.major.pad': 7,
               'ytick.major.pad': 7,
               'ytick.color': nearly_black,
               'xtick.color': nearly_black,
               'legend.facecolor': 'white',
               'legend.edgecolor': 'grey',
               'legend.fancybox': False,
               'legend.framealpha': 0,
               'patch.linewidth': 1.0,
               'pdf.fonttype': 42,
               'ps.fonttype': 42,
               'mathtext.fontset': 'custom',
               'font.size': fontsize,
               'text.usetex': True,
               #'font.family': ['sans-serif'], 
               'font.serif': ['Palatino Linotype'],
               'savefig.bbox':'tight',
               'lines.markersize': 5.0,
               'grid.color': 'lightgrey',
               'grid.linestyle': '--',
               'lines.scale_dashes': True }

def set_formatting( formatting ):
    for k, v in formatting.items():
        rcParams[k] = v

def set_colors( color_cycle ):
    rcParams['axes.prop_cycle'] = cycler(color=color_cycle)

set_formatting( formatting )
set_colors( tableau10_16.values() )

# LaTeX packages
rc('text.latex', preamble=r'\usepackage{amssymb}')