from matplotlib import rc, rcParams
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

fontsize=20
nearly_black = '#161616'
light_grey = '#EEEEEE'
lighter_grey = '#F5F5F5'
white = '#ffffff'
grey = '#7F7F7F'

master_formatting = {'axes.formatter.limits': (-3,3),
                          'xtick.major.pad': 7,
                          'ytick.major.pad': 7,
                          'ytick.color': nearly_black,
                          'xtick.color': nearly_black,
                          'axes.labelcolor': nearly_black,
                          'axes.spines.bottom': False,
                          'axes.spines.left': False,
                          'axes.spines.right': False,
                          'axes.spines.top': False,
                          'axes.axisbelow': True,
                          'legend.frameon': False,
                          'pdf.fonttype': 42,
                          'ps.fonttype': 42,
                          'mathtext.fontset': 'custom',
                          'font.size': fontsize,
                          'font.family' : 'serif',
                          #'font.serif' : 'Source Serif Pro',
                          'text.usetex': True,
                          'savefig.bbox':'tight',
                          'axes.facecolor': lighter_grey,
                          'axes.labelpad': 10.0,
                          'axes.labelsize': fontsize * 0.8,
                          'axes.titlepad': 30,
                          'axes.titlesize': fontsize,
                          'axes.grid': False,
                          #'grid.color': white,
                          'lines.markersize': 7.0,
                          'lines.scale_dashes': False,
                          'xtick.labelsize': fontsize * 0.8,
                          'ytick.labelsize': fontsize * 0.8,
                          'legend.fontsize': fontsize * 0.8}

stab_region_formatting = {'axes.formatter.limits': (-3,3),
                          'xtick.major.pad': 7,
                          'ytick.major.pad': 7,
                          'ytick.color': nearly_black,
                          'xtick.color': nearly_black,
                          'axes.labelcolor': nearly_black,
                          'axes.spines.top': True,
                          'axes.spines.right': True,
                          'pdf.fonttype': 42,
                          'ps.fonttype': 42,
                          'mathtext.fontset': 'custom',
                          'font.size': fontsize,
                          'text.usetex' : True,
                          'text.usetex': True,
                          'savefig.bbox':'tight',
                          'axes.facecolor': white,
                          'axes.labelpad': 10.0,
                          'axes.labelsize': fontsize,
                          'axes.titlepad': 30,
                          'lines.markersize': 5.0,
                          'lines.scale_dashes': False,
                          'xtick.labelsize': fontsize * 0.8,
                          'ytick.labelsize': fontsize * 0.8}

transitionfig_formatting = {'axes.formatter.limits': (-3,3),
                          'xtick.major.pad': 7,
                          'ytick.major.pad': 7,
                          'ytick.color': nearly_black,
                          'xtick.color': nearly_black,
                          'axes.labelcolor': nearly_black,
			              'axes.spines.bottom': True,
                          'axes.spines.left': True,
                          'axes.spines.right': True,
                          'axes.spines.top': True,
                          'pdf.fonttype': 42,
                          'ps.fonttype': 42,
                          'mathtext.fontset': 'custom',
                          'font.size': fontsize,
                          'text.usetex' : True,
                          'text.usetex': True,
                          'savefig.bbox':'tight',
                          'axes.facecolor': white,
                          'axes.labelpad': 10.0,
                          'axes.labelsize': fontsize,
                          'axes.titlepad': 30,
                          'lines.markersize': 7.0,
                          'lines.scale_dashes': False,
                          'xtick.labelsize': fontsize * 0.8,
                          'ytick.labelsize': fontsize * 0.8}

space_charge_profile   = {'axes.formatter.limits': (-3,3),
                          'xtick.major.pad': 7,
                          'ytick.major.pad': 7,
                          'ytick.color': nearly_black,
                          'xtick.color': nearly_black,
                          'axes.labelcolor': nearly_black,
			              'axes.spines.bottom': True,
                          'axes.spines.right': False,
                          'axes.spines.top': False,
                          'axes.spines.left': True,
                          'pdf.fonttype': 42,
                          'ps.fonttype': 42,
                          'mathtext.fontset': 'custom',
                          'font.size': fontsize,
                          'text.usetex' : True,
                          'text.usetex': True,
                          'savefig.bbox':'tight',
                          'axes.facecolor': white,
                          'axes.labelpad': 10.0,
		                  'axes.grid.axis': 'y',
			              'axes.grid': True,
                          'grid.color': light_grey,
                          'axes.labelsize': fontsize,
                          'axes.titlepad': 30,
                          'lines.markersize': 7.0,
                          'lines.scale_dashes': False,
                          'xtick.labelsize': fontsize * 0.8,
                          'ytick.labelsize': fontsize * 0.8}

for k, v in master_formatting.items():
    rcParams[k] = v

color_cycle = tableau.values()
try:
    from matplotlib import cycler
    rcParams['axes.prop_cycle'] = cycler(color=color_cycle)
except Exception:
    raise