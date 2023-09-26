
_colors = [
        '#d60000', '#e2afaf', '#018700', '#a17569', '#e6a500', '#004b00',
        '#6b004f', '#573b00', '#005659', '#5e7b87', '#0000dd', '#00acc6',
        '#bcb6ff', '#bf03b8', '#645472', '#790000', '#0774d8', '#729a7c',
        '#8287ff', '#ff7ed1', '#8e7b01', '#9e4b00', '#8eba00', '#a57bb8',
        '#5901a3', '#8c3bff', '#a03a52', '#a1c8c8', '#f2007b', '#ff7752',
        '#bac389', '#15e18c', '#60383b', '#546744', '#380000', '#e252ff',
        ]

def get_cmap(cmap=None, N=None, category=False, theme='light'):

    import re
    import matplotlib as mpl

    if theme == 'light':
        facecolor = 'w'
    elif theme == 'dark':
        facecolor = 'k'

    is_hex_color = lambda x: re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', x)
    if isinstance(cmap, str) and is_hex_color(cmap):
        cmap = [facecolor, cmap]

    if isinstance(cmap, list):
        if category:
            if N is None:
                N = len(cmap)
            cmap = mpl.colors.ListedColormap(cmap[:N], N=N)
        else:
            cmap = mpl.colors.LinearSegmentedColormap.from_list('cmap', 
                    cmap, N=256)
    elif isinstance(cmap, str):
        if not category:
            N = 256
        else:
            assert N is not None
        cmap = mpl.cm.get_cmap(cmap, N)
        if category and isinstance(cmap, mpl.colors.LinearSegmentedColormap):
            cmap = [mpl.colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
            cmap = mpl.colors.ListedColormap(cmap, N=N)
    return cmap

from .gem import Gem
from .stoarr import Stoarr
from .segplot import seg_spatial
from .coords3d import Coords
from .mipplot import mip_spatial

