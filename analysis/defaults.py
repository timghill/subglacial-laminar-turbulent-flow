"""

Set default parameters and plot options

"""

import numpy as np
import cmocean
from palettes.code import palettes


# colors = np.array([[0.579, 0.677, 0.781, 1],
#                    [0.199, 0.328, 0.492, 1],
#                    [0.500, 0.500, 0.500, 1],
#                    [0.929, 0.835, 0.408, 1],
#                    [0.859, 0.683, 0.275, 1]])

colors = np.array([[0.579, 0.677, 0.781, 1],
                   [0.199, 0.328, 0.492, 1],
                   [0.250, 0.250, 0.250, 1],
                   [0.929, 0.835, 0.408, 1],
                   [0.836, 0.590, 0.160, 1]])

linewidths = np.array([2.5, 1.25, 0.5, 2.5, 1.25])
linestyles = ['solid', 'solid', 'dashed', 'solid', 'solid']
zorders = [2, 2, 3, 2, 2]

tslice = 365 + 190

x_bands = [15, 30, 70]
band_width = 5

labels = ['Turbulent 5/4', 'Turbulent 3/2', 'Laminar', 'Transition 5/4', 'Transition 3/2']

cmaps = {   'floatation':cmocean.cm.dense,
            'Q':palettes.get_cmap('BrownYellow'),
        }
