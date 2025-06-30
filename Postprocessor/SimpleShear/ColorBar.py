# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 18:28:41 2022

@author: ELEVEN
"""

import matplotlib.pyplot as plt

import matplotlib as mpl

figure, axes = plt.subplots(figsize =(5, 0.4))

figure.subplots_adjust(bottom = 0.5)

color_map = mpl.cm.jet

normlizer = mpl.colors.Normalize(vmin = 0, vmax = 14)

figure.colorbar(mpl.cm.ScalarMappable(norm = normlizer,

cmap = color_map),

cax = axes, orientation ='horizontal',

label ='Average sand particle rotation (rad)')

plt.savefig('GRAPH/ColorBar.tif',format='tif',dpi=300)
