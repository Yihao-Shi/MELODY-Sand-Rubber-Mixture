# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 18:28:41 2022

@author: ELEVEN
"""

import matplotlib.pyplot as plt

import matplotlib as mpl

figure, axes = plt.subplots(figsize =(10, 1.2))

figure.subplots_adjust(bottom = 0.5)

color_map = mpl.cm.RdBu

normlizer = mpl.colors.Normalize(vmin = -1, vmax = 1)

a = figure.colorbar(mpl.cm.ScalarMappable(norm = normlizer,

cmap = color_map),

cax = axes, orientation ='horizontal')

a.ax.tick_params(labelsize=16)

a.set_label("Accmulated sand particle rotation (rad)",fontdict={"size": 16})

plt.savefig('GRAPH/ColorBar.tif',format='tif',dpi=600)
