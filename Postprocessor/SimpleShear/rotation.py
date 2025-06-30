# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 13:23:23 2022

@author: ELEVEN
"""

import numpy as np
import math
from matplotlib import cm, colors
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

content = '40'
timestep = '241'
Pnum = 797
Rnum = math.ceil((1-float(content)/100)*Pnum)
Snum = Pnum-Rnum
ax = plt.subplot(111,aspect='equal')
loc = np.loadtxt('../../Simulation_DATA/DirectShearTest/'+content+'percent/location'+timestep+'.txt')
angle0 = np.loadtxt('../../Simulation_DATA/DirectShearTest/0percent/rotation'+timestep+'.txt')
angle = np.loadtxt('../../Simulation_DATA/DirectShearTest/'+content+'percent/rotation'+timestep+'.txt')

Spatches, Rpatches, X, Y = [], [], [], []
count = 0
for p in loc:
    poly_coords = []
    for i in range(0,63):
        poly_coords.append((p[i],p[i+64]))
    polygon = Polygon(xy=poly_coords)
    if count<Snum:
        Spatches.append(polygon)
        count+=1
    else:
        Rpatches.append(polygon)
        count+=1 
        X.append(np.mean(p[0:64]))
        Y.append(np.mean(p[64:128]))

particledata = (angle-np.min(angle0))/(np.max(angle0)/2-np.min(angle0))
color = cm.jet(particledata)

p = PatchCollection(Spatches,edgecolors='#6699cc',facecolors='#B2CCE5',linewidth=0.25)
ax.add_collection(p)
p = PatchCollection(Rpatches,edgecolors='#6699cc',facecolors=color,linewidth=0.25)
ax.add_collection(p)

ax.set_xlim([-0.4,0.4])
ax.set_ylim([-0.4,0.4])
ax.axison=False
plt.savefig('GRAPH/'+content+'%rveRotation'+timestep+'.tif',format='tif',dpi=300)