import numpy as np
import math
from matplotlib import cm, colors
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

content = '10'
timestep = '9'
Pnum = 5025
Rnum = math.ceil((1-float(content)/100)*Pnum)
Snum = Pnum-Rnum
ax = plt.subplot(111,aspect='equal')
loc = np.loadtxt('../../Simulation_DATA/BiaxialShearTest/'+content+'percent/location'+timestep+'.txt')
angle = np.loadtxt('../../Simulation_DATA/BiaxialShearTest/'+content+'percent/rotation'+timestep+'.txt')
angle0 = np.loadtxt('../../Simulation_DATA/BiaxialShearTest/'+content+'percent/rotation1.txt')

Spatches, Rpatches, X, Y = [], [], [], []
count = 0
for p in loc:
    poly_coords = []
    nodes=int(p[-1])
    for i in range(0,nodes-1):
        poly_coords.append((p[i],p[i+nodes]))
    polygon = Polygon(xy=poly_coords)
    if count<Snum:
        Spatches.append(polygon)
        count+=1
    else:
        Rpatches.append(polygon)
        count+=1 
        X.append(np.mean(p[0:64]))
        Y.append(np.mean(p[64:128]))

normlizer = plt.Normalize(vmin = -1, vmax = 1)
cmap=cm.get_cmap('RdBu')
sm=cm.ScalarMappable(norm=normlizer,cmap=cmap)
particledata = (angle-angle0)
color = sm.to_rgba(particledata)

p = PatchCollection(Spatches,edgecolors='#6699cc',facecolors='#B2CCE5',linewidth=0.25)
ax.add_collection(p)
p = PatchCollection(Rpatches,edgecolors='#6699cc',facecolors=color,linewidth=0.25)
ax.add_collection(p)

ax.set_xlim([-0.8,0.8])
ax.set_ylim([-1.5,1.5])
ax.axison=False
plt.savefig('GRAPH/'+content+'%rveRotation'+timestep+'.tif',format='tif',dpi=1200)
