import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection

content = '40'
timestep1 = '1'
timestep2 = '241'
Pnum = 797
Rnum = math.ceil((1-float(content)/100)*Pnum)
Snum = Pnum-Rnum
ax = plt.subplot(111,aspect='equal')
loc_ini = np.loadtxt('../../Simulation_DATA/DirectShearTest/'+content+'percent/location'+timestep1+'.txt')
lwidth_ini = np.loadtxt('../../Simulation_DATA/DirectShearTest/'+content+'percent/linewidth'+timestep1+'.txt')
loc_end = np.loadtxt('../../Simulation_DATA/DirectShearTest/'+content+'percent/location'+timestep2+'.txt')
lwidth_end = np.loadtxt('../../Simulation_DATA/DirectShearTest/'+content+'percent/linewidth'+timestep2+'.txt')

def fc(loc,lwidth,ch,content):
    Spatches = []
    Rpatches = []
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
            
    p = PatchCollection(Spatches,edgecolors='#6699cc',facecolors='#B2CCE5',linewidth=0.25)
    ax.add_collection(p)
    p = PatchCollection(Rpatches,edgecolors='#A8A8A8',facecolors='#CDCDCD',linewidth=0.25)
    ax.add_collection(p)
    
    for f in lwidth:
        if f[5]==0:
            ssline = Line2D([f[0],f[2]],[f[1],f[3]],lw=f[4]/8000.,color='black')
            ax.add_line(ssline)
        elif f[5]==1:
            srline = Line2D([f[0],f[2]],[f[1],f[3]],lw=f[4]/8000.,color='blue')
            ax.add_line(srline)
        elif f[5]==2:
            rrline = Line2D([f[0],f[2]],[f[1],f[3]],lw=f[4]/8000.,color='red')
            ax.add_line(rrline)
   
    ax.set_xlim([-0.4,0.4])
    ax.set_ylim([-0.4,0.4])
    ax.axison=False
    plt.savefig('GRAPH/'+content+'%rve'+ch+'.tif',format='tif',dpi=300)

#fc(loc_ini,lwidth_ini,timestep1,content)
fc(loc_end,lwidth_end,timestep2,content)
