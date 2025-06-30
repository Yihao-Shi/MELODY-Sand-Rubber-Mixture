import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection

def fc(timestep, content):
    Pnum = 5025
    Rnum = math.ceil((1-float(content)/100)*Pnum)
    Snum = Pnum-Rnum
    ax = plt.subplot(111,aspect='equal')
    loc = np.loadtxt('../../Simulation_DATA/BiaxialShearTest/'+content+'percent/STRlocation'+timestep+'.txt')
    lwidth = np.loadtxt('../../Simulation_DATA/BiaxialShearTest/'+content+'percent/STRlinewidth'+timestep+'.txt')
    
    Spatches = []
    Rpatches = []
    count = 0
    for p in loc:
        poly_coords = []
        nodes=int(p[-2])
        for i in range(0,nodes-1):
            poly_coords.append((p[i],p[i+nodes]))
        polygon = Polygon(xy=poly_coords)
        if p[-1]<Snum:
            Spatches.append(polygon)
            count+=1
        else:
            Rpatches.append(polygon)
            count+=1
            
    p = PatchCollection(Spatches,edgecolors='#6699cc',facecolors='#B2CCE5',linewidth=0.25)
    ax.add_collection(p)
    p = PatchCollection(Rpatches,edgecolors='#A8A8A8',facecolors='#CDCDCD',linewidth=0.25)
    ax.add_collection(p)
    
    '''for f in lwidth:
        if f[5]==0:
            ssline = Line2D([f[0],f[2]],[f[1],f[3]],lw=f[4]/56000.,color='black')
            ax.add_line(ssline)
        elif f[5]==1:
            srline = Line2D([f[0],f[2]],[f[1],f[3]],lw=f[4]/56000.,color='blue')
            ax.add_line(srline)
        elif f[5]==2:
            rrline = Line2D([f[0],f[2]],[f[1],f[3]],lw=f[4]/56000.,color='red')
            ax.add_line(rrline)'''
   
    ax.set_xlim([-0.7,0.7])
    ax.set_ylim([-1.3,1.3])
    ax.axison=False
    plt.savefig('GRAPH/'+content+'STR%rve'+timestep+'.tif',format='tif',dpi=600)

for i in range(1,37):
    fc(str(i),'0')
