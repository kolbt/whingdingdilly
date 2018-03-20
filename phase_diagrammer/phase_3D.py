import sys
import matplotlib
#matplotlib.use('Agg')
#matplotlib.use('macosx')
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
#from mayavi import mlab

#mlab.test_plot3d()
#mlab.show()

my_cols = ['b', 'g', 'r', 'c', 'm', 'y', 'b', 'g', 'r', 'c', 'm']

# Some code from stack to plot voxels
def cuboid_data(pos, size=(1,1,1)):
    # code taken from
    # https://stackoverflow.com/a/35978146/4124317
    # suppose axis direction: x: to left; y: to inside; z: to upper
    # get the (left, outside, bottom) point
    o = [a - b / 2 for a, b in zip(pos, size)]
    # get the length, width, and height
    l, w, h = size
    x = [[o[0], o[0] + l, o[0] + l, o[0], o[0]],
         [o[0], o[0] + l, o[0] + l, o[0], o[0]],
         [o[0], o[0] + l, o[0] + l, o[0], o[0]],
         [o[0], o[0] + l, o[0] + l, o[0], o[0]]]
    y = [[o[1], o[1], o[1] + w, o[1] + w, o[1]],
         [o[1], o[1], o[1] + w, o[1] + w, o[1]],
         [o[1], o[1], o[1], o[1], o[1]],
         [o[1] + w, o[1] + w, o[1] + w, o[1] + w, o[1] + w]]
    z = [[o[2], o[2], o[2], o[2], o[2]],
         [o[2] + h, o[2] + h, o[2] + h, o[2] + h, o[2] + h],
         [o[2], o[2], o[2] + h, o[2] + h, o[2]],
         [o[2], o[2], o[2] + h, o[2] + h, o[2]]]
    return x, y, z

def plotCubeAt(pos=(0,0,0),ax=None, z_ind=0):
    # Plotting a cube element at position pos
    if ax !=None:
        X, Y, Z = cuboid_data( pos )
        ax.plot_surface(X, Y, Z, color=my_cols[z_ind], rstride=1, cstride=1, alpha=0.8)

# start with 1D array

'''
    The data you are drawing from is stacked in set Pe_B planes. So,
    rows 1-11 (that's 11 rows) are all binary values which indicate
    phase separation ( 0 = no, 1 = yes) at Pe_B = 0.  Thus, rows 
    12-22 are for Pe_B = 10 etc. etc. until Pe_B = 150 (16 total
    planes).
    
    Ultimately, Mathematica was better at this, check out 
    "exp_theory_overlay.nb"
'''

file = str(sys.argv[1])
data = np.loadtxt(file, dtype=np.int8)

peb = np.zeros((16), dtype=np.ndarray)
var = 0
for i in range(0,16):
    peb[i] = data[var:var+11, :16]
    var += 11
    print(peb[i])

# for loop the fuck outta this
#phase_3D = np.zeros((11,16,16), dtype=np.int8)
#for j in range(0,16):
#    r=0
#    c=0
#    for i in range(0,len(data)):
#        phase_3D[r][c][j] = data[i][j]
#        r += 1
#        if r == 11:
#            r = 0
#            c += 1

phase_3D = np.zeros((16,16,11), dtype=np.int8)
for j in range(0,16):
    r=0
    c=0
    for i in range(0,len(data)):
        phase_3D[c][j][r] = data[i][j]
        r += 1
        if r == 11:
            r = 0
            c += 1

# start plotting the data as voxels
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.set_aspect('equal')

#for idx in range(0,11):
#    for idy in range(0,16):
#        for idz in range(0,16):
#            if phase_3D[idx,idy,idz] != 0:
##                ax.scatter(idx, idy, idz, c='b', marker='s', s=20)
#                plotCubeAt(pos=(idx-0.5, idy-0.5, idz-0.5), ax=ax)
##            else:
##                ax.scatter(idx, idy, idz, c='w')

#my_cols = np.arange(11.0)
#my_cols /= 10.0
#colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
#my_cm = LinearSegmentedColormap.from_list('my_map', my_cols, N=10)
for idx in range(0,16):
    for idy in range(0,16):
        for idz in range(0,11):
            if phase_3D[idx,idy,idz] != 0:
                ax.scatter(idx, idy, idz, c=my_cols[idz], s=40, alpha=0.8, depthshade=True)

#                if idz % 2 == 0:
#                    ax.scatter(idx, idy, idz, marker='.', s=40, c='b', depthshade=True)
#                else:
#                    ax.scatter(idx, idy, idz, marker='.', s=40, c='k', depthshade=True)
#                plotCubeAt(pos=(idx-0.5, idy-0.5, idz-0.5), ax=ax, z_ind=idz)
#            else:
#                ax.scatter(idx, idy, idz, c='w')

#ax.set_xlabel('Particle Fraction of A')
#ax.set_xlim3d(0,10)
#ax.set_ylabel('Pe_A')
#ax.set_ylim3d(0,15)
#ax.set_zlabel('Pe_B')
#ax.set_zlim3d(0,15)
#plt.show()

#ax.controuf()
ax.set_xlabel('Pe_A')
ax.set_xlim3d(0,15)
ax.set_ylabel('Pe_B')
ax.set_ylim3d(0,15)
ax.set_zlabel('Particle Fraction of A')
ax.set_zlim3d(0,10)
#plt.colorbar()
plt.show()

# WORD. Now plot that shit
#plt.plot(phase_3D, cmap='Blues')
fig.savefig('phase_3D.png', dpi=1000)
plt.close()

#mlab.contour3d(phase_3D, contours=5)

