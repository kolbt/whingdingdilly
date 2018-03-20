import sys
import matplotlib
#matplotlib.use('Agg')
#matplotlib.use('macosx')
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from colorsys import hls_to_rgb

#https://stackoverflow.com/questions/37299142/how-to-set-a-colormap-which-can-give-me-over-20-distinct-colors-in-matplotlib
# User Neocortex
def get_distinct_colors(n):
    colors = []
    for i in np.arange(0., 360., 360. / n):
        h = i / 360.
        l = (50 + np.random.rand() * 10) / 100.
        s = (90 + np.random.rand() * 10) / 100.
        colors.append(hls_to_rgb(h, l, s))
    return colors

my_cols = get_distinct_colors(17)
print(my_cols)

'''
    The data you are drawing from is stacked in set Pe_B planes. So,
    rows 1-11 (that's 11 rows) are all binary values which indicate
    phase separation ( 0 = no, 1 = yes) at Pe_B = 0.  Thus, rows 
    12-22 are for Pe_B = 10 etc. etc. until Pe_B = 150 (16 total
    planes).
    
    Ultimately, Mathematica was better at this, check out 
    "exp_theory_overlay.nb"
'''

#file = str(sys.argv[1])
file = "/Users/kolbt/Desktop/phase_3D.txt"
data = np.loadtxt(file, dtype=np.int8)
# start with 1D array
peb = np.zeros((16), dtype=np.ndarray)
var = 0
for i in range(0,16):
    peb[i] = data[var:var+11, :16]
    var += 11

# for loop the fuck outta this
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

# add together set points over pb (min = 0, max = 16)
plane_heatmap = np.zeros((16, 11), dtype = np.int)
for iii in range(0,16):
    for jjj in range(0,11):
        for lll in range(0,16):
            plane_heatmap[iii][jjj] += phase_3D[iii][lll][jjj]

print(plane_heatmap)

plt.imshow(plane_heatmap.T, origin='lower', cmap='gnuplot_r', vmin=0, vmax=16)
plt.colorbar()
plt.show()

## start plotting the data as voxels
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#
#for idx in range(0,16):
#    for idy in range(0,16):
#        for idz in range(0,11):
#            if phase_3D[idx,idy,idz] != 0:
#                ax.scatter(idx, idy, idz, c=my_cols[idz], s=40, alpha=0.8, depthshade=True)

