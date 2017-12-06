'''
What we really want to obtain is the divergence...

Perhaps, with small enough bins, a net force of 0
would correspond to point sources and sinks...

Let's try it!
    1. Mesh the space, assign points to boxes
    2. Compute net active force tensor in each box
    3. Plot as matrix (using below discussion as reference)
https://stackoverflow.com/questions/16492830/colorplot-of-2d-array-matplotlib
'''

# Imports and loading the .gsd file
import sys

hoomd_path = str(sys.argv[4])
gsd_path = str(sys.argv[5])

# need to extract values from filename (pa, pb, xa) for naming
part_perc_a = int(sys.argv[3])
part_frac_a = float(part_perc_a) / 100.0
pe_a = int(sys.argv[1])
pe_b = int(sys.argv[2])

sys.path.append(hoomd_path)

import hoomd
from hoomd import md
from hoomd import deprecated
from mypy import load_bar
sys.path.append(gsd_path)
import gsd
from gsd import hoomd
from gsd import pygsd
import numpy as np
import math

def quat_to_theta(quat):
    x = quat[1]
    y = quat[2]
    rad = math.atan2(y,x) # gives values from [-pi,pi]
    return rad

boltz = 1
temp = 1
trans_diff = 1
drag = boltz * temp / trans_diff

# READ IN .GSD AND RAW DATA
myfile = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"
f = hoomd.open(name=myfile, mode='rb')
start = 990
dumps = f.__len__()
position_array = np.zeros((dumps), dtype=np.ndarray)        # array of position arrays
type_array = np.zeros((dumps), dtype=np.ndarray)            # particle types
orientations = np.zeros((dumps), dtype=np.ndarray)          # orientations
box_data = np.zeros((dumps), dtype=np.ndarray)              # box dimensions

with hoomd.open(name=myfile, mode='rb') as t:               # open for reading
    snap = t[0]                                             # snap oth snapshot
    box_data = snap.configuration.box                       # get box dimensions
    for i in range(start, dumps):
        snap = t[i]
        type_array[i] = snap.particles.typeid
        position_array[i] = snap.particles.position         # store all particle positions
        orientations[i] = snap.particles.orientation        # store the orientation of all particles

part_num = len(position_array[start])
points = np.zeros((part_num), dtype=np.ndarray)
active_pressure = np.zeros((dumps, part_num), dtype=np.float32)
rads = np.zeros((dumps, part_num), dtype=np.float32)    # orientations as angles
diameter = 1

for iii in range(start,dumps):
    for jjj in range(0,part_num):
        rads[iii][jjj] = quat_to_theta(orientations[iii][jjj])

# CREATE MESH
float_side = box_data[0] / 2.0
side = float((int(box_data[0])+1)) / 2
box_width = 10
while side % box_width != 0:                            # grid boxes of width 1 sigma
    side += 1                                           # make sure side length is divisible

spacer = int(side * 2 / (box_width * diameter))                     # spacing between indiviudal boxes
occ = 100                                                           # max occupancy of a grid box
show_bins = np.zeros((dumps, spacer, spacer), dtype = np.float32)   # array of each grid

# CREATE MESH, POPULATE WITH NET ACTIVE FORCE FOR ALL TSTEPS
load_bar.printLoadBar(0, dumps, prefix = "Progress:", suffix = "Complete")
for iii in range(start, dumps):
    points = position_array[iii]
    
    # Recreate the array (re-zero)
    mesh = np.zeros((spacer, spacer), dtype = np.ndarray)   # array of each grid
    test_occ = np.zeros((occ, 5))                           # occupation test index
    for j in range(0, spacer):
        for k in range(0, spacer):
            mesh[j][k] =  np.zeros_like(test_occ)

    # PLACE PARTICLES IN MESH
    for jjj in range(part_num):
        # get the index of the mesh the particle belongs in
        loc_x = int((points[jjj][0] + float_side) / (box_width * diameter))
        loc_y = int((points[jjj][1] + float_side) / (box_width * diameter))
        # place the particle in the first unoccupied space in this quadrant list
        for s in range(1, occ):                             # index 0 holds net active force
            if mesh[loc_x][loc_y][s][3] == 0:               # test occupancy of list
                mesh[loc_x][loc_y][s][0] = points[jjj][0]   # x coord
                mesh[loc_x][loc_y][s][1] = points[jjj][1]   # y coord
                mesh[loc_x][loc_y][s][2] = jjj              # particle id
                mesh[loc_x][loc_y][s][3] = 1                # switch flag to occupied
                break

    # COMPUTE ACTIVE FORCE TENSOR FOR EACH PARTICLE
    for jjj in range(part_num):
        if type_array[iii][jjj] == 0:
            sigma_x = drag * pe_a * math.cos(rads[iii][jjj])# * points[jjj][0]
            sigma_y = drag * pe_a * math.sin(rads[iii][jjj])# * points[jjj][1]
            active_pressure[iii][jjj] = -0.5 * (sigma_x + sigma_y)
        else:
            sigma_x = drag * pe_b * math.cos(rads[iii][jjj])# * points[jjj][0]
            sigma_y = drag * pe_b * math.sin(rads[iii][jjj])# * points[jjj][1]
            active_pressure[iii][jjj] = -0.5 * (sigma_x + sigma_y)

    # SUM ACTIVE FORCE TENSOR IN EACH BIN
    for jjj in range(part_num):
        loc_x = int((points[jjj][0] + float_side) / (box_width * diameter))
        loc_y = int((points[jjj][1] + float_side) / (box_width * diameter))
        show_bins[iii][loc_x][loc_y] += active_pressure[iii][jjj]
#        mesh[loc_x][loc_y][0][0] += active_pressure[iii][jjj]
    load_bar.printLoadBar(iii+1, dumps, prefix = "Progress:", suffix = "Complete")


# PLOT BINNED DATA

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

min = -(float(box_data[0]/2))
max = (float(box_data[0]/2))

# Make a better colormap
import matplotlib.colors as col

current_mesh = np.zeros((spacer, spacer), dtype=np.float32)
for iii in range(start,dumps):
    current_mesh = show_bins[iii]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.imshow(current_mesh)
    ax.set_aspect('equal')
    plt.colorbar(orientation='vertical')
    plt.savefig('test_'+
                str(iii)+
                '.png', dpi=1000)
    plt.close()


#for mmm in range(0,dumps):
#    xs = np.zeros((part_num), dtype = np.float32)
#    ys = np.zeros((part_num), dtype = np.float32)
#    for iii in range(0,part_num):
#        xs[iii] = position_array[mmm][iii][0]
#        ys[iii] = position_array[mmm][iii][1]
#
#    plt.scatter(xs, ys, s=0.20, c=active_pressure, edgecolors='none')
#    plt.colorbar()
#    plt.xlim(min, max)
#    plt.ylim(min, max)
#    plt.savefig('active_pressure_pa'+
#                str(pe_a) +
#                '_pb'+
#                str(pe_b)+
#                '_xa'+
#                str(part_perc_a)
#                +'_mvout_'+
#                str(mmm)+
#                '.png',
#                dpi=1000)
#    plt.close()







