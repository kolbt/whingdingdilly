'''
A simple way of looking at LOCAL distribution of TYPE:
    1. Create a mesh of the plane
    2. Populate mesh with particles
    3. Count
        a. number of a type
        b. number of b type
    4. Report as...
        raw number of A/B? (doesn't account for dense interior... )
        percentage of A/B normalized by total number of particles in bin? (must ignore gas)
        percentage of A/B normalized by total number of A/B? (doesn't ignore gas) *** THIS ONE ***
    ? Additional requirements ?
        - Bin must have greater than min_number particles... else = 0
          /// normalize = ((bin_part > min_number) ? 1 : 0)     /// ternary c++
          /// normalize = (1 if (bin_part > min_number) else 0) /// ternary python
          /// if normalize: ...                                 ///
'''

# Imports and loading the .gsd file
import sys

pe_a = int(sys.argv[1])
pe_b = int(sys.argv[2])
part_perc_a = int(sys.argv[3])
hoomd_path = str(sys.argv[4])
gsd_path = str(sys.argv[5])

part_frac_a = float(part_perc_a) / 100.0

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

boltz = 1                           # lj boltzmann constant
temp = 1                            # temperature
trans_diff = 1                      # translational diffusion coefficient
drag = boltz * temp / trans_diff    # computed drag coefficient

# READ IN .GSD AND RAW DATA
myfile = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"
f = hoomd.open(name=myfile, mode='rb')
start = 0
dumps = f.__len__()
#dumps = 72
position_array = np.zeros((dumps), dtype=np.ndarray)        # array of position arrays
type_array = np.zeros((dumps), dtype=np.ndarray)            # particle types
box_data = np.zeros((dumps), dtype=np.ndarray)              # box dimensions
my_dt = 0.000001

with hoomd.open(name=myfile, mode='rb') as t:               # open for reading
    snap = t[0]                                             # snap oth snapshot
    box_data = snap.configuration.box                       # get box dimensions
    for i in range(start, dumps):
        snap = t[i]
        type_array[i] = snap.particles.typeid
        position_array[i] = snap.particles.position         # store all particle positions

part_num = len(position_array[start])
points = np.zeros((part_num), dtype=np.ndarray)
diameter = 1

float_side = box_data[0] / 2.0
side = float((int(box_data[0])+1)) / 2
box_width = 5                                           # width of bins
while side % box_width != 0:                            # grid boxes of width 1 sigma
    side += 1                                           # make sure side length is divisible

spacer = int(side * 2 / (box_width * diameter))                 # spacing between boxes
occ = 100                                                       # max occupancy of an index
bin_part = np.zeros((dumps, spacer, spacer), dtype = np.int)    # total particles in bin
bin_A = np.zeros((dumps, spacer, spacer), dtype = np.int)       # total A type in bin
bin_B = np.zeros((dumps, spacer, spacer), dtype = np.int)       # total B type in bin
norm_A = np.zeros((dumps, spacer, spacer), dtype = np.float64)  # normalized number of A
norm_B = np.zeros((dumps, spacer, spacer), dtype = np.float64)  # normalized number of B

for iii in range(start, dumps):
    
    points = position_array[iii]
    points[:][:][2] = 0

    for jjj in range(part_num):
        
        loc_x = int((points[jjj][0] + float_side) / (box_width * diameter)) # x index in mesh
        loc_y = int((points[jjj][1] + float_side) / (box_width * diameter)) # y index in mesh
        
        bin_part[iii][loc_x][loc_y] += 1    # total number of particles in bin
        
        if type_array[iii][jjj] == 0:
            bin_A[iii][loc_x][loc_y] += 1   # total number of A particles in bin

    bin_B[iii] = bin_part[iii] - bin_A[iii]

# Get total number of A particles
tot_A = 0
tot_B = 0
for iii in range(0,part_num):
    if type_array[start][iii] == 0:
        tot_A += 1
    else:
        tot_B += 1

for iii in range(start, dumps):
    for jjj in range(spacer):
        for kkk in range(spacer):
            norm_A[iii][jjj][kkk] = float(bin_A[iii][jjj][kkk]) / float(tot_A)
            norm_B[iii][jjj][kkk] = float(bin_B[iii][jjj][kkk]) / float(tot_B)

# Let's plot it and take a look
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

min = -(float(box_data[0]/2))
max = (float(box_data[0]/2))

current_mesh = np.zeros((spacer, spacer), dtype=np.float32)
for iii in range(start,dumps):
    current_mesh = norm_A[iii].T
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.imshow(current_mesh, origin='lower')
    ax.set_aspect('equal')
    plt.colorbar(orientation='vertical')
    #    plt.clim(0,8)
    plt.savefig('test_dense_A_'+
                str(iii)+
                '.png', dpi=1000)
    plt.close()

current_mesh = np.zeros((spacer, spacer), dtype=np.float32)
for iii in range(start,dumps):
    current_mesh = norm_B[iii].T
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.imshow(current_mesh, origin='lower')
    ax.set_aspect('equal')
    plt.colorbar(orientation='vertical')
    #    plt.clim(0,8)
    plt.savefig('test_dense_B_'+
                str(iii)+
                '.png', dpi=1000)
    plt.close()


