'''
#                           This is an 80 character line                       #
This file is intended to look at the average distance between particles in the
dense phase (delta). This should directly result from the forces applied on
particles in the dense phase (that's how potentials work), but, we want to make
sure we are using an appropriate method of implementing the well depth
(epsilon).

What this file does:
    1. Read in your .gsd file
    2. Mesh the system
    3. Compute distance between 'close' particles (use r_cut=1.122)
    4. Store as lists of three types: AA, AB, BB
    5. Plot lists as population histograms
'''

# Imports and loading the .gsd file
import sys

pe_a = int(sys.argv[1])                     # activity A
pe_b = int(sys.argv[2])                     # activity B
part_perc_a = int(sys.argv[3])              # percentage A particles
part_frac_a = float(part_perc_a) / 100.0    # fraction A particles
hoomd_path = str(sys.argv[4])               # local path to hoomd-blue
gsd_path = str(sys.argv[5])                 # local path to gsd

sys.path.append(hoomd_path)     # ensure hoomd is in your python path
sys.path.append(gsd_path)       # ensure gsd is in your python path

import hoomd
from hoomd import md
from hoomd import deprecated

import gsd
from gsd import hoomd
from gsd import pygsd

import freud
from freud import parallel
from freud import box
from freud import density
from freud import cluster

import numpy as np
from scipy.interpolate import griddata

import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
sns.set(color_codes=True)

import math

# Any functions I'll need go here
def getDistance(point1, point2x, point2y):
    "Find the distance between two points"
    distance = np.sqrt((point2x - point1[0])**2 + (point2y - point1[1])**2)
    return distance

# File to read from
in_file = "pa"+str(pe_a)+\
"_pb"+str(pe_b)+\
"_xa"+str(part_perc_a)+\
".gsd"

# File base
b_file = "pa"+str(pe_a)+\
"_pb"+str(pe_b)+\
"_xa"+str(part_perc_a)

f = hoomd.open(name=in_file, mode='rb') # open gsd file with hoomd
dumps = f.__len__()                     # get number of timesteps dumped

start = 0       # gives first frame to read
end = dumps     # gives last frame to read
end = 2

start = 400
end = 402

positions = np.zeros((end), dtype=np.ndarray)       # array of positions
types = np.zeros((end), dtype=np.ndarray)           # particle types
box_data = np.zeros((1), dtype=np.ndarray)          # box dimensions
timesteps = np.zeros((end), dtype=np.float64)       # timesteps

# Get relevant data from .gsd file
with hoomd.open(name=in_file, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    for iii in range(start, end):
        snap = t[iii]                               # snapshot of frame
        types[iii] = snap.particles.typeid          # get types
        positions[iii] = snap.particles.position    # get positions
        timesteps[iii] = snap.configuration.step    # get timestep

timesteps -= timesteps[0]       # get rid of brownian run time

# Get number of each type of particle
part_num = len(types[start])
part_A = int(part_num * part_frac_a)
part_B = part_num - part_A

# Feed data into freud analysis software
l_box = box_data[0]
h_box = l_box / 2.0
a_box = l_box * l_box
f_box = box.Box(Lx = l_box, Ly = l_box, is2D = True)    # make freud box

# Make the mesh
nBins = 200
sizeBin = l_box / nBins
r_cut = 1.122
#sizeBin = r_cut
#nBins = int(l_box / sizeBin)

print(sizeBin)

for iii in range(start, end):
    
    # Mesh array
    binParts = [ [ [] for b in range(nBins) ] for a in range(nBins) ]

    # Easier accessors
    pos = positions[iii]
    pos = np.delete(pos, 2, 1)
    typ = types[iii]
    
    # Instantiate pair-wise lists
    ALL = []
    AA = []
    BB = []
    AB = []

    # Put particles in their respective bins
    for jjj in range(0, part_num):
        # Get mesh indices
        tmp_posX = pos[jjj][0] + h_box
        tmp_posY = pos[jjj][1] + h_box
        x_ind = int(tmp_posX / sizeBin)
        y_ind = int(tmp_posY / sizeBin)
        # Append particle id to appropriate bin
        binParts[x_ind][y_ind].append(jjj)

    countALL = 0
    countAA = 0
    countBB = 0
    countAB = 0

    # Compute distance, each pair will be counted twice
    for jjj in range(0, part_num):
        # Get mesh indices
        tmp_posX = pos[jjj][0] + h_box
        tmp_posY = pos[jjj][1] + h_box
        x_ind = int(tmp_posX / sizeBin)
        y_ind = int(tmp_posY / sizeBin)
        # Get index of surrounding bins
        l_bin = x_ind - 1   # index of left bins
        r_bin = x_ind + 1   # index of right bins
        b_bin = y_ind - 1   # index of bottom bins
        t_bin = y_ind + 1   # index of top bins
        if r_bin == nBins:
            r_bin -= nBins  # adjust if wrapped
        if t_bin == nBins:
            t_bin -= nBins  # adjust if wrapped
        h_list = [l_bin, x_ind, r_bin]  # list of horizontal bin indices
        v_list = [b_bin, y_ind, t_bin]  # list of vertical bin indices

        # Loop through all bins
        for h in range(0, len(h_list)):
            for v in range(0, len(v_list)):
                # Take care of periodic wrapping for position
                wrapX = 0.0
                wrapY = 0.0
                if h == 0 and h_list[h] == -1:
                    wrapX -= l_box
                if h == 2 and h_list[h] == 0:
                    wrapX += l_box
                if v == 0 and v_list[v] == -1:
                    wrapY -= l_box
                if v == 2 and v_list[v] == 0:
                    wrapY += l_box
                # Compute distance between particles
                for b in range(0, len(binParts[h_list[h]][v_list[v]])):
                    ref = binParts[h_list[h]][v_list[v]][b]
                    r = getDistance(pos[jjj],
                                    pos[ref][0] + wrapX,
                                    pos[ref][1] + wrapY)
                    
                    # If LJ potential is on, store into a list (omit self)
                    if 0.1 < r <= r_cut:
                        ALL.append(r)                           # All particles
                        countALL += 1
                        if typ[jjj] == 0 and typ[ref] == 0:     # AA distance
                            AA.append(r)
                            countAA += 1
                        elif typ[jjj] == 1 and typ[ref] == 1:   # BB distance
                            BB.append(r)
                            countBB += 1
                        else:                                   # AB distance
                            AB.append(r)
                            countAB += 1

    popALL = len(ALL) / 2
    popAA = len(AA) / 2
    popAB = len(AB) / 2
    popBB = len(BB) / 2

    # Plot histogram of data
    sns.distplot(ALL, norm_hist=True, bins=100)
    plt.savefig('ALL_' + b_file + "_" + str(iii) + '.png', dpi=1000)
    plt.close()

    # This normalizes the plot so that the area under the curve is = 1
    sns.distplot(AA, norm_hist=True, bins=100)
    plt.savefig('AA_' + b_file + "_" + str(iii) + '.png', dpi=1000)
    plt.close()

    sns.distplot(AB, norm_hist=True, bins=100)
    plt.savefig('AB_' + b_file + "_" + str(iii) + '.png', dpi=1000)
    plt.close()

    sns.distplot(BB, norm_hist=True, bins=100)
    plt.savefig('BB_' + b_file + "_" + str(iii) + '.png', dpi=1000)
    plt.close()

#    plt.hist(ALL, bins=100)
#    plt.savefig('ALL_' + b_file + "_" + str(iii) + '.png', dpi=1000)
#    plt.close()
#    
#    plt.hist(AA, bins=100)
#    plt.savefig('AA_' + b_file + "_" + str(iii) + '.png', dpi=1000)
#    plt.close()
#    
#    plt.hist(AB, bins=100)
#    plt.savefig('AB_' + b_file + "_" + str(iii) + '.png', dpi=1000)
#    plt.close()
#    
#    plt.hist(BB, bins=100)
#    plt.savefig('BB_' + b_file + "_" + str(iii) + '.png', dpi=1000)
#    plt.close()









