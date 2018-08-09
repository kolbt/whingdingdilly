'''
#                           This is an 80 character line                       #
PURE:  Read in simulation, output extrinsic data for each timestep

                    File   :   Column 1    Column 2    Column 3    Column 4
     gas_pa#_pb#_xa#.txt   :   time        gas_A       gas_B       gas_total
   dense_pa#_pb#_xa#.txt   :   time        dense_A     dense_B     dense_total
  lclust_pa#_pb#_xa#.txt   :   time        l_clust
         
Designations Explained:
         
            gas_A : number of A-type particles in gas phase
            gas_B : number of B-type particles in gas phase
        gas_total : total number of gaseous partilces
          dense_A : number of A-type particles in dense phase
          dense_B : number of B-type particles in dense phase
      dense_total : total number of dense phase particles
          l_clust : largest individual cluster
              mcs : mean cluster size
              ALL : mode radius of all particles
               AA : mode radius of AA interactions
               AB : mode radius of AB interactions
               BB : mode radius of BB interactions
              lcA : largest cluster area
              tcA : total cluster area
       dp_density : density of dense phase
       gp_density : density of gas phase
'''

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

def computeR(part1, part2):
    return np.sqrt(((part2[0]-part1[0])**2)+((part2[1]-part1[1])**2))

# File to read from
# in_file = "pa"+str(pe_a)+\
# "_pb"+str(pe_b)+\
# "_xa"+str(part_perc_a)+\
# "_ep1"+\
# ".gsd"

# Epsilon not in filename
in_file = "pa"+str(pe_a)+\
"_pb"+str(pe_b)+\
"_xa"+str(part_perc_a)+\
".gsd"

# File to write all data to
all_file = "all_pa"+str(pe_a)+\
"_pb"+str(pe_b)+\
"_xa"+str(part_perc_a)+\
"_ep1"+\
".txt"

f = hoomd.open(name=in_file, mode='rb') # open gsd file with hoomd
dumps = f.__len__()                     # get number of timesteps dumped

start = 0       # gives first frame to read
end = dumps     # gives last frame to read

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
my_clust = cluster.Cluster(box = f_box,                 # init cluster class
                              rcut = 1.0)               # distance to search
c_props = cluster.ClusterProperties(box = f_box)        # compute properties

# Parameters for sorting dense from dilute
min_size_abso = 1000
min_size_perc = 0.05 * part_num  # minimum cluster size 5% of all particles
min_size = min_size_abso if min_size_abso < min_size_perc else min_size_perc

f = open(all_file, 'w') # write file headings
f.write('Timestep'.center(10) + ' ' +\
        'Gas_A'.center(10) + ' ' +\
        'Gas_B'.center(10) + ' ' +\
        'Gas_tot'.center(10) + ' ' +\
        'Dense_A'.center(10) + ' ' +\
        'Dense_B'.center(10) + ' ' +\
        'Dense_tot'.center(10) + ' ' +\
        'Lg_clust'.center(10) + ' ' +\
        'MCS'.center(10) + ' ' +\
        'sigALL'.center(10) + ' ' +\
        'sigAA'.center(10) + ' ' +\
        'sigAB'.center(10) + ' ' +\
        'sigBB'.center(10) + ' ' +\
        'phiEff'.center(10) + ' ' +\
        'lg_clustA'.center(10) + ' ' +\
        'tot_clustA'.center(10) + ' ' +\
        'DP_density'.center(10) + ' ' +\
        'GP_density'.center(10) + '\n')
f.close()

# Make the mesh
r_cut = 1.122
sizeBin = r_cut
nBins = int(l_box / sizeBin)
nBins += 1  # account for integer rounding

for j in range(start, end):

    # Mesh array
    binParts = [[[] for b in range(nBins)] for a in range(nBins)]

    # Easier accessors
    pos = positions[j]
    pos = np.delete(pos, 2, 1)
    typ = types[j]

    # Lists for effective diameter
    ALL = []
    AA = []
    AB = []
    BB = []

    # Put particles in their respective bins
    for k in range(0, partNum):
        # Get mesh indices
        tmp_posX = pos[k][0] + h_box
        tmp_posY = pos[k][1] + h_box
        x_ind = int(tmp_posX / sizeBin)
        y_ind = int(tmp_posY / sizeBin)
        # Append particle id to appropriate bin
        binParts[x_ind][y_ind].append(k)

    # Compute distance, each pair will be counted twice
    for k in range(0, partNum):
        # Get mesh indices
        tmp_posX = pos[k][0] + h_box
        tmp_posY = pos[k][1] + h_box
        x_ind = int(tmp_posX / sizeBin)
        y_ind = int(tmp_posY / sizeBin)
        # Get index of surrounding bins
        l_bin = x_ind - 1  # index of left bins
        r_bin = x_ind + 1  # index of right bins
        b_bin = y_ind - 1  # index of bottom bins
        t_bin = y_ind + 1  # index of top bins
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
                    r = getDistance(pos[k],
                                    pos[ref][0] + wrapX,
                                    pos[ref][1] + wrapY)
                    r = round(r, 4)  # round value to 4 decimal places

                    # If LJ potential is on, store into a list (omit self)
                    if 0.1 < r <= r_cut:
                        ALL.append(r)  # All particles
                        if typ[k] == 0 and typ[ref] == 0:  # AA distance
                            AA.append(r)
                        elif typ[k] == 1 and typ[ref] == 1:  # BB distance
                            BB.append(r)
                        else:  # AB distance
                            AB.append(r)

    # You can compute the density from the intended radius and cluster area

    # Values have been set, write to text files
    f = open(all_file, 'a')
    f.write(('{0:.2f}'.format(tst)).center(10) + ' ' +\
            str(gas_A).center(10) + ' ' +\
            str(gas_B).center(10) + ' ' +\
            str(gas_num).center(10) + ' ' +\
            str(dense_A).center(10) + ' ' +\
            str(dense_B).center(10) + ' ' +\
            str(dense_num).center(10) + ' ' +\
            str(l_clust).center(10) + ' ' +\
            str(mcs).center(10) + ' ' +\
            str(modeALL).center(10) + ' ' +\
            str(modeAA).center(10) + ' ' +\
            str(modeAB).center(10) + ' ' +\
            str(modeBB).center(10) + ' ' +\
            str(phiEff).center(10) + ' ' +\
            str(lcA).center(10) + ' ' +\
            str(tcA).center(10) + ' ' +\
            '{0:.2f}'.format(dp_density).center(10) + ' ' +\
            '{0:.2f}'.format(gp_density).center(10) + '\n')
    f.close()
