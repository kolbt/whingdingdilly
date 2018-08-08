'''
#                           This is an 80 character line                       #
PURE:  The intent of this file is to get out a few basic values as text files

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
a_box = l_box * l_box
f_box = box.Box(Lx = l_box, Ly = l_box, is2D = True)    # make freud box
my_clust = cluster.Cluster(box = f_box,                 # init cluster class
                              rcut = 1.0)               # distance to search
c_props = cluster.ClusterProperties(box = f_box)        # compute properties

# Parameters for sorting dense from dilute
min_size = 1000         # minimum cluster size
dist_min = 2.0          # bin density for avg d.p. density
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

# Create mesh
float_side = box_data[0] / 2.0
side = float((int(box_data[0])+1)) / 2
box_width = 5                                           # bin width
diameter = 1.0                                          # sigma
while int(side+1) % box_width != 0:                       # must be divisible
    side += 1                                           # make divisible by bin

spacer = int(side * 2 / (box_width * diameter))         # number of bins
mesh = np.zeros((spacer, spacer), dtype = np.int)       # array of each bin
reset_mesh = np.zeros_like(mesh)                        # zero out mesh

for iii in range(start, end):
    
    # Get data from arrays
    pos = positions[iii]
    typ = types[iii]
    tst = timesteps[iii]

    # Run freud computations
    my_clust.computeClusters(pos)           # feed in position
    ids = my_clust.getClusterIdx()          # get id of each cluster
    
    c_props.computeProperties(pos, ids)     # find cluster properties
    clust_size = c_props.getClusterSizes()  # find cluster sizes
    
    # Querry array, that tells whether cluster ID is of sufficient size
    q_clust = np.zeros((len(clust_size)), dtype = np.int)
    clust_num = 0   # number of clusters
    l_clust = 0     # largest cluster
    
    for jjj in range(0, len(clust_size)):
        if clust_size[jjj] > min_size:
            q_clust[jjj] = 1
            clust_num += 1
        if clust_size[jjj] > l_clust:
            l_clust = clust_size[jjj]

    # Values to write out after computations are complete
    dense_num = 0   # number of particles in dense phase
    dense_A = 0     # number of A-particles in dense phase
    dense_B = 0     # number of B-particles in dense phase
    gas_num = 0     # number of particles in gas phase
    gas_A = 0       # number of A-particles in gas phase
    gas_B = 0       # number of B-particles in gas phase
    dp_density = 0  # density of dense phase
    gp_density = 0  # density of gas phase
    mcs = 0         # mean cluster size

    for jjj in range(0, part_num):
        
        # Get the index of the mesh the particle belongs in
        loc_x = int((pos[jjj][0] + float_side) / (box_width * diameter))
        loc_y = int((pos[jjj][1] + float_side) / (box_width * diameter))
        
        # Add the particle to it's bin
        mesh[loc_x][loc_y] += 1
        
        if q_clust[ids[jjj]] == 1:          # it's in the dense phase
            dense_num += 1
            if typ[jjj] == 0:
                dense_A += 1
            else:
                dense_B += 1
        else:                               # it's in the gas phase
            gas_num += 1
            if typ[jjj] == 0:
                gas_A += 1
            else:
                gas_B += 1

    for jjj in range(0, part_num):
        if q_clust[ids[jjj]] == 1:          # is in dense phase
            # Indices of particle's bin
            loc_x = int((pos[jjj][0] + float_side) / (box_width * diameter))
            loc_y = int((pos[jjj][1] + float_side) / (box_width * diameter))
            # Add it and surrounding particles to dp_density
            dp_density += mesh[loc_x][loc_y]

    if dense_num != 0:
        dp_density /= float(dense_num)                  # avg number in each bin
        dp_density /= (float(box_width*diameter)**2)    # normalize by bin area
        mcs = dense_num / clust_num                     # compute mean cluster size

    a_clust = 0.0
    if dp_density != 0:
        a_clust = float(dense_num) / float(dp_density)  # area of cluster
    a_gas = float(a_box) - a_clust                      # area of gas
    gp_density = gas_num / a_gas                        # number density in gas

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
            '{0:.2f}'.format(dp_density).center(10) + ' ' +\
            '{0:.2f}'.format(gp_density).center(10) + ' ' +\
            '{0:.2f}'.format(a_clust).center(10) + ' ' +\
            '{0:.2f}'.format(a_gas).center(10) + ' ' +\
            str(mcs).center(10) + '\n')
    f.close()

    mesh[:] = 0     # zero out the mesh
