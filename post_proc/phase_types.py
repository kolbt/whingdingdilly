'''
#                           This is an 80 character line                       #
PURPOSE:  The intent of this file is to get out a few basic values as text files

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

in_file = "pa"+str(pe_a)+\
"_pb"+str(pe_b)+\
"_xa"+str(part_perc_a)+\
".gsd"

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
        positions[iii] = snap.particles.positions   # get positions
        timesteps[iii] = snap.configuration.step    # get timestep

timesteps -= timesteps[0]       # get rid of brownian run time

# Get number of each type of particle
part_num = len(types[0])
part_A = int(part_num * part_frac_a)
part_B = part_num - part_A

# Feed data into freud analysis software
l_box = box_data[0]
f_box = box.Box(Lx = l_box, Ly = l_box, is2D = True)    # make freud box
my_clust = cluster.Cluster(box = f_box,                 # init cluster class
                              rcut = 1.0)               # distance to search
c_props = cluster.ClusterProperties(box = f_box)        # compute properties

# Parameters for sorting dense from dilute
min_size = 1000

for iii in range(start, end):
    
    # Get data from arrays
    pos = positions[iii]
    typ = types[iii]
    tst = timesteps[iii]

    # Run freud computations
    my_clust.computeClusters(pos)           # feed in position
    ids = my_clust.getClusterIdx()          # get id of each cluster
    sort_id = np.sort(ids)                  # sorted ids, low-high
    
    c_props.computeProperties(pos, ids)     # find cluster properties
    clust_size = c_props.getClusterSizes()  # find cluster sizes
    
    # Querry array, that tells whether cluster ID is of sufficient size
    q_clust = np.zeros((len(clust_size)), dtype = np.int)
    clust_num = 0   # number of clusters
    
    for jjj in range(0, len(clust_size)):
        if clust_size[jjj] > size_min:
            q_clust[jjj] = 1
            clust_num += 1

    
    # Values to write out after computations are complete
    dense_num = 0   # number of particles in dense phase
    dense_A = 0     # number of A-particles in dense phase
    dense_B = 0     # number of B-particles in dense phase
    gas_num = 0     # number of particles in gas phase
    gas_A = 0       # number of A-particles in gas phase
    gas_B = 0       # number of B-particles in gas phase

    for jjj in range(0, part_num):
        if q_clust[ids[jjj]] == 1:          # it's in the dense phase
            dense_num
        else:                               # it's in the gas phase



# CLUSTER_COM GIVES CENTER OF MASS!!!! YES!!!!



















































