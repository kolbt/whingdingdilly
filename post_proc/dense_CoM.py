'''
#                           This is an 80 character line                       #
PURPOSE:  This script will take many forms in it's life.  First I want to out a
          .png of my system with particles belonging to the largest cluster 
          highlighted.  I will plot the center of mass on top of this .png and
          create a movie to see how this script does in guessing the CoM.
          
          Should it work well, I'll come back, edit, and output the CoM MSD
          and the average distance from the CoM by type.

          Output :  mvy_pa$_pb$_xa$_$.png
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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

# File to read from
in_file = "pa"+str(pe_a)+\
"_pb"+str(pe_b)+\
"_xa"+str(part_perc_a)+\
".gsd"
# File to write all data to
out_file = "mvy_pa"+str(pe_a)+\
"_pb"+str(pe_b)+\
"_xa"+str(part_perc_a)+\
"_"

f = hoomd.open(name=in_file, mode='rb') # open gsd file with hoomd
dumps = f.__len__()                     # get number of timesteps dumped

start = 0       # gives first frame to read
end = dumps     # gives last frame to read

positions = np.zeros((end), dtype=np.ndarray)       # array of positions
types = np.zeros((end), dtype=np.ndarray)           # particle types
box_data = np.zeros((1), dtype=np.ndarray)          # box dimensions
timesteps = np.zeros((end), dtype=np.float64)       # timesteps
CoMs = np.zeros((end), dtype=np.ndarray)            # holds CoM

# Get relevant data from .gsd file
with hoomd.open(name=in_file, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    for iii in range(start, end):
        snap = t[iii]                               # snapshot of frame
        types[iii] = snap.particles.typeid          # get types
        positions[iii] = snap.particles.position    # get positions
        timesteps[iii] = snap.configuration.step    # get timestep

timesteps -= timesteps[start]   # get rid of brownian run time

# Get number of each type of particle
part_num = len(types[start])
part_A = int(part_num * part_frac_a)
part_B = part_num - part_A

# Feed data into freud analysis software
l_box = box_data[0]
f_box = box.Box(Lx = l_box, Ly = l_box, is2D = True)    # make freud box
my_clust = cluster.Cluster(box = f_box,                 # init cluster class
                              rcut = 1.0)               # distance to search
c_props = cluster.ClusterProperties(box = f_box)        # compute properties

# Parameters for sorting dense from dilute
min_size = 1000                 # minimum cluster size
min = -(float(box_data[0]/2))   # smallest axes value, plotting
max = (float(box_data[0]/2))    # largest axes value, plotting

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
    
    # Find largest cluster, store it's size and ID
    l_clust = 0     # largest cluster
    l_id = 0        # id of largest cluster
    l_flag = 0      # switched if a cluster is sufficiently big
    
    for jjj in range(0, len(clust_size)):
        if clust_size[jjj] > min_size and clust_size[jjj] > l_clust:
            l_clust = clust_size[jjj]
            l_id = ids[jjj]
            l_flag = 1

    # Write len(part_num) array where particle sin l_clust are 1's
    l_color = np.zeros((part_num), dtype=np.int)
    if l_flag == 1:
        for jjj in range(0, part_num):
            if ids[jjj] == l_id:
                l_color[jjj] = 1
        # Compute CoM of largest cluster
        getCoM = c_props.getClusterCOM()
        CoMs[iii] = getCoM[l_id]

    # Plot all points (using l_clust colormap)
    plt.scatter(pos[:, 0],
                pos[:, 1],
                s=0.25,
                c=l_color,
                edgecolors='none')
    # If there is a sufficient cluster, plot it's CoM
    if CoMs[iii][0] != 0:
        plt.scatter(CoMs[iii][0],
                    CoMs[iii][1],
                    s=4.0,
                    c='k',
                    edgecolors='none')

    plt.xlim(min, max)
    plt.ylim(min, max)
    plt.axes().get_xaxis().set_visible(False)
    plt.axes().get_yaxis().set_visible(False)
    plt.axes().set_aspect('equal')
    plt.savefig(out_file + str(iii) + '.png', dpi=1000)
    plt.close()
