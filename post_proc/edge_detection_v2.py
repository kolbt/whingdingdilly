'''
#                           This is an 80 character line                       #
INTENT: Pick out the cluster edges so that you can analyze edge statistics

    1.  Cluster algorithm
    2.  Count nearest neighbors in largest cluster only
    --- Break here and plot this data ---
    3.  Look at contiguous set of low neighbor points (cluster them)
    4.  Perform usual statistics on this set (maybe the MSD?)

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

# This idea from the following blog post, Kevin Dwyer, YOU ARE THE MAN
# http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import shapely.geometry as geometry
from shapely.ops import cascaded_union, polygonize
from descartes import PolygonPatch
from scipy.spatial import Delaunay

def plot_polygon(polygon):
    fig = pl.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    margin = .3
    x_min, y_min, x_max, y_max = polygon.bounds
    ax.set_xlim([x_min-margin, x_max+margin])
    ax.set_ylim([y_min-margin, y_max+margin])
    patch = PolygonPatch(polygon, fc='#999999',
                         ec='#000000', fill=True,
                         zorder=-1)
    ax.add_patch(patch)
    return fig

# Load the file, get output name style
try:
    in_file = "pa" + str(pe_a) + \
              "_pb" + str(pe_b) + \
              "_xa" + str(part_perc_a) + \
              ".gsd"
    # File to write all data to
    out_file = "diam_pa" + str(pe_a) + \
               "_pb" + str(pe_b) + \
               "_xa" + str(part_perc_a)
    f = hoomd.open(name=in_file, mode='rb') # open gsd file with hoomd
except:
    in_file = "pa" + str(pe_a) + \
              "_pb" + str(pe_b) + \
              "_xa" + str(part_perc_a) + \
              "_ep1" + \
              ".gsd"
    # File to write all data to
    out_file = "diam_pa" + str(pe_a) + \
               "_pb" + str(pe_b) + \
               "_xa" + str(part_perc_a) + \
               "_ep1"
    f = hoomd.open(name=in_file, mode='rb')  # open gsd file with hoomd

dumps = f.__len__()                     # get number of timesteps dumped

start = 0       # gives first frame to read
start = dumps - 1
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
partNum = len(types[start])
part_A = int(partNum * part_frac_a)
part_B = partNum - part_A

# Feed data into freud analysis software
l_box = box_data[0]
h_box = l_box / 2.0
a_box = l_box * l_box
f_box = box.Box(Lx = l_box, Ly = l_box, is2D = True)    # make freud box
my_clust = cluster.Cluster(box = f_box,                 # init cluster class
                           rcut = 5.0)                # distance to search
c_props = cluster.ClusterProperties(box = f_box)        # compute properties

# Parameters for sorting dense from dilute
min_size_abso = 1000
min_size_perc = 0.05 * partNum  # minimum cluster size 5% of all particles
min_size = min_size_abso if min_size_abso < min_size_perc else min_size_perc

# Make the mesh
r_cut = 1.122
sizeBin = r_cut
nBins = int(l_box / sizeBin)
nBins += 1  # account for integer rounding

### LOOP THROUGH GSD FRAMES AND ANALYZE ###
for j in range(start, end):

    # Easier accessors
    pos = positions[j]
    # pos = np.delete(pos, 2, 1)
    typ = types[j]
    tst = timesteps[j]

    # Feed in the reduced particle set to the cluster algorithm
    my_clust.computeClusters(pos)  # feed in position
    ids = my_clust.getClusterIdx()  # get id of each cluster
    c_props.computeProperties(pos, ids)  # find cluster properties
    clust_size = c_props.getClusterSizes()  # find cluster sizes

    # Querry array, that tells whether cluster ID is of sufficient size
    lcIndex = 0     # particle index of largest cluster
    l_clust = 0     # any cluster is larger than this
    lcID = 0        # id of largest cluster

    for k in range(0, len(clust_size)):
        if clust_size[k] > l_clust:
            l_clust = clust_size[k]
            lcIndex = k
            lcID = ids[k]

    lc_posX = []
    lc_posY = []
    for k in range(0, len(ids)):
        if ids[k] == lcID:
            lc_posX.append(pos[k][0])
            lc_posY.append(pos[k][1])

    # Now I have an x,y,z array that is suitable for shapely
    lc_pos = np.zeros((len(lc_posX), 2), dtype=np.float64)
    lc_pos[:, 0] = lc_posX[:]
    lc_pos[:, 1] = lc_posY[:]
