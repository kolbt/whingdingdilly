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

import numpy as np
import matplotlib.pyplot as plt

def getDistance(point1, point2x, point2y):
    """Find the distance between two points"""
    distance = np.sqrt((point2x - point1[0])**2 + (point2y - point1[1])**2)
    return distance

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

    # Mesh array
    binParts = [[[] for b in range(nBins)] for a in range(nBins)]

    # Neighbor array
    neighbors = np.zeros(partNum, dtype=np.int)

    # Easier accessors
    pos = positions[j]
    # pos = np.delete(pos, 2, 1)
    typ = types[j]
    tst = timesteps[j]

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

                    # If particles are touching, they are neighbors
                    if 0.1 < r <= 1.0:
                        neighbors[k] += 1

    # # Sanity check -- plot heatmap position w/ neighbors overlaid
    # plt.scatter(pos[:, 0], pos[:, 1], c=neighbors, s=0.5, edgecolors='none')
    # plt.xticks(())
    # plt.yticks(())
    # plt.xlim(-h_box, h_box)
    # plt.ylim(-h_box, h_box)
    # plt.colorbar()
    # plt.title('Number of Nearest Neighbors')
    # plt.savefig(out_file + '_frame' + str(j) + '.png', dpi=1000)

    # Let's keep everything with 3,4,5 nearest neighbors and plot it
    edgeX = []
    edgeY = []
    neigh = []

    nMin = 2
    nMax = 3

    for k in range(0, partNum):
        if (nMin <= neighbors[k] and neighbors[k] <= nMax):
            edgeX.append(pos[k][0])
            edgeY.append(pos[k][1])
            neigh.append(neighbors[k])

    neigh_pos = np.zeros((len(edgeX), 3), dtype=np.float64)
    neigh_pos[:, 0] = edgeX[:]
    neigh_pos[:, 1] = edgeY[:]

    # # Sanity check -- plot only particles with specific number of neighbors
    # plt.scatter(edgeX, edgeY, c=neigh, s=0.5, edgecolors='none')
    # plt.xticks(())
    # plt.yticks(())
    # plt.xlim(-h_box, h_box)
    # plt.ylim(-h_box, h_box)
    # plt.colorbar()
    # plt.title('3 < Neighbors < 5')
    # plt.savefig('3-5neigh_' + out_file + '_frame' + str(j) + '.png', dpi=1000)

    # Feed in the reduced particle set to the cluster algorithm
    my_clust.computeClusters(neigh_pos)  # feed in position
    ids = my_clust.getClusterIdx()  # get id of each cluster
    c_props.computeProperties(neigh_pos, ids)  # find cluster properties
    clust_size = c_props.getClusterSizes()  # find cluster sizes

    # Querry array, that tells whether cluster ID is of sufficient size
    q_clust = np.zeros((len(clust_size)), dtype=np.int)
    clust_num = 0   # number of clusters
    lcIndex = 0     # id of largest cluster
    lc_numA = 0     # number of A particles in dense phase
    lc_numB = 0     # number of B particles in dense phase
    l_clust = 0     # any cluster is larger than this

    for k in range(0, len(clust_size)):
        if clust_size[k] > min_size:
            q_clust[k] = 1
            clust_num += 1
        if clust_size[k] > l_clust:
            l_clust = clust_size[k]
            lcIndex = k
            lcID = ids[k]

    lc_posX = []
    lc_posY = []
    for k in range(0, len(ids)):
        if ids[k] == lcID:
            lc_posX.append(neigh_pos[k][0])
            lc_posY.append(neigh_pos[k][1])

    plt.scatter(lc_posX, lc_posY, s=0.5, edgecolors='none')
    plt.xticks(())
    plt.yticks(())
    plt.xlim(-h_box, h_box)
    plt.ylim(-h_box, h_box)
    plt.title('Largest edge cluster')
    plt.savefig('lcEdge_' + out_file + '_frame' + str(j) + '.png', dpi=1000)
