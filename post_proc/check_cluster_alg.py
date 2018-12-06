'''
#                           This is an 80 character line                       #
Summary:  Read in simulation gsds, output image of clustered particles
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
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def computeR(part1, part2):
    """Computes distance"""
    return np.sqrt(((part2[0]-part1[0])**2)+((part2[1]-part1[1])**2))

def computeA(diameter):
    """Computes area of circle"""
    radius = diameter / 2.0
    return np.pi * (radius**2)

def getDistance(point1, point2x, point2y):
    """Find the distance between two points"""
    distance = np.sqrt((point2x - point1[0])**2 + (point2y - point1[1])**2)
    return distance

def slowSort(array):
    """Sort an array the slow (but certain) way"""
    cpy = np.copy(array)
    ind = np.arange(0, len(array))
    for i in xrange(len(cpy)):
        for j in xrange(len(cpy)):
            if cpy[i] > cpy[j] and i < j:
                # Swap the copy array values
                tmp = cpy[i]
                cpy[i] = cpy[j]
                cpy[j] = tmp
                # Swap the corresponding indices
                tmp = ind[i]
                ind[i] = ind[j]
                ind[j] = tmp
    return ind

def indSort(arr1, arr2):
    """Take sorted index array, use to sort array"""
    # arr1 is array to sort
    # arr2 is index array
    cpy = np.copy(arr1)
    for i in xrange(len(arr1)):
        arr1[i] = cpy[arr2[i]]

def chkSort(array):
    """Make sure sort actually did its job"""
    for i in xrange(len(array)-2):
        if array[i] > array[i+1]:
            print("{} is not greater than {} for indices=({},{})").format(array[i+1], array[i], i, i+1)
            return False
    return True


try:
    # This is for the long timescale data
    long_file = "pa" + str(pe_a) + \
                "_pb" + str(pe_b) + \
                "_xa" + str(part_perc_a) + \
                ".gsd"
    # This is for the fine timescale data
    short_file = "log_pa" + str(pe_a) + \
                 "_pb" + str(pe_b) + \
                 "_xa" + str(part_perc_a) + \
                 ".gsd"
    # File to write all data to
    out_file = "pa" + str(pe_a) + \
               "_pb" + str(pe_b) + \
               "_xa" + str(part_perc_a) + \
               "_frame"
    f = hoomd.open(name=long_file, mode='rb')  # open gsd file with hoomd
    g = hoomd.open(name=short_file, mode='rb')  # open gsd file with hoomd
except:
    try:
        eps = str(sys.argv[6])
    except:
        eps = 1

    long_file = "pa" + str(pe_a) + \
                "_pb" + str(pe_b) + \
                "_xa" + str(part_perc_a) + \
                "_ep" + str(eps) + \
                ".gsd"
    short_file = "log_pa" + str(pe_a) + \
                 "_pb" + str(pe_b) + \
                 "_xa" + str(part_perc_a) + \
                 "_ep" + str(eps) + \
                 ".gsd"
    # File to write all data to
    out_file = "pa" + str(pe_a) + \
               "_pb" + str(pe_b) + \
               "_xa" + str(part_perc_a) + \
               "_ep" + str(eps) + \
               "_frame"
    f = hoomd.open(name=long_file, mode='rb')  # open gsd file with hoomd
    g = hoomd.open(name=short_file, mode='rb')  # open gsd file with hoomd

dump_long = int(f.__len__())  # get number of timesteps dumped
dump_short = int(g.__len__())  # get number of timesteps dumped

start = 0  # gives first frame to read
end = dump_long + dump_short  # gives last frame to read

positions = np.zeros((end), dtype=np.ndarray)  # array of positions
types = np.zeros((end), dtype=np.ndarray)  # particle types
box_data = np.zeros((1), dtype=np.ndarray)  # box dimensions
timesteps = np.zeros((end), dtype=np.float64)  # timesteps

# Get relevant data from short.gsd file
with hoomd.open(name=short_file, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    for iii in range(start, dump_short):
        snap = t[iii]  # snapshot of frame
        types[iii] = snap.particles.typeid  # get types
        positions[iii] = snap.particles.position  # get positions
        timesteps[iii] = snap.configuration.step  # get timestep
# Get relevant data from long.gsd file
with hoomd.open(name=long_file, mode='rb') as t:
    snap = t[0]
    for iii in range(dump_short, end):
        snap = t[iii - dump_short]  # snapshot of frame
        types[iii] = snap.particles.typeid  # get types
        positions[iii] = snap.particles.position  # get positions
        timesteps[iii] = snap.configuration.step  # get timestep

# Get index order for chronological timestep sorting
newInd = slowSort(timesteps)
# Use these indexes to reorder other arrays
indSort(timesteps, newInd)
indSort(positions, newInd)
indSort(types, newInd)

if chkSort(timesteps):
    print("Array succesfully sorted")
else:
    print("Array not sorted")
timesteps -= timesteps[0]

# Get number of each type of particle
partNum = len(types[start])
part_A = int(partNum * part_frac_a)
part_B = partNum - part_A

# Feed data into freud analysis software
l_box = box_data[0]
h_box = l_box / 2.0
a_box = l_box * l_box
f_box = box.Box(Lx=l_box, Ly=l_box, is2D=True)  # make freud box
my_clust = cluster.Cluster(box=f_box,  # init cluster class
                           rcut=1.005)  # distance to search
c_props = cluster.ClusterProperties(box=f_box)  # compute properties

# Parameters for sorting dense from dilute
min_size_abso = 1000
min_size_perc = 0.05 * partNum  # minimum cluster size 5% of all particles
min_size = min_size_abso if min_size_abso < min_size_perc else min_size_perc

den_pinkA = '#ff00ff'
gas_pinkA = '#ff7fff'
den_greenB = '#2a621c'
gas_greenB = '#6bf648'

for j in range(start, end):

    # Easier accessors
    pos = positions[j]
    typ = types[j]
    tst = timesteps[j]

    # Run freud computations
    my_clust.computeClusters(pos)           # feed in position
    ids = my_clust.getClusterIdx()          # get id of each cluster
    c_props.computeProperties(pos, ids)     # find cluster properties
    clust_size = c_props.getClusterSizes()  # find cluster sizes

    # Querry array, that tells whether cluster ID is of sufficient size
    q_clust = np.zeros((len(clust_size)), dtype=np.int)
    clust_num = 0   # number of clusters
    lcIndex = 0     # id of largest cluster

    # This sorts clusters that are too small
    l_clust = 1
    for k in range(0, len(clust_size)):
        if clust_size[k] > min_size:
            q_clust[k] = 1
            clust_num += 1
        if clust_size[k] > l_clust:
            l_clust = clust_size[k]
            lcIndex = k

    # Various lists for quick plot
    gas_Ax = []
    gas_Ay = []
    gas_Bx = []
    gas_By = []
    lc_Ax = []
    lc_Ay = []
    lc_Bx = []
    lc_By = []
    den_Ax = []
    den_Ay = []
    den_Bx = []
    den_By = []

    # Put points into lists
    for k in range(0, partNum):
        # This is in gas phase
        if q_clust[ids[k]] == 0:
            # A
            if typ[k] == 0:
                gas_Ax.append(pos[k][0])
                gas_Ay.append(pos[k][1])
            # B
            else:
                gas_Bx.append(pos[k][0])
                gas_By.append(pos[k][1])
        # This is in dense phase
        else:
            # Part of largest cluster
            if ids[k] == lcIndex:
                if typ[k] == 0:
                    lc_Ax.append(pos[k][0])
                    lc_Ay.append(pos[k][1])
                else:
                    lc_Bx.append(pos[k][0])
                    lc_By.append(pos[k][1])
            # Other dense phase
            else:
                if typ[k] == 0:
                    den_Ax.append(pos[k][0])
                    den_Ay.append(pos[k][1])
                else:
                    den_Bx.append(pos[k][0])
                    den_By.append(pos[k][1])

    # It's easier to plot this way?
    plt.scatter(gas_Ax, gas_Ay, s=0.75, facecolor=gas_pinkA, edgecolor='none')
    plt.scatter(gas_Bx, gas_By, s=0.75, facecolor=gas_greenB, edgecolor='none')
    plt.scatter(lc_Ax, lc_Ay, s=0.75, facecolor=den_pinkA, edgecolor='none')
    plt.scatter(lc_Bx, lc_By, s=0.75, facecolor=den_greenB, edgecolor='none')
    plt.scatter(den_Ax, den_Ay, s=0.75, facecolor=den_pinkA, edgecolor='none')
    plt.scatter(den_Bx, den_By, s=0.75, facecolor=den_greenB, edgecolor='none')

    plt.xlim(-h_box, h_box)
    plt.ylim(-h_box, h_box)
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()
    plt.axes().set_aspect('equal')
    plt.savefig(out_file + str(j) + '.png', dpi=1000)
    plt.close()
