'''
#                           This is an 80 character line                       #
PURE:  Read in simulation, output extrinsic data for each timestep
         
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
    out_file = "heatType_pa" + str(pe_a) + \
               "_pb" + str(pe_b) + \
               "_xa" + str(part_perc_a) + \
                "_fm"
    f = hoomd.open(name=long_file, mode='rb') # open gsd file with hoomd
    g = hoomd.open(name=short_file, mode='rb') # open gsd file with hoomd
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
    out_file = "heatType_pa" + str(pe_a) + \
               "_pb" + str(pe_b) + \
               "_xa" + str(part_perc_a) + \
               "_ep" + str(eps) +\
               "_fm"
    f = hoomd.open(name=long_file, mode='rb')  # open gsd file with hoomd
    g = hoomd.open(name=short_file, mode='rb')  # open gsd file with hoomd

dump_long = int(f.__len__())                # get number of timesteps dumped
dump_short = int(g.__len__())               # get number of timesteps dumped

start = 0                       # gives first frame to read
end = dump_long + dump_short    # gives last frame to read

positions = np.zeros((end), dtype=np.ndarray)       # array of positions
types = np.zeros((end), dtype=np.ndarray)           # particle types
box_data = np.zeros((1), dtype=np.ndarray)          # box dimensions
timesteps = np.zeros((end), dtype=np.float64)       # timesteps

# These make testing the script faster
stopOne = start + 2
stopTwo = dump_short + 2


# Get relevant data from short.gsd file
with hoomd.open(name=short_file, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    for iii in range(start, dump_short):
    # for iii in range(start, stopOne):
        snap = t[iii]                               # snapshot of frame
        types[iii] = snap.particles.typeid          # get types
        positions[iii] = snap.particles.position    # get positions
        timesteps[iii] = snap.configuration.step    # get timestep
# Get relevant data from long.gsd file
with hoomd.open(name=long_file, mode='rb') as t:
    snap = t[0]
    for iii in range(dump_short, end):
    # for iii in range(dump_short, stopTwo):
        snap = t[iii - dump_short]                  # snapshot of frame
        types[iii] = snap.particles.typeid          # get types
        positions[iii] = snap.particles.position    # get positions
        timesteps[iii] = snap.configuration.step    # get timestep

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
f_box = box.Box(Lx = l_box, Ly = l_box, is2D = True)    # make freud box
my_clust = cluster.Cluster(box = f_box,                 # init cluster class
                           rcut = 1.005)                 # distance to search
c_props = cluster.ClusterProperties(box = f_box)        # compute properties

# Parameters for sorting dense from dilute
min_size_abso = 1000
min_size_perc = 0.05 * partNum  # minimum cluster size 5% of all particles
min_size = min_size_abso if min_size_abso < min_size_perc else min_size_perc

# Make the mesh
r_cut = 1.122
sizeBin = 5
nBins = int(l_box / sizeBin)
nBins += 1  # account for integer rounding

for j in range(start, end):

    # Mesh array
    binParts = [[[] for b in range(nBins)] for a in range(nBins)]
    binSlow = np.zeros((nBins, nBins), dtype=np.int)
    binFast = np.zeros((nBins, nBins), dtype=np.int)
    binTotal = np.zeros((nBins, nBins), dtype=np.int)
    binMinus = np.zeros((nBins, nBins), dtype=np.int)

    # Easier accessors
    pos = positions[j]
    typ = types[j]
    tst = timesteps[j]

    # Loop through particles add to bins by type
    for k in range(0, partNum):
        # Get mesh indices
        tmp_posX = pos[k][0] + h_box
        tmp_posY = pos[k][1] + h_box
        x_ind = int(tmp_posX / sizeBin)
        y_ind = int(tmp_posY / sizeBin)
        binTotal[x_ind][y_ind] += 1
        if typ[k] == 0:
            binSlow[x_ind][y_ind] -= 1
            # Slow particles subtract
            binMinus[x_ind][y_ind] -= 1
        else:
            binFast[x_ind][y_ind] += 1
            # Fast particles add
            binMinus[x_ind][y_ind] += 1

    # Plot this timestep
    fig, axs = plt.subplots(2, 2)
    # im = axs[0, 0].imshow(binSlow.T, cmap='viridis', origin='lower', vmin=-25, vmax=25) # top left
    # im = axs[0, 1].imshow(binFast.T, cmap='viridis', origin='lower', vmin=-25, vmax=25) # top right
    # im = axs[1, 0].imshow(binTotal.T, cmap='viridis', origin='lower', vmin=-25, vmax=25) # bottom left
    # im = axs[1, 1].imshow(binMinus.T, cmap='viridis', origin='lower', vmin=-25, vmax=25) # bottom right
    im = axs[0, 0].imshow(binSlow.T, cmap='bwr', origin='lower', vmin=-25, vmax=25) # top left
    im = axs[0, 1].imshow(binFast.T, cmap='bwr', origin='lower', vmin=-25, vmax=25) # top right
    im = axs[1, 0].imshow(binTotal.T, cmap='bwr', origin='lower', vmin=-25, vmax=25) # bottom left
    im = axs[1, 1].imshow(binMinus.T, cmap='bwr', origin='lower', vmin=-25, vmax=25) # bottom right
    axs[0, 0].set_aspect('equal')
    axs[0, 1].set_aspect('equal')
    axs[1, 0].set_aspect('equal')
    axs[1, 1].set_aspect('equal')
    axs[0, 0].set_title('Slow only')
    axs[0, 1].set_title('Fast only')
    axs[1, 0].set_title('Total')
    axs[1, 1].set_title('Fast - Slow')
    axs[0, 0].get_xaxis().set_ticks([])
    axs[0, 0].get_yaxis().set_ticks([])
    axs[0, 1].get_xaxis().set_ticks([])
    axs[0, 1].get_yaxis().set_ticks([])
    axs[1, 0].get_xaxis().set_ticks([])
    axs[1, 0].get_yaxis().set_ticks([])
    axs[1, 1].get_xaxis().set_ticks([])
    axs[1, 1].get_yaxis().set_ticks([])
    fig.colorbar(im, ax=axs.ravel().tolist(), orientation='vertical')
    plt.savefig(out_file +
                str(j) +
                '.png', dpi=1000)
    plt.close()
