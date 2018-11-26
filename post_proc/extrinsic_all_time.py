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
    all_file = "diam_pa" + str(pe_a) + \
               "_pb" + str(pe_b) + \
               "_xa" + str(part_perc_a) + \
               ".txt"
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
    all_file = "diam_pa" + str(pe_a) + \
               "_pb" + str(pe_b) + \
               "_xa" + str(part_perc_a) + \
               "_ep" + str(eps) +\
               ".txt"
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

# Get relevant data from short.gsd file
with hoomd.open(name=short_file, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    for iii in range(start, dump_short):
        snap = t[iii]                               # snapshot of frame
        types[iii] = snap.particles.typeid          # get types
        positions[iii] = snap.particles.position    # get positions
        timesteps[iii] = snap.configuration.step    # get timestep
# Get relevant data from long.gsd file
with hoomd.open(name=long_file, mode='rb') as t:
    snap = t[0]
    for iii in range(dump_short, end):
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
                           rcut = 1.0)                  # distance to search
c_props = cluster.ClusterProperties(box = f_box)        # compute properties

# Parameters for sorting dense from dilute
min_size_abso = 1000
min_size_perc = 0.05 * partNum  # minimum cluster size 5% of all particles
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
        'LC_density'.center(10) + ' ' +\
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
    # pos = np.delete(pos, 2, 1)
    typ = types[j]
    tst = timesteps[j]

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

    # Now compute the mode radii for each type of interaction
    # ALL
    modeALL = stats.mode(ALL)
    modeALL = round(modeALL[0][0], 4)
    # fALL = computeLJForce(modeALL, epsHS[i])
    phiEff = modeALL * 0.6  # Effective area fraction phi=0.6

    # AA
    try:
        modeAA = stats.mode(AA)
        modeAA = round(modeAA[0][0], 4)
        # fAA = computeLJForce(modeAA, epsHS[i])
    except:
        modeAA = 0.0
        fAA = 0.0

    # AB
    try:
        modeAB = stats.mode(AB)
        modeAB = round(modeAB[0][0], 4)
        # fAB = computeLJForce(modeAB, epsHS[i])
    except:
        modeAB = 0.0
        fAB = 0.0

    # BB
    try:
        modeBB = stats.mode(BB)
        modeBB = round(modeBB[0][0], 4)
        # fBB = computeLJForce(modeBB, epsHS[i])
    except:
        modeBB = 0.0
        fBB = 0.0

    # Values to write out after computations are complete
    gas_A = 0       # number of A-particles in gas phase
    gas_B = 0       # number of B-particles in gas phase
    gas_num = 0     # number of particles in gas phase
    dense_A = 0     # number of A-particles in dense phase
    dense_B = 0     # number of B-particles in dense phase
    dense_num = 0       # number of particles in dense phase
    l_clust = 0         # largest cluster size
    mcs = 0             # mean cluster size
    lcA = 0.0           # largest cluster area
    tcA = 0.0           # total cluster area
    dp_density = 0.0    # density of dense phase
    gp_density = 0.0    # density of gas phase

    # Run freud computations
    my_clust.computeClusters(pos)           # feed in position
    ids = my_clust.getClusterIdx()          # get id of each cluster
    c_props.computeProperties(pos, ids)     # find cluster properties
    clust_size = c_props.getClusterSizes()  # find cluster sizes

    # Querry array, that tells whether cluster ID is of sufficient size
    q_clust = np.zeros((len(clust_size)), dtype=np.int)
    clust_num = 0   # number of clusters
    lcIndex = 0     # id of largest cluster
    lc_numA = 0     # number of A particles in dense phase
    lc_numB = 0     # number of B particles in dense phase

    for k in range(0, len(clust_size)):
        if clust_size[k] > min_size:
            q_clust[k] = 1
            clust_num += 1
        if clust_size[k] > l_clust:
            l_clust = clust_size[k]
            lcIndex = k

    for k in range(0, partNum):
        # Data pertaining to largest cluster
        if ids[k] == lcIndex:
            if typ[k] == 0:
                lc_numA += 1
            else:
                lc_numB += 1
        # This is in gas phase
        if q_clust[ids[k]] == 0:
            gas_num += 1
            if typ[k] == 0:
                gas_A += 1
            else:
                gas_B += 1
        # This is in dense phase
        else:
            dense_num += 1
            if typ[k] == 0:
                dense_A += 1
            else:
                dense_B += 1

    # Compute some things...
    if clust_num != 0:
        mcs = dense_num / clust_num
    lcA = (lc_numA * computeA(modeAA)) + (lc_numB * computeA(modeBB))
    tcA = (dense_A * computeA(modeAA)) + (dense_B * computeA(modeBB))
    # Number densities (number per unit area)
    if lcA != 0.0:
        lc_density = l_clust / lcA
    if tcA != 0.0:
        dp_density = dense_num / tcA
    gp_density = gas_num / (a_box - tcA)

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
            str(mcs).center(10) + ' ' + \
            '{0:.4f}'.format(modeALL).center(10) + ' ' + \
            '{0:.4f}'.format(modeAA).center(10) + ' ' + \
            '{0:.4f}'.format(modeAB).center(10) + ' ' + \
            '{0:.4f}'.format(modeBB).center(10) + ' ' + \
            '{0:.2f}'.format(phiEff).center(10) + ' ' + \
            '{0:.1f}'.format(lcA).center(10) + ' ' + \
            '{0:.1f}'.format(tcA).center(10) + ' ' +\
            '{0:.2f}'.format(lc_density).center(10) + ' ' +\
            '{0:.2f}'.format(dp_density).center(10) + ' ' +\
            '{0:.2f}'.format(gp_density).center(10) + '\n')
    f.close()
