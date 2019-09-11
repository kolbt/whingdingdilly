'''
#                           This is an 80 character line                       #
What does this file do?

1.) Read in .gsd file of particle positions
2.) Mesh the space
3.) Loop through tsteps and ...
3a.) Place all particles in appropriate mesh grid
3b.) Loop through all particles ...
3b.i.) Compute distance to every particle in adjacent grids
3.b.ii.) If distance is equal to particle diameter, increment neighbor count
3c.) Plot particle position colored by nearest neighbors (0-6)


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
    gsd_file = "pa" + str(pe_a) + \
                "_pb" + str(pe_b) + \
                "_xa" + str(part_perc_a) + \
                ".gsd"
    # File to write all data to
    out_file = "defects_pa" + str(pe_a) + \
               "_pb" + str(pe_b) + \
               "_xa" + str(part_perc_a) + \
               "_frame"
    f = hoomd.open(name=gsd_file, mode='rb') # open gsd file with hoomd
except:
    try:
        eps = str(sys.argv[6])
    except:
        eps = 1
        
    gsd_file = "pa" + str(pe_a) + \
                "_pb" + str(pe_b) + \
                "_xa" + str(part_perc_a) + \
                "_ep" + str(eps) + \
                ".gsd"
    # File to write all data to
    out_file = "defects_pa" + str(pe_a) + \
               "_pb" + str(pe_b) + \
               "_xa" + str(part_perc_a) + \
               "_ep" + str(eps) + \
               "_frame"
    f = hoomd.open(name=gsd_file, mode='rb')  # open gsd file with hoomd

dumps = int(f.__len__())                # get number of timesteps dumped

start = 0                       # gives first frame to read
end = dumps                     # gives last frame to read

positions = np.zeros((end), dtype=np.ndarray)       # array of positions
types = np.zeros((end), dtype=np.ndarray)           # particle types
box_data = np.zeros((1), dtype=np.ndarray)          # box dimensions
timesteps = np.zeros((end), dtype=np.float64)       # timesteps

# Get relevant data from long.gsd file
with hoomd.open(name=gsd_file, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    for iii in range(start, end):
        snap = t[iii]                               # snapshot of frame
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

# Make the mesh
r_cut = 1.122
# Make size of bin divisible by l_box:
divBox = round(l_box, 4)    # round l_box
if divBox - l_box < 0:
    divBox += 0.0001        # make sure rounded > actual

# Adjust sizeBin so that it divides into divBox:
convert = 100000.0
intBinSize = int(r_cut * convert)
intDivBox = int(divBox * convert)
while (intDivBox % intBinSize) != 0:
    intBinSize += 1
sizeBin = (intBinSize / convert)    # divisible bin size
nBins = int(divBox / sizeBin)       # must be an integer

# Enlarge the box to include the periodic images
buff = float(int(r_cut * 2.0) + 1)

# Image rendering options
drawBins = False
myCols = 'viridis'
#myCols = 'jet_r'
#myCols = 'gist_rainbow_r'

#for j in range(start, end):
for j in range(start, end):

    # Mesh array
    binParts = [[[] for b in range(nBins)] for a in range(nBins)]

    # Easier accessors
    pos = positions[j]
    # pos = np.delete(pos, 2, 1)
    typ = types[j]
    tst = timesteps[j]
    
    imageParts = []                                 # append (x, y) tuples
    # Duplicate replicate positions
    for k in range(0, partNum):
        
        # x-coordinate image creation
        if pos[k][0] + h_box < buff:                # against left box edge
            imageParts.append((pos[k][0] + l_box, pos[k][1]))
            # Create image of corners of periodic box
            if pos[k][1] + h_box < buff: # bottom left corner
                imageParts.append((pos[k][0] + l_box, pos[k][1] + l_box))
            if (pos[k][1] + h_box - l_box) > -buff: # top left corner
                imageParts.append((pos[k][0] + l_box, pos[k][1] - l_box))
    
        if (pos[k][0] + h_box - l_box) > -buff:     # against right box edge
            imageParts.append((pos[k][0] - l_box, pos[k][1]))
            # Create image of corners of periodic box
            if pos[k][1] + h_box < buff: # bottom right corner
                imageParts.append((pos[k][0] - l_box, pos[k][1] + l_box))
            if (pos[k][1] + h_box - l_box) > -buff: # top right corner
                imageParts.append((pos[k][0] - l_box, pos[k][1] - l_box))
        
        # y-coordinate image creation
        if pos[k][1] + h_box < buff:                # against bottom box edge
            imageParts.append((pos[k][0], pos[k][1] + l_box))
        
        if (pos[k][1] + h_box - l_box) > -buff:     # against top box edge
            imageParts.append((pos[k][0], pos[k][1] - l_box))

    # Put particles in their respective bins
    for k in range(0, partNum):
        # Get mesh indices
        tmp_posX = pos[k][0] + h_box
        tmp_posY = pos[k][1] + h_box
        x_ind = int(tmp_posX / sizeBin)
        y_ind = int(tmp_posY / sizeBin)
        # Append particle id to appropriate bin
        binParts[x_ind][y_ind].append(k)
    
    # Make an array that will hold the number of nearest neighbors
    near_neigh = [0] * partNum
    A_neigh = [0] * partNum
    B_neigh = [0] * partNum

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
                        near_neigh[k] += 1
                        if typ[ref] == 0:   # neighbor is A particle
                            A_neigh[k] += 1
                        else:               # neighbor is B particle
                            B_neigh[k] += 1

    # Plot position colored by neighbor number
    sz = 0.75
    scatter = plt.scatter(pos[:,0], pos[:,1],
                          c=near_neigh[:], cmap=myCols,
                          s=sz, edgecolors='none')
    xImages, yImages = zip(*imageParts)
    periodicIm = plt.scatter(xImages, yImages, c='#DCDCDC', s=sz, edgecolors='none')
    # Get colorbar
    plt.clim(0, 6)  # limit the colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label('# of Neighbors', rotation=270, labelpad=15)

    # Limits and ticks
    viewBuff = buff / 2.0
    plt.xlim(-h_box - viewBuff, h_box + viewBuff)
    plt.ylim(-h_box - viewBuff, h_box + viewBuff)
    plt.tick_params(axis='both', which='both',
                    bottom=False, top=False, left=False, right=False,
                    labelbottom=False, labeltop=False, labelleft=False, labelright=False)
                    
    pad = str(j).zfill(4)
    
    if drawBins:
        # Add the bins as vertical and horizontal lines:
        for binInd in xrange(nBins):
            coord = (r_cut * binInd) - h_box
            plt.axvline(x=coord, c='k', lw=1.0, zorder=0)
            plt.axhline(y=coord, c='k', lw=1.0, zorder=0)
    
    plt.tight_layout()
    plt.savefig(out_file + pad + '.png', dpi=500)
    plt.close()
