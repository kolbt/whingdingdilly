'''
#                           This is an 80 character line                       #
What does this file do?
(Reads single argument, .gsd file name)

1.) Read in .gsd file of particle positions
2.) Mesh the space
3.) Loop through tsteps and ...
3a.) Place all particles in appropriate mesh grid
3b.) Loop through all particles ...
3b.i.) Compute distance to every particle in adjacent grids
3.b.ii.) If distance is less than LJ cutoff, store as effective diameter
3c.) Plot particle with effective diameter as patch
4.) Generate movie from frames

'''

import sys
import os

sys.path.append("/Users/kolbt/Desktop/compiled/gsd/build")
sys.path.append("/nas/longleaf/home/kolbt/programs/gsd/build")
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
import matplotlib.collections

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
    
# Grab files
infile = str(sys.argv[1])                               # gsd file
print("Creating images for: {}").format(infile)
file_name = os.path.basename(infile)
outfile, file_extension = os.path.splitext(file_name)   # get base name
out = outfile + "_frame_"

f = hoomd.open(name=infile, mode='rb')  # open gsd file with hoomd

dumps = int(f.__len__())                # get number of timesteps dumped

start = 0                       # gives first frame to read
end = dumps                     # gives last frame to read

positions = np.zeros((end), dtype=np.ndarray)       # array of positions
types = np.zeros((end), dtype=np.ndarray)           # particle types
box_data = np.zeros((1), dtype=np.ndarray)          # box dimensions
timesteps = np.zeros((end), dtype=np.float64)       # timesteps

# Get relevant data from long.gsd file
with hoomd.open(name=infile, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    for iii in range(start, end):
        snap = t[iii]                               # snapshot of frame
        types[iii] = snap.particles.typeid          # get types
        positions[iii] = snap.particles.position    # get positions
        timesteps[iii] = snap.configuration.step    # get timestep

# Get simulation paramaters
partNum = len(positions[-1])
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

# Feed data into freud analysis software
l_box = box_data[0]
h_box = l_box / 2.0
a_box = l_box * l_box

# Make the mesh
r_cut = 2**(1./6.)

def getNBins(length, minSz=(2**(1./6.))):
    "Given box size, return number of bins"
    initGuess = int(length) + 1
    nBins = initGuess
    # This loop only exits on function return
    while True:
        if length / nBins > minSz:
            return nBins
        else:
            nBins -= 1
            
# Compute mesh
nBins = (getNBins(l_box, r_cut))
sizeBin = (l_box / nBins)

# Enlarge the box to include the periodic images
buff = float(int(r_cut * 2.0) + 1)

# Image rendering options
drawBins = False

#for j in range(start, end):
for j in range(start, end):

    # Mesh array
    binParts = [[[] for b in range(nBins)] for a in range(nBins)]

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
    
    # Make an array that will hold the effective diameter
    effSigma = [1.0] * partNum
    nearNeigh = [0] * partNum

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

                    # Store if shortest interparticle distance
                    if 0.1 < r < effSigma[k]:
                        effSigma[k] = r

    outDPI = 72.
    fig, ax = plt.subplots(figsize=(2000./outDPI, 2000./outDPI))
    xy = np.delete(pos, 2, 1)
    coll = matplotlib.collections.EllipseCollection(effSigma, effSigma,
                                                    np.zeros_like(effSigma),
                                                    offsets=xy, units='xy',
                                                    transOffset=ax.transData)
    ax.add_collection(coll)

    # Limits and ticks
    viewBuff = buff / 2.0
    viewBuff = 0.0
    ax.set_xlim(-h_box - viewBuff, h_box + viewBuff)
    ax.set_ylim(-h_box - viewBuff, h_box + viewBuff)
    ax.tick_params(axis='both', which='both',
                   bottom=False, top=False, left=False, right=False,
                   labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    pad = str(j).zfill(4)

    if drawBins:
        # Add the bins as vertical and horizontal lines:
        for binInd in xrange(nBins):
            coord = (sizeBin * binInd) - h_box
            plt.axvline(x=coord, c='k', lw=1.0, zorder=0)
            plt.axhline(y=coord, c='k', lw=1.0, zorder=0)

    # Change width of box edges
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(7.5)
        
    ax.set_aspect('equal')
    plt.tight_layout(pad=5.)
    plt.savefig(out + pad + '.png', dpi=outDPI)
    plt.close()
