'''
#                           This is an 80 character line                       #
Read in:
    -file name
    -bin size
Output (for each timestep):
    -histogram of local density (number density, area fraction)
    -need to process this into text file that gives just the one-two-three?
'''

import sys

# Run locally
sys.path.append('/Users/kolbt/Desktop/compiled/hoomd-blue/build')
sys.path.append('/Users/kolbt/Desktop/compiled/gsd/build')
# Run on the cpu
sys.path.append('/nas/longleaf/home/kolbt/programs/cpu-hoomd/hoomd-blue/build')
# Run on the gpu
sys.path.append('/nas/longleaf/home/kolbt/programs/hoomd_2.2.1/hoomd-blue/build')
sys.path.append('/nas/longleaf/home/kolbt/programs/gsd/build')

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
import math
import random
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.collections
from matplotlib import colors
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
#plt.rcParams.update({'font.size': 4})
#plt.rcParams['axes.linewidth'] = 0.5
#plt.rcParams["xtick.major.size"] = 1
#plt.rcParams["xtick.major.width"] = 0.5
#plt.rcParams["xtick.minor.size"] = 1
#plt.rcParams["xtick.minor.width"] = 0.5
#plt.rcParams["ytick.major.size"] = 1
#plt.rcParams["ytick.major.width"] = 0.5
#plt.rcParams["ytick.minor.size"] = 1
#plt.rcParams["ytick.minor.width"] = 0.5

def getDistance(point1, point2x, point2y):
    """Find the distance between two points"""
    distance = np.sqrt((point2x - point1[0])**2 + (point2y - point1[1])**2)
    return distance
    
def computeTauPerTstep(epsilon, mindt=0.000001):
    '''Read in epsilon, output tauBrownian per timestep'''
    kBT = 1.0
    tstepPerTau = float(epsilon / (kBT * mindt))
    return 1. / tstepPerTau

def roundUp(n, decimals=0):
    '''Round up size of bins to account for floating point inaccuracy'''
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier
    
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

def findBins(lookN, currentInd, maxInds):
    maxInds -= 1
    left = currentInd - lookN
    right = currentInd + lookN
    binsList = []
    for i in range(left, right):
        ind = i
        if i > maxInds:
            ind -= maxInds
        binsList.append(ind)
    return binsList

# Get infile and open
inFile = str(sys.argv[1])
f = hoomd.open(name=inFile, mode='rb')

# Inside and outside activity from command line
peSlow = int(sys.argv[2])
peFast = int(sys.argv[3])

start = 0                   # first frame to process
dumps = int(f.__len__())    # get number of timesteps dumped
end = dumps                 # final frame to process
start = dumps - 2

box_data = np.zeros((1), dtype=np.ndarray)  # box dimension holder
r_cut = 2**(1./6.)                          # potential cutoff
tauPerDT = computeTauPerTstep(epsilon=1.)   # brownian time per timestep

lookDist = [1., 1.5]

# Set the colormap
myCols = plt.cm.viridis
fast = '#d8b365'
slow = '#5ab4ac'
colorsList = [slow, fast]
my_cmap = colors.ListedColormap(colorsList)

# Access and read .gsd data
with hoomd.open(name=inFile, mode='rb') as t:
    snap = t[0]
    first_tstep = snap.configuration.step
    box_data = snap.configuration.box
    # Box dimensions
    x_box = box_data[0]
    y_box = box_data[1]
    l_box = x_box
    hx_box = x_box / 2.0
    hy_box = y_box / 2.0
    h_box = hx_box
    a_box = x_box * y_box
    # Get the number particles
    pos = snap.particles.position
    partNum = len(pos)

    # Compute each mesh
    NBins = getNBins(x_box, r_cut)
    sizeBin = roundUp((x_box / NBins), 6)
    print("Number of bins: {}").format(NBins)
    print("Size of bins: {}").format(sizeBin)
    # Number of bins to search in each direction
    lookNBins = int( max(lookDist) / sizeBin ) + 1

    for j in range(start, end):
        # Set the system snapshot
        snap = t[j]
        # Empyt mesh to store particle IDs
        binParts = [[[] for b in range(NBins)] for a in range(NBins)]
        
        # Defaults to diameter of 1
        effSigma = [1.] * partNum
        # Number of neighbors depends on search distance
        nearNeigh = []
        for k in xrange(len(lookDist)):
            # Defaults to zero neighbors
            nearNeigh.append([0] * partNum)

        # Easier accessors
        pos = snap.particles.position               # position
        xy = np.delete(pos, 2, 1)
        typ = snap.particles.typeid                 # type
        tst = snap.configuration.step               # timestep
        tst -= first_tstep                          # normalize by first timestep
        tst *= tauPerDT                             # convert to Brownian time

        # Put particles in their respective bins
        for k in range(0, partNum):
            # Convert position to be > 0 to place in list mesh
            tmp_posX = pos[k][0] + h_box
            tmp_posY = pos[k][1] + h_box
            x_ind = int(tmp_posX / sizeBin)
            y_ind = int(tmp_posY / sizeBin)
            # Append particle id to appropriate bin
            binParts[x_ind][y_ind].append(k)
        
        # Get the number density
        for k in xrange(partNum):
            # Get mesh indices
            tmp_posX = pos[k][0] + h_box
            tmp_posY = pos[k][1] + h_box
            x_ind = int(tmp_posX / sizeBin)
            y_ind = int(tmp_posY / sizeBin)
            # Number of bins to search in each direction
            lookBinsX = findBins(lookNBins, x_ind, NBins)
            lookBinsY = findBins(lookNBins, y_ind, NBins)
            
            # Loop through horizontal bins
            for h in xrange(len(lookBinsX)):
                wrapX = 0.
                # Finding negative indices
                if lookBinsX[h] < 0:
                    wrapX = l_box
                # Searching bins to the right, finding smaller index value
                if h > lookNBins and lookBinsX[h] < x_ind:
                    wrapX = -l_box
                
                # Loop through vertical bins
                for v in xrange(len(lookBinsY)):
                    wrapY = 0.
                    # Finding negative indices
                    if lookBinsY[v] < 0:
                        wrapY = l_box
                    # Searching above bins, finding smaller index value
                    if v > lookNBins and lookBinsY[v] < y_ind:
                        wrapY = -l_box
        
                    # Compute distance between particles: effective radius
                    for b in range(0, len(binParts[lookBinsX[h]][lookBinsY[v]])):
                        # Reference index of particle we are comparing to
                        ref = binParts[lookBinsX[h]][lookBinsY[v]][b]
                        r = getDistance(pos[k],
                                        pos[ref][0] + wrapX,
                                        pos[ref][1] + wrapY)
                        r = round(r, 4)  # round value to 4 decimal places
                        for z in xrange(len(lookDist)):
                            if 0.1 < r < lookDist[z]:
                                nearNeigh[z][k] += 1
                        
        print(nearNeigh[0])
        
# This is VERY sensitive to the search distance
# r=1 -> 30s
# r=2 -> 1m53s
# r=3 -> 4m13s
