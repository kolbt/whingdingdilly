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
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.collections
from matplotlib import colors
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.tri as tri

def getDistance(point1, point2x, point2y):
    """Find the distance between two points"""
    distance = np.sqrt((point2x - point1[0])**2 + (point2y - point1[1])**2)
    return distance
    
def computeTauPerTstep(epsilon, mindt=0.000001):
    '''Read in epsilon, output tauBrownian per timestep'''
    if epsilon != 1.:
        mindt=0.00001
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
    
def findMinBins(lookDistance, xInd, yInd, maxInds, binSize):
    '''Try to pass fewer indicies to speed this process up'''
    lookN = int( lookDistance / sizeBin ) + 1
    # Entire square of indices
    left = xInd - lookN
    right = xInd + lookN
    bottom = yInd - lookN
    top = yInd + lookN
    # List to hold bins to check [ [x_inds], [y_inds] ]
    binsList = [ [], [] ]
    for i in range(left, right):
        horiz = (np.abs(x_ind - i) * binSize) - binSize
        for j in range(bottom, top):
            vert = (np.abs(y_ind - j) * binSize) - binSize
            binDist = np.sqrt(horiz**2 + vert**2)
            if binDist <= lookDistance:
                indX = i
                indY = j
                if i >= maxInds:
                    indX -= maxInds
                if j >= maxInds:
                    indY -= maxInds
                binsList[0].append(indX)
                binsList[1].append(indY)
    return binsList

# Get infile and open
inFile = str(sys.argv[1])
f = hoomd.open(name=inFile, mode='rb')
# Inside and outside activity from command line
peA = int(sys.argv[2])
peB = int(sys.argv[3])
parFrac = int(sys.argv[4])
eps = float(sys.argv[5])

start = 0                   # first frame to process
dumps = int(f.__len__())    # get number of timesteps dumped
end = dumps                 # final frame to process
start = dumps - 1

box_data = np.zeros((1), dtype=np.ndarray)  # box dimension holder
r_cut = 2**(1./6.)                          # potential cutoff
tauPerDT = computeTauPerTstep(epsilon=eps)  # brownian time per timestep

# Area of a particle (pi * r^2)
a_particle = np.pi * 0.25
# Search distance for local number density
lookDist = [1.0]
factor = 1.1

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
        nearNeigh = [0] * partNum

        # Easier accessors
        pos = snap.particles.position               # position
        xy = np.delete(pos, 2, 1)
        typ = snap.particles.typeid                 # type
        tst = snap.configuration.step               # timestep
        tst -= first_tstep                          # normalize by first timestep
        tst *= tauPerDT                             # convert to Brownian time
        
        plt.hexbin(pos[:,0], pos[:,1], gridsize=int(l_box/3.), linewidths=0.05)
        ax = plt.gca()
        ax.set_aspect('equal')
        ax.set_xlim(-hx_box, hx_box)
        ax.set_ylim(-hy_box, hy_box)
        ax.axis('off')
        plt.savefig('hexbin.png', dpi=1000)
        plt.close()

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
            # Get bins to search
            lookBins = findMinBins(max(lookDist), x_ind, y_ind, NBins, sizeBin)
            # Loop through all bins to search
            for w in xrange(len(lookBins[0])):
                wrapX = 0
                wrapY = 0
                # Finding negative horizontal indices
                if lookBins[0][w] < 0:
                    wrapX = l_box
                # Searching bins to the right, finding smaller index value
                if 0 <= lookBins[0][w] < (x_ind - lookNBins):
                    wrapX = -l_box

                # Finding negative vertical indices
                if lookBins[1][w] < 0:
                    wrapY = l_box
                # Searching above bins, finding smaller index value
                if 0 <= lookBins[1][w] < (y_ind - lookNBins):
                    wrapY = -l_box
                
                # Compute distance between particles: effective radius
                for b in range(0, len(binParts[lookBins[0][w]][lookBins[1][w]])):
                    # Reference index of particle we are comparing to
                    ref = binParts[lookBins[0][w]][lookBins[1][w]][b]
                    r = getDistance(pos[k],
                                    pos[ref][0] + wrapX,
                                    pos[ref][1] + wrapY)
                    r = round(r, 4)  # round value to 4 decimal places
                    if 0.1 < r < effSigma[k]:
                        effSigma[k] = r
                
                # Compute distance between particles: effective radius
                for b in range(0, len(binParts[lookBins[0][w]][lookBins[1][w]])):
                    # Reference index of particle we are comparing to
                    ref = binParts[lookBins[0][w]][lookBins[1][w]][b]
                    r = getDistance(pos[k],
                                    pos[ref][0] + wrapX,
                                    pos[ref][1] + wrapY)
                    r = round(r, 4)  # round value to 4 decimal places
                    if 0.1 < r < (factor * effSigma[k]):
                        nearNeigh[k] += 1
                        
        # Plot a countour plot with nearest neighbors
        ngrid = int(l_box/2.)
        xi = np.linspace(-hx_box, hx_box, ngrid)
        yi = np.linspace(-hx_box, hx_box, ngrid)
        triang = tri.Triangulation(pos[:,0], pos[:,1])
        interpolator = tri.LinearTriInterpolator(triang, nearNeigh)
        Xi, Yi = np.meshgrid(xi, yi)
        zi = interpolator(Xi, Yi)

        plt.contourf(xi, yi, zi)
        plt.savefig('myContour.png', dpi=1000)
