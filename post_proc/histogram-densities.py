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
        horiz = (np.abs(x_ind - i) * binSize)
        for j in range(bottom, top):
            vert = (np.abs(y_ind - j) * binSize)
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
try:
    phi = float(sys.argv[6])
    intPhi = int(phi)
    phi /= 100.
except:
    phi = 0.6
    intPhi = 60

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
lookDist = np.arange(1., 5.5, 0.5)
#lookDist = [1., 2., 3., 3.5, 4.0, 4.5, 5.0]
# Area of each search distance
lookArea = []
# This gives the number of bins to plot
plotBins = []
# This gives the limits of the x-axis
minDense = []
maxDense = []
for z in xrange(len(lookDist)):
    area = np.pi * (lookDist[z]**2)
    lookArea.append(area)
    maxParts = int(area / a_particle) + 1
    plotBins.append(maxParts)
#    minDense.append(a_particle / area)
    minDense.append(0.)
    maxDense.append(2.0)

# Set the width of each bin in the histogram (in units of sigma)
histWidth = 0.005

# Write the low density and high density peaks to a text file
txtFile = 'phase_density_pa' + str(peA) +\
          '_pb' + str(peB) +\
          '_xa' + str(parFrac) +\
          '_phi' + str(intPhi) +\
          '.txt'
f = open(txtFile, 'w') # write file headings
f.write('Timestep'.center(10) + ' ')
for z in xrange(len(lookDist)):
    myString = 'Gas-r=' + str(lookDist[z])
    f.write(myString.center(10) + ' ')
    myString = 'Liq-r=' + str(lookDist[z])
    f.write(myString.center(10))
    if z < (len(lookDist) - 1):
        f.write(' ')
    else:
        f.write('\n')
f.close()

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
    # Let's get the simulation length
    snap = t[1]
    dt = snap.configuration.step - first_tstep
    snap = t[-1]
    last_tstep = snap.configuration.step
    tot_tstep = last_tstep - first_tstep
    space = tot_tstep / 10
    outFrames = np.arange(first_tstep, tot_tstep, space)
    # Make sure these timesteps exist
    for j in xrange(len(outFrames)):
        while (outFrames[j] - first_tstep) % dt != 0:
            outFrames[j] += 1

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
            # Defaults to one neighbor (itself)
            nearNeigh.append([1] * partNum)

        # Easier accessors
        pos = snap.particles.position               # position
        xy = np.delete(pos, 2, 1)
        typ = snap.particles.typeid                 # type
        tst = snap.configuration.step               # timestep
        tst -= first_tstep                          # normalize by first timestep
        tst *= tauPerDT                             # convert to Brownian time
        
        # Write the timestep to the text output
        f = open(txtFile, 'a')
        f.write('{0:.1f}'.format(tst).center(10) + ' ')
        f.close()

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
                    for z in xrange(len(lookDist)):
                        if 0.1 < r < lookDist[z]:
                            nearNeigh[z][k] += 1
                        
        # Now plot a histogram of the number density 0 -> HS
        pad = str(j).zfill(4)
        for z in xrange(len(lookDist)):
            # Number of particles * area of a particle / area of bin
            nearNeigh[z] = [ (y * a_particle / lookArea[z]) for y in nearNeigh[z] ]
            # Weight each particle equally
            weight = np.ones_like(nearNeigh[z]) / float(len(nearNeigh[z]))
            # Plot histogram
#            plt.hist(nearNeigh[z], weights=weight,
#                     bins=plotBins[z], range=(minDense[z], maxDense[z]),
#                     edgecolor='k',
#                     alpha=1.0,
#                     zorder=0)
                     
            # N is the count in each bin, bins is the lower-limit of the bin
            N, bins, patches = plt.hist(nearNeigh[z],
                                        bins=int(maxDense[z] / histWidth),
                                        range=(minDense[z], maxDense[z]),
                                        density=True)
            
            diluteNMax = 0.
            dilutePhi = 0.
            denseNMax = 0.
            densePhi = 0.
            for g in xrange(len(bins) - 1):
                if N[g] > diluteNMax and bins[g] < phi:
                    diluteNMax = N[g]
                    dilutePhi = bins[g]
                if N[g] > denseNMax and bins[g] > phi:
                    denseNMax = N[g]
                    densePhi = bins[g]
                    
            # If there is a large population in the dense phase, we know there is a dilute phase
            if denseNMax > 10.:
                diluteNMax = 0
                for g in xrange(len(bins) - 1):
                    if N[g] > diluteNMax and bins[g] < 0.4:
                        diluteNMax = N[g]
                        dilutePhi = bins[g]
            
            plt.axvline(dilutePhi, c='r')
            plt.axvline(densePhi, c='g')
                     
            plt.xlim(minDense[z], maxDense[z])
            plt.xlabel(r'Local Area Fraction $(\phi_{local})$')
            plt.ylabel(r'Population')
            plt.savefig('look_distance' + str(lookDist[z]) +\
                        '_pa' + str(peA) + '_pb' + str(peB) + '_xa' + str(parFrac) + '_phi' + str(intPhi) +\
                        '_fm' + str(pad) +'.jpg',
                        dpi=1000, bbox_inches='tight', pad_inches=0.05)
            plt.close()
            
            # Write the data to the textfile
            f = open(txtFile, 'a')
            f.write('{0:.4f}'.format(dilutePhi).center(10) + ' ')
            f.write('{0:.4f}'.format(densePhi).center(10))
            if z < (len(lookDist) - 1):
                f.write(' ')
            else:
                f.write('\n')
            f.close()
