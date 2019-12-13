'''
#                           This is an 80 character line                       #

Make something artisitc!

Pro tip: it can be a random distribution of data, no need to make it a
         gsd snapshot.

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

import numpy as np
import math
import random
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.collections
from matplotlib import colors
from scipy.spatial import Voronoi, voronoi_plot_2d

def getDistance(point1, point2x, point2y):
    """Find the distance between two points"""
    distance = np.sqrt((point2x - point1[0])**2 + (point2y - point1[1])**2)
    return distance

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
        for j in range(bottom, top):
            indX = i
            indY = j
            if i >= maxInds:
                indX -= maxInds
            if j >= maxInds:
                indY -= maxInds
            binsList[0].append(indX)
            binsList[1].append(indY)
    return binsList

inFile = str(sys.argv[1])
snapList = [100]
r_cut = 2**(1./6.)
factor = r_cut
with hoomd.open(name=inFile, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    l_box = box_data[0]
    h_box = l_box / 2.
    # Mesh stuff
    NBins = getNBins(l_box, r_cut)
    sizeBin = roundUp((l_box / NBins), 6)
    print("Number of bins: {}").format(NBins)
    print("Size of bins: {}").format(sizeBin)
    lookNBins = int( 1.5 / sizeBin ) + 1
    
    # Loop through frames to plot
    for i in snapList:
        snap = t[i]
        pos = snap.particles.position
        partNum = len(pos)
        xy = np.delete(pos, 2, 1)
        typ = snap.particles.typeid
        
        # Bins for particles (by position)
        binParts = [[[] for b in range(NBins)] for a in range(NBins)]
        # Defaults to diameter of 1
        effSigma = [1.] * partNum
        # Number of nearest neighbors
        nearNeigh = [0] * partNum
        
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
            lookBins = findMinBins(1.5, x_ind, y_ind, NBins, sizeBin)
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
        
        # Normalize color map for nearest neighbor
        norm = mpl.colors.Normalize(vmin=1, vmax=6, clip=True)
        mapper = cm.ScalarMappable(norm=norm, cmap='viridis')
        # Make voronoi diagram
        vor = Voronoi(xy)
        voronoi_plot_2d(vor, show_points=False, show_vertices=False,
                        line_colors='k', line_width=0.2)
        for r in range(len(vor.point_region)):
            region = vor.regions[vor.point_region[r]]
            if not -1 in region:
                polygon = [vor.vertices[s] for s in region]
                plt.fill(*zip(*polygon), edgecolor='none', linewidth=0,
                         color=mapper.to_rgba(nearNeigh[r]))
        plt.xticks([])
        plt.yticks([])
        plt.xlim(-h_box, h_box)
        plt.ylim(-h_box, h_box)
        ax = plt.gca()
        ax.set_aspect('equal')
        plt.savefig('voronoi_scipy_art.png', dpi=1000, bbox_inches='tight', pad_inches=0.05)
        plt.close
