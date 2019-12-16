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

import freud
from freud import parallel
from freud import box
from freud import voronoi

import numpy as np
import math
import random
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.collections
from matplotlib import colors
from scipy.spatial import Voronoi, voronoi_plot_2d

cmaps = ['viridis', 'plasma', 'inferno', 'magma',
'cividis', 'Greys', 'Purples', 'Blues',
'Greens', 'Oranges', 'Reds', 'YlOrBr',
'YlOrRd', 'OrRd', 'PuRd', 'RdPu',
'BuPu', 'GnBu', 'PuBu', 'YlGnBu',
'PuBuGn', 'BuGn', 'YlGn', 'binary',
'gist_yarg', 'gist_gray', 'gray', 'bone',
'pink', 'spring', 'summer', 'autumn',
'winter', 'cool', 'Wistia', 'hot',
'afmhot', 'gist_heat', 'copper', 'PiYG',
'PRGn', 'BrBG', 'PuOr', 'RdGy',
'RdBu', 'RdYlBu', 'RdYlGn', 'Spectral',
'coolwarm', 'bwr', 'seismic', 'hsv',
'Pastel1', 'Pastel2', 'Paired', 'Accent',
'Dark2', 'Set1', 'Set2', 'Set3',
'tab10', 'tab20', 'tab20b', 'tab20c',
'flag', 'prism', 'ocean', 'gist_earth',
'terrain', 'gist_stern', 'gnuplot', 'gnuplot2',
'CMRmap', 'cubehelix', 'brg', 'gist_rainbow',
'rainbow', 'jet', 'nipy_spectral', 'gist_ncar']

cmaps_i_like = ['Blues', 'magma', 'ocean', 'viridis_r', 'winter', 'YlGnBu', 'plasma_r']

# Set the colormap (neighbors = [1, 2, 3, 4, 5, 6])
colorsList = ['#36413E', '#BEB2C8', '#757083', '#D7D6D6', '#70ABAF', '#66A182']
colorsList = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33']
colorsList = ['#FF595E', '#FFCA3A', '#8AC926', '#1982C4', '#6A4C93', '#E56399']
custom_cmap = colors.ListedColormap(colorsList)

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
    
def draw_voronoi(box, points, cells, colorArray, nlist=None, inMap='viridis'):
    
    # Set up the figure
    fig = plt.figure()
    fig.patch.set_facecolor('none')
    ax = fig.add_subplot(111)
    ax.patch.set_facecolor('none')
    ax.set_aspect('equal')
    ax.axis('off')
    # Draw Voronoi cells
    normColor = np.divide(colorArray, max(colorArray))
    patches = [plt.Polygon(cell[:, :2]) for cell in cells]
    patch_collection = matplotlib.collections.PatchCollection(patches, edgecolors='black', alpha=1.0, lw=0.15)

    # Set the color array
    patch_collection.set_array(np.ravel(nearNeigh))
    patch_collection.set_cmap(inMap)
    patch_collection.set_clim([min(nearNeigh), max(nearNeigh)])
    ax.add_collection(patch_collection)

    # Draw points
#    ax.scatter(points[:,0], points[:,1], c='k', s=0.5, edgecolors='none')
    ax.set_xlim((-box.Lx/2., box.Lx/2.))
    ax.set_ylim((-box.Ly/2., box.Ly/2.))
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    if not isinstance(inMap, str):
        inMap = 'custom_cmap'
    plt.savefig('voronoi_art_freud_' + inMap + '.png', dpi=1000, bbox_inches='tight', pad_inches=0)
    plt.close()

inFile = str(sys.argv[1])
snapList = [2]
r_cut = 2**(1./6.)
factor = r_cut
with hoomd.open(name=inFile, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    l_box = box_data[0]
    h_box = l_box / 2.
    f_box = box.Box(Lx = l_box, Ly = l_box, is2D = True)
    voro = voronoi.Voronoi(box=f_box, buff=(l_box / 2.0))
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
        
        voro.compute(box=f_box, positions=pos, buff=h_box)
#        draw_voronoi(f_box, pos, voro.polytopes, nearNeigh, inMap=custom_cmap)
        # This loops through all available colormaps
        for k in cmaps:
            draw_voronoi(f_box, pos, voro.polytopes, nearNeigh, inMap=k)
            draw_voronoi(f_box, pos, voro.polytopes, nearNeigh, inMap=k+'_r')
