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

# Run locally
sys.path.append('/Users/kolbt/Desktop/compiled/hoomd-blue/build')
sys.path.append('/Users/kolbt/Desktop/compiled/gsd/build')
# Run on the cpu
sys.path.append('/nas/longleaf/home/kolbt/programs/cpu-hoomd/hoomd-blue/build')
# Run on the gpu
sys.path.append('/nas/longleaf/home/kolbt/programs/hoomd_2.2.1/hoomd-blue/build')
sys.path.append('/nas/longleaf/home/kolbt/programs/gsd/build')

import gsd
from gsd import hoomd
from gsd import pygsd

import freud
from freud import parallel
from freud import box
from freud import density
from freud import cluster

import math
import numpy as np
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.collections
from matplotlib import collections  as mc
from matplotlib import lines

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
    
def getR(x, y):
    return np.sqrt(x**2 + y**2)
    
# Grab files
slowCol = '#d8b365'
fastCol = '#5ab4ac'

# Command line arguments
infile = str(sys.argv[1])                               # gsd file
peA = float(sys.argv[2])
peB = float(sys.argv[3])
parFrac = float(sys.argv[4])
eps = float(sys.argv[5])
try:
    phi = float(sys.argv[6])
    intPhi = int(phi)
    phi /= 100.
except:
    phi = 0.6
    intPhi = 60
try:
    dtau = float(sys.argv[7])
except:
    dtau = 0.000001
    
out = "final_pe" + "{:.0f}".format(peA) +\
      "_phi" + "{:.0f}".format(intPhi) +\
      "_eps" + "{:.5f}".format(eps) +\
      "_fm"
    
# Create outfile name from infile name
file_name = os.path.basename(infile)
outfile, file_extension = os.path.splitext(file_name)   # get base name
out = outfile + "_frame_"

# Get dumps to output
f = hoomd.open(name=infile, mode='rb')  # open gsd file with hoomd
dumps = int(f.__len__())                # get number of timesteps dumped
start = 0
#start = dumps - 1                       # gives first frame to read
end = dumps                             # gives last frame to read
#end = 20

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

# Round up size of bins to account for floating point inaccuracy
def roundUp(n, decimals=0):
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier

# Compute mesh
r_cut = 2**(1./6.)

# List to hold initial distances
distances = []

# Access file frames
with hoomd.open(name=infile, mode='rb') as t:

    # Take first snap for box
    snap = t[0]
    first_tstep = snap.configuration.step
    box_data = snap.configuration.box
    # Get box dimensions
    l_box = box_data[0]
    h_box = l_box / 2.
    a_box = l_box * l_box
    nBins = (getNBins(l_box, r_cut))
    sizeBin = roundUp((l_box / nBins), 6)
    partNum = len(snap.particles.typeid)
    pos = snap.particles.position
    
    # Get the largest x-position in the largest cluster
    f_box = box.Box(Lx=l_box, Ly=l_box, is2D=True)
    my_clust = cluster.Cluster()
    c_props = cluster.ClusterProperties()
    # Compute clusters for this timestep
    system = freud.AABBQuery(f_box, f_box.wrap(pos))
    # Compute neighbor list for only largest cluster
    my_clust.compute(system, neighbors={'r_max': 1.0})
    ids = my_clust.cluster_idx              # get id of each cluster
    c_props.compute(system, ids)            # find cluster properties
    clust_size = c_props.sizes              # find cluster sizes
    lcID = np.where(clust_size == np.amax(clust_size))
    
    maxx = 0.0
    # Loop through all particles
    for j in range(0, partNum):
        # Check if in cluster
        if ids[j] == lcID:
            if pos[j][0] > maxx:
                maxx = pos[j][0]
    
    # Compute the colormap
    for j in range(0, partNum):
        r = getR(pos[j][0], pos[j][1])
#        r /= float(h_box/2.)
#        if r > maxx:
#            r = 0.0
        r /= maxx
        distances.append(r)

    # Loop through snapshots
    for j in range(start, end):
    
        # Get the current snapshot
        snap = t[j]
        # Easier accessors
        pos = snap.particles.position               # position
        pos[:,-1] = 0.0
        xy = np.delete(pos, 2, 1)
        typ = snap.particles.typeid                 # type
        tst = snap.configuration.step               # timestep
        tst -= first_tstep                          # normalize by first timestep
        tst *= dtau                                 # convert to Brownian time
        
        # Create frame pad for images
        pad = str(j).zfill(4)
        
        # Let's do this as a patch collections instead (set diameter to 1.)
        diams = [0.9 for i in range(partNum)]
        outDPI = 150.
        fig, ax = plt.subplots(1, 1)
        coll = matplotlib.collections.EllipseCollection(diams, diams,
                                                        np.zeros_like(diams),
                                                        offsets=xy, units='xy',
                                                        cmap=plt.cm.viridis_r,
                                                        transOffset=ax.transData)
        coll.set_array(np.ravel(distances))
        minCol = 0.0
        coll.set_clim([minCol, 1.0])
        ax.add_collection(coll)
        cbar = plt.colorbar(coll, format='%.1f')
        cbar.set_label(r'Initial radial distance $(r)$', labelpad=20, rotation=270)

        
        # Add time as well
        ax.text(0.95, 0.025, s=r'$\tau_{r}=$' + '{:0.1f}'.format(tst*3.),
                horizontalalignment='right', verticalalignment='bottom',
                transform=ax.transAxes,
                fontsize=18,
                bbox=dict(facecolor=(1,1,1,0.5), edgecolor=(0,0,0,1), boxstyle='round, pad=0.1'))
        # Adjust limits and plotting parameters
        ax.set_xlim(-h_box, h_box)
        ax.set_ylim(-h_box, h_box)
        ax.axes.set_xticks([])
        ax.axes.set_yticks([])
        ax.set_aspect('equal')
        plt.tight_layout()
        plt.savefig("radial_" + out + pad + ".png", dpi=outDPI)
        plt.close()
