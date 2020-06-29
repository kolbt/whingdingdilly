'''
#                           This is an 80 character line                       #
What does this file do?
(Reads single argument, .gsd file name)

1.) Get center of mass

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
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.collections
from matplotlib import collections  as mc
from matplotlib import lines

import matplotlib.gridspec as gridspec
    
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
    
# Create outfile name from infile name
file_name = os.path.basename(infile)
outfile, file_extension = os.path.splitext(file_name)   # get base name
out = outfile

# Get dumps to output
f = hoomd.open(name=infile, mode='rb')  # open gsd file with hoomd
dumps = int(f.__len__())                # get number of timesteps dumped
start = 0
end = dumps                             # gives last frame to read
start = dumps - 1

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
    
def quatToAngle(quat):
    "Take vector, output angle between [-pi, pi]"
    x = quat[1]
    y = quat[2]
    rad = math.atan2(y, x)
    return rad
    
def distComps(point1, point2x, point2y):
    '''Given points output x, y and r'''
    dx = point2x - point1[0]
    dy = point2y - point1[1]
    r = np.sqrt((dx**2) + (dy**2))
    return dx, dy, r
    
def computeFLJ(r, dx, dy, eps):
    sig = 1.
    f = (24. * eps / r) * ( (2*((sig/r)**12)) - ((sig/r)**6) )
    fx = f * (dx) / r
    fy = f * (dy) / r
    return fx, fy
    
# Get lattice spacing for particle size
def ljForce(r, eps, sigma=1.):
    div = (sigma/r)
    dU = (24. * eps / r) * ((2*(div**12)) - (div)**6)
    return dU
    
def avgCollisionForce(pe, power=1.):
    '''Computed from the integral of possible angles'''
    magnitude = np.sqrt(28.)
    return (magnitude * (pe)) / (np.pi)
    
def conForRClust(peA, peB, xA, eps):
    if xA > 1.:
        xA /= 100.
    penet = (peA * xA) + (peB * (1.-xA))
    out = []
    r = 1.112
    skip = [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001]
    for j in skip:
        while ljForce(r, eps) < avgCollisionForce(penet):
            r -= j
        r += j
    out = r
    return out

# Compute mesh
r_cut = 2**(1./6.)

lat = conForRClust(peA, peB, parFrac, eps)
myCols = plt.cm.plasma

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
    density = freud.density.LocalDensity(r_max=10., diameter=1.0)

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
        ori = snap.particles.orientation            # orientation
        ang = np.array(list(map(quatToAngle, ori))) # convert to [-pi, pi]
        
#        # Compute the center of mass
#        system = freud.AABBQuery(f_box, f_box.wrap(pos))
#        # Compute neighbor list for only largest cluster
#        my_clust.compute(system, neighbors={'r_max': r_cut})
#        ids = my_clust.cluster_idx              # get id of each cluster
#        c_props.compute(system, ids)            # find cluster properties
#        clust_size = c_props.sizes              # find cluster sizes
#        maxSize = np.amax(clust_size)
#        lcID = np.where(clust_size == maxSize)
        
        posA = []
        posB = []
        acount = 0
        bcount = 0
        for k in range(0, partNum):
            if typ[k] == 0:
                posA.append(pos[k])
                acount += 1
            else:
                posB.append(pos[k])
                bcount += 1
        print(acount)
        print(bcount)
                
        # Get diameter from prediction
        parts = [lat for i in range(0, len(posA))]
        system = freud.AABBQuery(f_box, f_box.wrap(posA))
        
        # Compute density around largest cluster points
        phi_locs = density.compute(system, query_points=posA)
        phi_locA = phi_locs.density * np.pi * 0.25
        
        # Get diameter from prediction
        parts = [lat for i in range(0, len(posA))]
        system = freud.AABBQuery(f_box, f_box.wrap(posB))
        
        # Compute density around largest cluster points
        phi_locs = density.compute(system, query_points=posB)
        phi_locB = phi_locs.density * np.pi * 0.25

    outDPI = 500.
#    fig, ax = plt.subplots(1, 2)
    fig = plt.figure(figsize=(11.5, 5))
    pgs = gridspec.GridSpec(2, 1, figure=fig, hspace=0.05)
    ax = []
    ax.append(fig.add_subplot(pgs[0]))
    ax.append(fig.add_subplot(pgs[1]))
    # Plot A particles
    xy = np.delete(posA, 2, 1)
    coll = matplotlib.collections.EllipseCollection(parts, parts,
                                                    np.zeros_like(parts),
                                                    offsets=xy, units='xy',
                                                    cmap=myCols,
                                                    transOffset=ax[0].transData)
    coll.set_array(np.ravel(phi_locA))
#    minCol = min(phi_loc)
#    maxCol = max(phi_loc)
    coll.set_clim([0.6, 1.0])
    ax[0].add_collection(coll)
#    cbar = plt.colorbar(coll, format='%.2f')
#    cbar.set_label(r'Local area fraction $(\phi_{loc})$', labelpad=20, rotation=270)

    # Limits and ticks
    viewBuff = 0.0
    ax[0].set_xlim(-h_box - viewBuff, h_box + viewBuff)
    ax[0].set_ylim(-h_box - viewBuff, h_box + viewBuff)
    ax[0].tick_params(axis='both', which='both',
                   bottom=False, top=False, left=False, right=False,
                   labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    pad = str(j).zfill(4)
    ax[0].set_aspect('equal')
    
    # Now plot B particles
    xy = np.delete(posB, 2, 1)
    coll = matplotlib.collections.EllipseCollection(parts, parts,
                                                    np.zeros_like(parts),
                                                    offsets=xy, units='xy',
                                                    cmap=myCols,
                                                    transOffset=ax[1].transData)
    coll.set_array(np.ravel(phi_locB))
#    minCol = min(phi_loc)
#    maxCol = max(phi_loc)
    coll.set_clim([0.6, 1.0])
    ax[1].add_collection(coll)
    cbar = fig.colorbar(coll, ax=ax, format='%.2f', fraction=0.05, pad=0.01)
#    cbar = plt.colorbar(coll, format='%.2f', ax=ax, fraction=0.05)
    cbar.set_label(r'Local area fraction $(\phi_{loc})$', labelpad=20, rotation=270)

    # Limits and ticks
    viewBuff = 0.0
    ax[1].set_xlim(-h_box - viewBuff, h_box + viewBuff)
    ax[1].set_ylim(-h_box - viewBuff, h_box + viewBuff)
    ax[1].tick_params(axis='both', which='both',
                   bottom=False, top=False, left=False, right=False,
                   labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    pad = str(j).zfill(4)
    ax[1].set_aspect('equal')
    
    # Plot text
    ax[0].text(0.01, 0.925, r'Slow-slow',  transform=ax[0].transAxes)
    ax[1].text(0.01, 0.925, r'Fast-fast',  transform=ax[1].transAxes)
    
    plt.savefig('typed_' + out + '_' + pad + '.png', bbox_inches='tight', pad_inches=0.02, dpi=outDPI)
    plt.close()
