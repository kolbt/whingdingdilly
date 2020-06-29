'''
#                           This is an 80 character line                       #
Compute the length of the cluster edge:
-Use Freud to find the complete system neighborlist
-Grab the largest cluster
-Mesh the system
-Compute which bins have largest cluster particles
-If adjacent bins are empty, the reference bin is an edge
-Multiply by bin size to get length
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

def computeDist(x1, y1, x2, y2):
    '''Compute distance between two points'''
    return np.sqrt( ((x2-x1)**2) + ((y2 - y1)**2) )
    
def computeFLJ(r, x1, y1, x2, y2, eps):
    sig = 1.
    f = (24. * eps / r) * ( (2*((sig/r)**12)) - ((sig/r)**6) )
    fx = f * (x2 - x1) / r
    fy = f * (y2 - y1) / r
    return fx, fy

def computeTauPerTstep(epsilon, mindt=0.000001):
    '''Read in epsilon, output tauBrownian per timestep'''
#    if epsilon != 1.:
#        mindt=0.00001
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
    NBins = initGuess
    # This loop only exits on function return
    while True:
        if length / NBins > minSz:
            return NBins
        else:
            NBins -= 1

def findBins(lookN, currentInd, maxInds):
    '''Get the surrounding bin indices'''
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
    
def quatToAngle(quat):
    "Take vector, output angle between [-pi, pi]"
    x = quat[1]
    y = quat[2]
    rad = math.atan2(y, x)
    return rad
    
def getMagnitude(vecF):
    "Take force vector, output magnitude"
    x = vecF[0]
    y = vecF[1]
    magF = np.sqrt((x**2)+(y**2))
    return magF

# Get infile and open
inFile = str(sys.argv[1])
if inFile[0:7] == "cluster":
    add = 'cluster_'
else:
    add = ''
    
f = hoomd.open(name=inFile, mode='rb')
# Inside and outside activity from command line
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
    
#start = end - 1
start = 260
end = 350
peA = 150.
peB = 500.

box_data = np.zeros((1), dtype=np.ndarray)  # box dimension holder
r_cut = 2**(1./6.)                          # potential cutoff
tauPerDT = computeTauPerTstep(epsilon=eps)  # brownian time per timestep

with hoomd.open(name=inFile, mode='rb') as t:
    snap = t[0]
    first_tstep = snap.configuration.step
    box_data = snap.configuration.box
    l_box = box_data[0]
    h_box = l_box / 2.
    typ = snap.particles.typeid
    partNum = len(typ)
#    # Set up cluster computation using box
#    f_box = box.Box(Lx=l_box, Ly=l_box, is2D=True)
#    my_clust = cluster.Cluster()
#    c_props = cluster.ClusterProperties()
    # Compute each mesh
    NBins = getNBins(l_box, r_cut) / 4
    sizeBin = roundUp((l_box / NBins), 6)
    
    # Loop through each timestep
    for j in range(start, end):
        snap = t[j]
        # Easier accessors
        pos = snap.particles.position               # position
        pos[:,-1] = 0.0
        pos[:, 0] -= 230.
        pos[:, 1] -= 130.
        for k in range(0, len(pos)):
            # Wrap x
            if pos[k][0] < -h_box:
                pos[k][0] += l_box
            if pos[k][0] > h_box:
                pos[k][0] -= l_box
            # Wrap y
            if pos[k][1] < -h_box:
                pos[k][1] += l_box
            if pos[k][1] > h_box:
                pos[k][1] -= l_box
                
        xy = np.delete(pos, 2, 1)
        typ = snap.particles.typeid                 # type
        tst = snap.configuration.step               # timestep
        tst -= first_tstep                          # normalize by first timestep
        tst *= dtau                                 # convert to Brownian time
        ori = snap.particles.orientation            # orientation
        ang = np.array(list(map(quatToAngle, ori))) # convert to [-pi, pi]
            
#        # Compute clusters for this timestep
#        system = freud.AABBQuery(f_box, f_box.wrap(pos))
#
#        # Compute neighbor list for only largest cluster
#        my_clust.compute(system, neighbors={'r_max': 1.0})
#        ids = my_clust.cluster_idx              # get id of each cluster
#        c_props.compute(system, ids)            # find cluster properties
#        clust_size = c_props.sizes              # find cluster sizes
#
#        # We can also grab all clusters over a set size
#        min_size = 5000
        
        # Get the positions of all particles in LC
        binParts = [[[] for b in range(NBins)] for a in range(NBins)]
        occParts = [[0 for b in range(NBins)] for a in range(NBins)]
        edgeBin = [[0 for b in range(NBins)] for a in range(NBins)]
        mesh = np.zeros((NBins, NBins, 2), dtype=np.float32)
        # Array to hold the magnitude of the binned force
        binnedF = np.zeros((NBins, NBins), dtype=np.float32)
        for k in range(0, partNum):
            # Convert position to be > 0 to place in list mesh
            tmp_posX = pos[k][0] + h_box
            tmp_posY = pos[k][1] + h_box
            x_ind = int(tmp_posX / sizeBin)
            y_ind = int(tmp_posY / sizeBin)
            # Append all particles to appropriate bin
            px = np.sin(ang[k])
            py = -np.cos(ang[k])
            binParts[x_ind][y_ind].append(k)
            if typ[k] == 0:
                px *= peA
                py *= peA
            else:
                px *= peB
                py *= peB
            mesh[x_ind][y_ind][0] += px
            mesh[x_ind][y_ind][1] += py
        # Take magnitude of each bin (removes directional component)
        for k in range(0, NBins):
            for l in range(0, NBins):
                binnedF[k][l] += getMagnitude(mesh[k][l])

        # Plot binned data using imshow
        pad = str(j).zfill(4)
        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(111)
#        plt.imshow(binnedF, origin='lower')
        im = plt.imshow(binnedF.T/2.,
                   extent=(0,NBins,0,NBins),
                   origin='lower',
                   vmin=0, vmax=5200)
        plt.xticks(())
        plt.yticks(())
        
        ax.text(0.95, 0.025, s=r'$\tau_{r}=$' + '{:0.1f}'.format(tst*3.),
                horizontalalignment='right', verticalalignment='bottom',
                transform=ax.transAxes,
                fontsize=18,
                bbox=dict(facecolor=(1,1,1,0.5), edgecolor=(0,0,0,1), boxstyle='round, pad=0.1'))
        
#        cb = fig.colorbar(im, ax=ax, fraction=0.045, orientation='vertical')
#        cb.set_label(r'$||\mathbf{F}_{act}||$', rotation=270, labelpad=15.)
        ax.set_aspect('equal')
        plt.savefig('f_act_bin_fm_' + pad + '.png', bbox_inches='tight', pad_inches=0.02, dpi=500)
        plt.close()
