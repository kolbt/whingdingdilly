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
    nBins = initGuess
    # This loop only exits on function return
    while True:
        if length / nBins > minSz:
            return nBins
        else:
            nBins -= 1

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

# Outfile to write data to
outFile = add + 'edge_pa' + str(peA) +\
          '_pb' + str(peB) +\
          '_xa' + str(parFrac) +\
          '_phi' + str(intPhi) +\
          '_ep' + '{0:.3f}'.format(eps) +\
          '.txt'
          
g = open(outFile, 'w') # write file headings
g.write('Timestep'.center(10) + ' ' + 'Length'.center(10) + '\n')
g.close()

start = 0                   # first frame to process
dumps = int(f.__len__())    # get number of timesteps dumped
end = dumps                 # final frame to process
start = end - 1

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
    # Set up cluster computation using box
    f_box = box.Box(Lx=l_box, Ly=l_box, is2D=True)
    my_clust = cluster.Cluster()
    c_props = cluster.ClusterProperties()
    # Compute each mesh
    NBins = getNBins(l_box, r_cut)
    sizeBin = roundUp((l_box / NBins), 6)
    
    # Loop through each timestep
    for j in range(start, end):
        snap = t[j]
        # Easier accessors
        pos = snap.particles.position               # position
        pos[:,-1] = 0.0
        xy = np.delete(pos, 2, 1)
        typ = snap.particles.typeid                 # type
        tst = snap.configuration.step               # timestep
        tst -= first_tstep                          # normalize by first timestep
        tst *= tauPerDT                             # convert to Brownian time
        
        # Compute clusters for this timestep
        system = freud.AABBQuery(f_box, f_box.wrap(pos))
        
        # Compute neighbor list for only largest cluster
        my_clust.compute(system, neighbors={'r_max': 1.0})
        ids = my_clust.cluster_idx              # get id of each cluster
        c_props.compute(system, ids)            # find cluster properties
        clust_size = c_props.sizes              # find cluster sizes

        # We can also grab all clusters over a set size
        suffIDs = []
        for k in range(0, len(clust_size)):
            # Cluster must be >= 5% of all particles
            if clust_size[k] >= (partNum * 0.05):
                # Grab sufficiently large cluster IDs
                suffIDs.append(k)
        
#        # Try and grab edge of largest cluster
#        lClust = max(clust_size)
#        # Get the id of the largest cluster
#        for k in range(0, len(clust_size)):
#            if clust_size[k] == lClust:
#                lcID = ids[k]

        # Get the positions of all particles in LC
        binParts = [[[] for b in range(NBins)] for a in range(NBins)]
        occParts = [[0 for b in range(NBins)] for a in range(NBins)]
        edgeBin = [[0 for b in range(NBins)] for a in range(NBins)]
        lcPos = []
        for k in range(0, len(ids)):
            if ids[k] in suffIDs:
#            if ids[k] == lcID:
                lcPos.append(pos[k])
                # Convert position to be > 0 to place in list mesh
                tmp_posX = pos[k][0] + h_box
                tmp_posY = pos[k][1] + h_box
                x_ind = int(tmp_posX / sizeBin)
                y_ind = int(tmp_posY / sizeBin)
                # Append particle id to appropriate bin
                binParts[x_ind][y_ind].append(k)
                occParts[x_ind][y_ind] = 1
        
        # Loop through indices, figure out which are empty
        thresh = 1.5
        for ix in range(0, len(occParts)):
            lookx = [ix-1, ix, ix+1]
            for iy in range(0, len(occParts[ix])):
                count = 0
                if occParts[ix][iy] == 0:
                    continue
                looky = [iy-1, iy, iy+1]
                for indx in lookx:
                    for indy in looky:
                        if occParts[indx][indy] == 0:
                            if indx != ix and indy != iy:
                                count += 0.5
                            else:
                                count += 1
                if count >= thresh:
                    edgeBin[indx][indy] = 1
        
        # Sum the resultant mesh
        Nedges = 0
        for ix in range(0, len(occParts)):
            for iy in range(0, len(occParts[ix])):
                Nedges += edgeBin[ix][iy]
        
        # The edge length of sufficiently large clusters
        lEdge = Nedges * sizeBin
        
        # Write this to a textfile with the timestep
        g = open(outFile, 'a')
        g.write('{0:.3f}'.format(tst).center(10) + ' ')
        g.write('{0:.1f}'.format(lEdge).center(10) + '\n')
        g.close()
        
#        # A sanity check on a perfect hcp circle
#        print(Nedges)
#        print(Nedges * sizeBin)
#        x = list(list(zip(*lcPos))[0])
#        y = list(list(zip(*lcPos))[1])
#        diam = max(x) - min(x)
#        circ = diam * np.pi
#        print(circ)
#        print(Nedges * sizeBin / circ)
#
#        # Let's plot imshow to make sure we're good thus far
#        fig, ax = plt.subplots()
#        ax.imshow(edgeBin, extent=[0, l_box, 0, l_box], aspect='auto', origin='lower')
#        ax.set_aspect('equal')
#        plt.show()
                
        
