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
outFile = add + 'pressure_pa' + str(peA) +\
          '_pb' + str(peB) +\
          '_xa' + str(parFrac) +\
          '_phi' + str(intPhi) +\
          '_ep' + '{0:.3f}'.format(eps) +\
          '.txt'

g = open(outFile, 'w') # write file headings
g.write('Timestep'.center(10) + ' ' +\
        'gasArea'.center(20) + ' ' +\
        'gasTrace'.center(20) + ' ' +\
        'gasPress'.center(20) + ' ' +\
        'bulkArea'.center(20) + ' ' +\
        'bulkTrace'.center(20) + ' ' +\
        'bulkPress'.center(20) + ' ' +\
        'SurfaceTense'.center(20) + ' ' +\
        'Length'.center(20) + '\n')
g.close()

start = 0                   # first frame to process
dumps = int(f.__len__())    # get number of timesteps dumped
end = dumps                 # final frame to process
#start = end - 1
end = start + 1

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
        min_size = 5000
        
        # Get the positions of all particles in LC
        binParts = [[[] for b in range(NBins)] for a in range(NBins)]
        occParts = [[0 for b in range(NBins)] for a in range(NBins)]
        edgeBin = [[0 for b in range(NBins)] for a in range(NBins)]
        liqPos = []
        gasPos = []
        for k in range(0, len(ids)):
            # Convert position to be > 0 to place in list mesh
            tmp_posX = pos[k][0] + h_box
            tmp_posY = pos[k][1] + h_box
            x_ind = int(tmp_posX / sizeBin)
            y_ind = int(tmp_posY / sizeBin)
            # Append all particles to appropriate bin
            binParts[x_ind][y_ind].append(k)
            
            # Get sufficient cluster mesh as well
            if clust_size[ids[k]] >= min_size:
                liqPos.append(pos[k])
                occParts[x_ind][y_ind] = 1
            # Get a gas particle list as well
            elif clust_size[ids[k]] <= 100:
                gasPos.append(pos[k])
        
        # If sufficient neighbor bins are empty, we have an edge
        thresh = 1.5
        # Loop through x index of mesh
        for ix in range(0, len(occParts)):
        
            # If at right edge, wrap to left
            if (ix + 1) != NBins:
                lookx = [ix-1, ix, ix+1]
            else:
                lookx = [ix-1, ix, 0]
                
            # Loop through y index of mesh
            for iy in range(0, len(occParts[ix])):
            
                # Reset neighbor counter
                count = 0
                # If the bin is not occupied, skip it
                if occParts[ix][iy] == 0:
                    continue
                # If at top edge, wrap to bottom
                if (iy + 1) != NBins:
                    looky = [iy-1, iy, iy+1]
                else:
                    looky = [iy-1, iy, 0]
                
                # Loop through surrounding x-index
                for indx in lookx:
                    # Loop through surrounding y-index
                    for indy in looky:
                    
                        # If neighbor bin is NOT occupied
                        if occParts[indx][indy] == 0:
                            # If neighbor bin shares a vertex
                            if indx != ix and indy != iy:
                                count += 0.5
                            # If neighbor bin shares a side
                            else:
                                count += 1
                                
                # If sufficient neighbors are empty, we found an edge
                if count >= thresh:
                    edgeBin[indx][indy] = 1
        
        blurBin = [[0 for b in range(NBins)] for a in range(NBins)]
        Nedges = 0
        # Blur the edge bins a bit
        for ix in range(0, len(occParts)):
            for iy in range(0, len(occParts)):
                Nedges += edgeBin[ix][iy]
                if edgeBin[ix][iy] == 1:
                    # If at right edge, wrap to left
                    if (ix + 1) != NBins:
                        lookx = [ix-1, ix, ix+1]
                    else:
                        lookx = [ix-1, ix, 0]
                    # If at top edge, wrap to bottom
                    if (iy + 1) != NBins:
                        looky = [iy-1, iy, iy+1]
                    else:
                        looky = [iy-1, iy, 0]
                    # Loop through surrounding x-index
                    for indx in lookx:
                        # Loop through surrounding y-index
                        for indy in looky:
                            # Make all surrounding bins 'edge' bins
                            blurBin[indx][indy] = 1
        
        # Now let's compute the pressure
        bulkSigXX = 0
        bulkSigYY = 0
        gasSigXX = 0
        gasSigYY = 0
        print("Computing the pressure")
        for k in range(0, len(ids)):
            # Convert position to be > 0 to place in list mesh
            tmp_posX = pos[k][0] + h_box
            tmp_posY = pos[k][1] + h_box
            x_ind = int(tmp_posX / sizeBin)
            y_ind = int(tmp_posY / sizeBin)
            # Only compute pressure for non-edge bins
            if blurBin[x_ind][y_ind] != 1:
                # If at right edge, wrap to left
                if (x_ind + 1) != NBins:
                    lookx = [x_ind-1, x_ind, x_ind+1]
                else:
                    lookx = [x_ind-1, x_ind, 0]
                # If at top edge, wrap to bottom
                if (y_ind + 1) != NBins:
                    looky = [y_ind-1, y_ind, y_ind+1]
                else:
                    looky = [y_ind-1, y_ind, 0]
                # Reference particle position
                refx = pos[k][0]
                refy = pos[k][1]
                # Loop through surrounding x-index
                for indx in lookx:
                    # Loop through surrounding y-index
                    for indy in looky:
                        # Loop through all particles in that bin
                        for comp in binParts[indx][indy]:
                            # Compute the distance
                            dist = computeDist(refx, refy, pos[comp][0], pos[comp][1])
                            # If potential is on ...
                            if 0.1 < dist <= r_cut:
                                # Compute the x and y components of force
                                fx, fy = computeFLJ(dist, refx, refy, pos[comp][0], pos[comp][1], eps)
                                # This will go into the bulk pressure
                                if clust_size[ids[k]] >= min_size:
                                    bulkSigXX += (fx * (pos[comp][0] - refx))
                                    bulkSigYY += (fy * (pos[comp][1] - refy))
                                # This goes into the gas pressure
                                else:
                                    gasSigXX += (fx * (pos[comp][0] - refx))
                                    gasSigYY += (fy * (pos[comp][1] - refy))
        
        # Now let's get the area of each phase (by summing bin areas)
        gasBins = 0
        bulkBins = 0
        testIDs = [[0 for b in range(NBins)] for a in range(NBins)]
        for ix in range(0, len(occParts)):
            for iy in range(0, len(occParts)):
                # Is the bin an edge?
                if blurBin[ix][iy] == 1:
                    testIDs[ix][iy] = 0
                    continue
                # Does the bin belong to the dense phase?
                if len(binParts[ix][iy]) != 0:
                    if clust_size[ids[binParts[ix][iy][0]]] > min_size:
                        bulkBins += 1
                        testIDs[ix][iy] = 1
                        continue
                gasBins += 1
                testIDs[ix][iy] = 2
                 
        # The edge length of sufficiently large clusters
        lEdge = Nedges * sizeBin
        # Divide by two because each pair is counted twice
        bulkTrace = (bulkSigXX + bulkSigYY)/2.
        gasTrace = (gasSigXX + gasSigYY)/2.
        # Area of a bin
        binArea = sizeBin * sizeBin
        # Area of each phase
        gasArea = binArea * gasBins
        bulkArea = binArea * bulkBins
        # Pressure of each phase
        gasPress = gasTrace / gasArea
        bulkPress = bulkTrace / bulkArea
        # Surface tension
        surfaceTense = (bulkPress - gasPress) * lEdge
        
        print("Number of gas bins: ", gasBins)
        print("Gas phase area: ", gasArea)
        print("Number of bulk bins: ", bulkBins)
        print("Bulk phase area: ", bulkArea)
        print("Trace of gas stress tensor: ", gasTrace)
        print("Trace of bulk stress tensor: ", bulkTrace)
        print("Gas phase pressure: ", gasPress)
        print("Bulk phase pressure: ", bulkPress)
        print("Surface tension: ", surfaceTense)
        
#        # A sanity check on a perfect hcp circle
#        print(Nedges)
#        print(Nedges * sizeBin)
#        x = list(list(zip(*liqPos))[0])
#        y = list(list(zip(*liqPos))[1])
#        diam = max(x) - min(x)
#        circ = diam * np.pi
#        print(circ)
#        print(Nedges * sizeBin / circ)

#        # Let's plot imshow to make sure we're good thus far
#        fig, ax = plt.subplots()
#        ax.imshow(testIDs, extent=[0, l_box, 0, l_box], aspect='auto', origin='lower')
#        ax.set_aspect('equal')
#        plt.show()
                
        # Write this to a textfile with the timestep
        g = open(outFile, 'a')
        g.write('{0:.3f}'.format(tst).center(10) + ' ')
        g.write('{0:.3f}'.format(gasArea).center(20) + ' ')
        g.write('{0:.3f}'.format(gasTrace).center(20) + ' ')
        g.write('{0:.3f}'.format(gasPress).center(20) + ' ')
        g.write('{0:.3f}'.format(bulkArea).center(20) + ' ')
        g.write('{0:.3f}'.format(bulkTrace).center(20) + ' ')
        g.write('{0:.3f}'.format(bulkPress).center(20) + ' ')
        g.write('{0:.3f}'.format(surfaceTense).center(20) + ' ')
        g.write('{0:.1f}'.format(lEdge).center(20) + '\n')
        g.close()
