'''
#                           This is an 80 character line                       #
What does this file do?
(Reads single argument, .gsd file name)

Track MSD of particles that stay within the dense phase.

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
    
    # Instantiate cluster computation
    f_box = box.Box(Lx=l_box, Ly=l_box, is2D=True)
    my_clust = cluster.Cluster()
    c_props = cluster.ClusterProperties()

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
        
        # Compute clusters for this timestep
        system = freud.AABBQuery(f_box, f_box.wrap(pos))
        # Compute neighbor list for only largest cluster
        my_clust.compute(system, neighbors={'r_max': 1.0})
        ids = my_clust.cluster_idx              # get id of each cluster
        c_props.compute(system, ids)            # find cluster properties
        clust_size = c_props.sizes              # find cluster sizes
        lcID = np.where(clust_size == np.amax(clust_size))
        
        if j == start:
            inClust = []
            for m in range(0, partNum):
                if ids[m] == lcID:
                    inClust.append(m)
                    
        # Loop over IDs that are still in the cluster
        for k in range(len(inClust)-1, -1, -1):
            # If the particle is no longer in the largest cluster
            if ids[inClust[k]] != lcID:
                # Remove it from the list
                del inClust[k]
    
    # Loop back through to get MSD of particles that remain in cluster
    msd = [ 0. for j in range(0, end) ]
    print(len(msd))
    # Get initial positions
    snap = t[0]
    pos = snap.particles.position               # position
    pos[:,-1] = 0.0
    
    initPos = []
    for j in range(0, partNum):
        initPos.append((pos[j][0], pos[j][1]))
        
    time = [0. for i in range(0, end)]
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
        time[j] = tst
        
        # Create frame pad for images
        pad = str(j).zfill(4)
        
        # Compute average MSD of all cluster-bound particles
        msd[j] += msd[j-1]
        for k in inClust:
            # (x_now - x_init)^2 + (y_now - y_init)^2
            msd[j] += (pos[k][0] - initPos[k][0])**2 + (pos[k][1] - initPos[k][1])**2

#
out = 'msd_' + infile
# Divide this by the number of particles
pltmsd = [ i / len(inClust) for i in msd ]
# Plot the MSD
plt.plot(time, pltmsd, lw=1.5)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Time $(\tau_{r})$')
plt.ylabel(r'MSD')
plt.xlim(time[1], time[-1])
plt.tight_layout()
plt.savefig(out + '_taur.png', dpi=1000)
plt.close()

btime = [ i * 3. for i in time ]
plt.plot(btime, pltmsd, lw=1.5)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Time $(\tau_{B})$')
plt.ylabel(r'MSD')
plt.xlim(btime[1], btime[-1])
plt.tight_layout()
plt.savefig(out + '_taub.png', dpi=1000)
plt.close()

# Write textfile
outTxt = 'msd_' + infile + '.txt'
g = open(outTxt, 'w') # write file headings
g.write('TauR'.center(30) + ' ' +\
        'MSD'.center(30) + '\n')
g.close()

# Append data to file
g = open(outTxt, 'a')
for j in range(start, end):
    g.write('{0:.3f}'.format(time[j]).center(30) + ' ')
    g.write('{0:.6f}'.format(pltmsd[j]).center(30) + '\n')
g.close()
