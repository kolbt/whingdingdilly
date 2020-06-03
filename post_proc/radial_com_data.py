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
out = outfile

# Get dumps to output
f = hoomd.open(name=infile, mode='rb')  # open gsd file with hoomd
dumps = int(f.__len__())                # get number of timesteps dumped
start = 0
#start = dumps - 1                       # gives first frame to read
end = dumps                             # gives last frame to read
#end = 20
start = dumps - 100

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
    
    # Get the largest x-position in the largest cluster
    f_box = box.Box(Lx=l_box, Ly=l_box, is2D=True)
    my_clust = cluster.Cluster()
    c_props = cluster.ClusterProperties()
    density = freud.density.LocalDensity(r_max=10., diameter=1.0)
    
    # You need one array of length nBins to hold the sum
    phi_sum = [ 0. for i in range(0, nBins) ]
    # You need one array of length nBins to hold the counts
    phi_num = [ 0. for i in range(0, nBins) ]
    # List to hold the average
    phi_avg = []
    # You should store the max distance of each bin as well
    r_bins = np.arange(sizeBin, sizeBin * nBins, sizeBin)
    
    # Accumulate data on figure to see if averaging is working
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))

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
        
        # Compute the center of mass
        system = freud.AABBQuery(f_box, f_box.wrap(pos))
        # Compute neighbor list for only largest cluster
        my_clust.compute(system, neighbors={'r_max': 1.0})
        ids = my_clust.cluster_idx              # get id of each cluster
        c_props.compute(system, ids)            # find cluster properties
        clust_size = c_props.sizes              # find cluster sizes
        com = c_props.centers[0]                # find cluster CoM
        lcID = np.where(clust_size == np.amax(clust_size))
        
        # Keep only positions of largest cluster (faster computations)
        lc_pos = []
        r_com = []
        lc_ids = []
        for k in range(0, partNum):
            if ids[k] == lcID:
                lc_pos.append(pos[k])
                lc_ids.append(k)
                # See if particle should be wrapped
                rx = np.abs(lc_pos[-1][0] - com[0])
                if rx >= h_box:
                    rx -= l_box
                    rx = np.abs(rx)
                ry = np.abs(lc_pos[-1][1] - com[1])
                if ry >= h_box:
                    ry -= l_box
                    ry = np.abs(ry)
                r_com.append(np.sqrt(rx**2 + ry**2))

        # Compute density around largest cluster points
        phi_locs = density.compute(system, query_points=lc_pos)
        phi_loc = phi_locs.density * np.pi * 0.25
        
        # Add/increment each particle in appropriate index
        count = 0
        for k in lc_ids:
            # Convert r to appropriate bin
            tmp_r = int(r_com[count] / sizeBin)
            phi_sum[tmp_r] += phi_loc[count]
            phi_num[tmp_r] += 1
            count += 1
        
#        # Accumulate data on plot
#        ax[0].scatter(r_com, phi_loc, s=0.25)
          
# Compute the average in each bin
for k in range(0, len(phi_sum)):
    if phi_num[k] > 0:
        phi_avg.append(phi_sum[k] / phi_num[k])
    else:
        phi_avg.append(0.)
        
# Write textfile
outTxt = 'CoM_' + out + '.txt'
g = open(outTxt, 'w') # write file headings
g.write('rCoM'.center(30) + ' ' +\
        'phiAvg'.center(30) + '\n')
g.close()
# Append data to file
g = open(outTxt, 'a')
for j in range(0, len(phi_avg)):
    g.write('{0:.5f}'.format(r_bins[j]).center(30) + ' ')
    g.write('{0:.6f}'.format(phi_avg[j]).center(30) + '\n')
g.close()

## Plot scatter of phi_loc vs r_com
#outDPI = 500
##plt.scatter(r_com, phi_loc, s=0.25)
#ax[1].scatter(r_bins, phi_avg, s=0.5, c='k')
#ax[0].set_xlim(0, 100)
#ax[1].set_xlim(0, 100)
#ax[0].set_ylim(0.5, 2.5)
#ax[1].set_ylim(0.5, 2.5)
#plt.tight_layout()
#plt.savefig("CoM_" + out + ".png", dpi=outDPI)
#plt.close()
