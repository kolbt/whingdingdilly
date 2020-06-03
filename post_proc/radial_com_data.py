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
    
    # You need a list of length nBins to hold the sum
    phi_sum = [ 0. for i in range(0, nBins) ]
    p_sum = [ 0. for i in range(0, nBins) ]
    pswim_sum = [ 0. for i in range(0, nBins) ]
    pint_sum = [ 0. for i in range(0, nBins) ]
    # You need one array of length nBins to hold the counts
    num = [ 0. for i in range(0, nBins) ]
    # List to hold the average
    phi_avg = []
    p_avg = []
    pint_avg = []
    pswim_avg = []
    # You should store the max distance of each bin as well
    r_bins = np.arange(sizeBin, sizeBin * nBins, sizeBin)

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
        aligns = []
        for k in range(0, partNum):
            if ids[k] == lcID:
                lc_pos.append(pos[k])
                lc_ids.append(k)
                # See if particle should be wrapped
                rrx = lc_pos[-1][0] - com[0]
                rx = np.abs(rrx)
                if rx >= h_box:
                    rx -= l_box
                    rx = np.abs(rx)
                    # How should rrx be adjusted
                    if rrx < -h_box:
                        rrx += l_box
                    else:
                        rrx -= l_box
                rry = lc_pos[-1][1] - com[1]
                ry = np.abs(rry)
                if ry >= h_box:
                    print("Wrapping")
                    ry -= l_box
                    ry = np.abs(ry)
                    # How should rrx be adjusted
                    if rry < -h_box:
                        rry += l_box
                    else:
                        rry -= l_box
                # We want the vector from ref to CoM
                rrx = -rrx
                rry = -rry
                mag = np.sqrt(rx**2 + ry**2)
                r_com.append(mag)
                # Now let's get the x and y components of the body axis
                px = np.sin(ang[k])
                py = -np.cos(ang[k])
                # Now compute the dot product
                r_dot_p = (rrx * px) + (rry * py)
                # A value of 1 is perfectly aligned
                r_dot_p /= mag
                # We don't need to normalize by p, it's magnitude is 1
                aligns.append(r_dot_p)
        
        # Compute interparticle pressure
        pressure = [0. for i in range(0, len(lc_pos))]
        # Create the neighborlist of the system
        lc = freud.locality.AABBQuery(box=f_box, points=f_box.wrap(lc_pos))
        nlist = lc.query(lc_pos, dict(r_min=0.1, r_max=r_cut)).toNeighborList()
        
        # Loop through and compute pressure
        pairs = set()
        for (m, n) in nlist:
            # Never compute things twice
            if (m, n) in pairs or (n, m) in pairs:
                continue
            # So we know we've computed it
            pairs.add( (m, n) )
            # Loops through each j neighbor of reference particle i
            xx, yy, rr = distComps(lc_pos[m], lc_pos[n][0], lc_pos[n][1])
            fx, fy = computeFLJ(rr, xx, yy, eps)
            # Compute the x force times x distance
            sigx = fx * (xx)
            # Likewise for y
            sigy = fy * (yy)
            pressure[m] += ((sigx + sigy) / 2.)
            pressure[n] += ((sigx + sigy) / 2.)
                
#        # Let's take a look at the max alignment
#        print(max(aligns))
#        print(min(aligns))
#        # Sanity check is good values are [-1, 1]
#        inXY = np.delete(lc_pos, 2, 1)
#        # Let's do this as a patch collections instead (set diameter to 1.)
#        diams = [0.9 for i in range(len(lc_pos))]
#        outDPI = 150.
#        fig, ax = plt.subplots(1, 1)
#        coll = matplotlib.collections.EllipseCollection(diams, diams,
#                                                        np.zeros_like(diams),
#                                                        offsets=inXY, units='xy',
#                                                        cmap=plt.cm.viridis,
#                                                        transOffset=ax.transData)
##        coll.set_array(np.ravel(aligns))
##        coll.set_array(np.ravel(rrys))
##        coll.set_array(np.ravel(pys))
##        coll.set_array(np.ravel(pys))
##        coll.set_array(np.ravel(lc_angs))
#        coll.set_array(np.ravel(pressure))
##        coll.set_clim([-np.pi, np.pi])
##        coll.set_clim([-1., 1,])
#        ax.add_collection(coll)
#        cbar = plt.colorbar(coll, format='%.1f')
#        cbar.set_label(r'Pressure', labelpad=20, rotation=270)
#
#
#        # Add time as well
#        ax.text(0.95, 0.025, s=r'$\tau_{r}=$' + '{:0.1f}'.format(tst*3.),
#                horizontalalignment='right', verticalalignment='bottom',
#                transform=ax.transAxes,
#                fontsize=18,
#                bbox=dict(facecolor=(1,1,1,0.5), edgecolor=(0,0,0,1), boxstyle='round, pad=0.1'))
#        # Adjust limits and plotting parameters
#        ax.set_xlim(-h_box, h_box)
#        ax.set_ylim(-h_box, h_box)
#        ax.axes.set_xticks([])
#        ax.axes.set_yticks([])
#        ax.set_aspect('equal')
#        plt.tight_layout()
#        plt.show()

        # Compute density around largest cluster points
        phi_locs = density.compute(system, query_points=lc_pos)
        phi_loc = phi_locs.density * np.pi * 0.25
        
        # Add/increment each particle in appropriate index
        for k in range(0, len(lc_ids)):
            # Convert r to appropriate bin
            tmp_r = int(r_com[k] / sizeBin)
            p_sum[tmp_r] += aligns[k]
            phi_sum[tmp_r] += phi_loc[k]
            pint_sum[tmp_r] += pressure[k]
            num[tmp_r] += 1
          
## Compute the average in each bin
#for k in range(0, len(phi_sum)):
#    if num[k] > 0:
#        phi_avg.append(phi_sum[k] / num[k])
#        p_avg.append(p_sum[k] / num[k])
#        pint_avg.append(pint_sum[k] / num[k])
#    else:
#        phi_avg.append(0.)
#        p_avg.append(0.)
#        pint_avg.append(0.)

# Write textfile
outTxt = 'CoM_' + out + '.txt'
g = open(outTxt, 'w') # write file headings
g.write('rCoM'.center(25) + ' ' +\
        'NinBin'.center(25) + ' ' +\
        'phiLoc'.center(25) + ' ' +\
        'align'.center(25) + ' ' +\
        'pInt'.center(25) + '\n')
g.close()
# Append data to file
g = open(outTxt, 'a')
for j in range(0, len(phi_avg)):
    g.write('{0:.6f}'.format(r_bins[j]).center(25) + ' ')
    g.write('{0:.0f}'.format(num[j]).center(25) + ' ')
    g.write('{0:.6f}'.format(phi_sum[j]).center(25) + ' ')
    g.write('{0:.6f}'.format(p_sum[j]).center(25) + ' ')
    g.write('{0:.1f}'.format(pint_sum[j]).center(25) + '\n')
g.close()

## Plot scatter of phi_loc vs r_com
#outDPI = 500
#plt.plot(r_bins, pint_avg, lw=1.5, c='r', zorder=0)
#plt.scatter(r_bins, pint_avg, s=5, c='k', zorder=1)
#plt.xlim(0,)
#plt.tight_layout()
#plt.savefig("alignment_" + out + ".png", dpi=outDPI)
#plt.close()
