'''
#                           This is an 80 character line                       #
Compute the per-particle active pressure
 (with both swim and interparticle contributions)
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
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.collections
from matplotlib import colors
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
    
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
            
# Get lattice spacing for particle size
def ljForce(r, eps, sigma=1.):
    div = (sigma/r)
    dU = (24. * eps / r) * ((2*(div**12)) - (div)**6)
    return dU
    
def avgCollisionForce(pe, power=1.):
    '''Computed from the integral of possible angles'''
    magnitude = np.sqrt(28.)
    return (magnitude * (pe)) / (np.pi)
    
def conForRClust(pe, eps):
    out = []
    r = 1.112
    skip = [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001]
    for j in skip:
        while ljForce(r, eps) < avgCollisionForce(pe):
            r -= j
        r += j
    out = r
    return out

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
    
lat = conForRClust(peA, eps)
    
# Get number of particles for textfile name
with hoomd.open(name=inFile, mode='rb') as t:
    snap = t[0]
    typ = snap.particles.typeid
    partNum = len(typ)

# Outfile to write data to
base = add + 'hmap_pressure_pa' + str(peA) +\
       '_pb' + str(peB) +\
       '_xa' + str(parFrac) +\
       '_phi' + str(intPhi) +\
       '_ep' + '{0:.3f}'.format(eps) +\
       '_N' + str(partNum)

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
    # Compute each mesh
    NBins = getNBins(l_box, r_cut)
    sizeBin = roundUp((l_box / NBins), 6)
    # Set some values for plotting
    parts = [lat for i in range(0, partNum)]
    myCols = plt.cm.jet
    
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
        tst *= dtau                                 # convert to Brownian time
        
        # Mesh and bin the particles by position
        binParts = [[[] for b in range(NBins)] for a in range(NBins)]
        for k in range(0, partNum):
            # Convert position to be > 0 to place in list mesh
            tmp_posX = pos[k][0] + h_box
            tmp_posY = pos[k][1] + h_box
            x_ind = int(tmp_posX / sizeBin)
            y_ind = int(tmp_posY / sizeBin)
            # Append all particles to appropriate bin
            binParts[x_ind][y_ind].append(k)
            
        # Loop through particles and compute pressure
        pressure = [0. for i in range(0, partNum)]
        for k in range(0, partNum):
            # Get mesh indices
            tmp_posX = pos[k][0] + h_box
            tmp_posY = pos[k][1] + h_box
            x_ind = int(tmp_posX / sizeBin)
            y_ind = int(tmp_posY / sizeBin)
            # Get index of surrounding bins
            l_bin = x_ind - 1  # index of left bins
            r_bin = x_ind + 1  # index of right bins
            b_bin = y_ind - 1  # index of bottom bins
            t_bin = y_ind + 1  # index of top bins
            if r_bin == NBins:
                r_bin -= NBins  # adjust if wrapped
            if t_bin == NBins:
                t_bin -= NBins  # adjust if wrapped
            h_list = [l_bin, x_ind, r_bin]  # list of horizontal bin indices
            v_list = [b_bin, y_ind, t_bin]  # list of vertical bin indices

            # Loop through all bins
            for h in range(0, len(h_list)):
                for v in range(0, len(v_list)):
                    # Take care of periodic wrapping for position
                    wrapX = 0.0
                    wrapY = 0.0
                    if h == 0 and h_list[h] == -1:
                        wrapX -= l_box
                    if h == 2 and h_list[h] == 0:
                        wrapX += l_box
                    if v == 0 and v_list[v] == -1:
                        wrapY -= l_box
                    if v == 2 and v_list[v] == 0:
                        wrapY += l_box
                    # Compute distance between particles
                    for b in range(0, len(binParts[h_list[h]][v_list[v]])):
                        ref = binParts[h_list[h]][v_list[v]][b]
                        dx, dy, r = distComps(pos[k],
                                              pos[ref][0] + wrapX,
                                              pos[ref][1] + wrapY)
                        r = round(r, 4)  # round value to 4 decimal places

                        # If LJ potential is on, compute pressure
                        if 0.1 < r <= r_cut:
                            # Compute the x and y components of force
                            fx, fy = computeFLJ(r, dx, dy, eps)
                            # Compute the x force times x distance
                            sigx = fx * (dx)
                            # Likewise for y
                            sigy = fy * (dy)
                            # Test the code
                            if sigx < 0 or sigy < 0:
                                print("Negative value for force times r!!!")
                            # Add the pressure from this neighbor
                            pressure[k] += ((sigx + sigy) / 2.)
        
        outDPI = 500.
        fig, ax = plt.subplots()
        xy = np.delete(pos, 2, 1)
        coll = matplotlib.collections.EllipseCollection(parts, parts,
                                                        np.zeros_like(parts),
                                                        offsets=xy, units='xy',
                                                        cmap=myCols,
                                                        transOffset=ax.transData)
        coll.set_array(np.ravel(pressure))
#        minCol = min(pressure)
#        maxCol = max(pressure)
#        coll.set_clim([minCol, 1.0])
        ax.add_collection(coll)
        cbar = plt.colorbar(coll, format='%.0f')
        cbar.set_label(r'Interparticle Pressure $(\Pi^{P})$', labelpad=20, rotation=270)

        # Limits and ticks
        viewBuff = 0.0
        ax.set_xlim(-h_box - viewBuff, h_box + viewBuff)
        ax.set_ylim(-h_box - viewBuff, h_box + viewBuff)
        ax.tick_params(axis='both', which='both',
                       bottom=False, top=False, left=False, right=False,
                       labelbottom=False, labeltop=False, labelleft=False, labelright=False)

        pad = str(j).zfill(4)
        ax.set_aspect('equal')
        plt.tight_layout()
        plt.savefig(base + '_' + pad + '.png', dpi=outDPI)
        plt.close()
