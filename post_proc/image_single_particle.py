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
# Enlarge the box to include the periodic images
buff = float(int(r_cut * 2.0) + 1)
# Image rendering options
drawBins = False

# Draw persistence length line on figure
x1 = 1.02
x2 = 1.02
y1 = 0.87
y2 = 0.97

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

    # Gather up positions as you go
    curPos = []
    dPos = []
    tauPos = []
    tauSegs = []
    my_cols = []

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
        
        # Check if timestep is a multiple of tau rot
        rounded = round(tst * 3., 2)
        if (100 * rounded) % 100 == 0:
            tauPos.append([pos[0][0], pos[0][1]])
            # Make line segments for persistence length
            if len(tauPos) > 1:
                # Now let's take care of periodic wrapping
                dx = tauPos[-2][0] - tauPos[-1][0]
                dy = tauPos[-2][1] - tauPos[-1][1]
                dist = np.sqrt( (dx**2) +  (dy**2) )
                # Crossed the periodic edge
                if dist > h_box:
                    # Store the coordinates
                    x = tauPos[-1][0]
                    y = tauPos[-1][1]
                    shiftx = 0.
                    shifty = 0.
                    # Crossed x
                    if np.abs(dx) > h_box:
                        # Crossed right
                        if dx > 0:
                            x += l_box
                            shiftx = -l_box
                        # Crossed left
                        elif dx < 0:
                            x -= l_box
                            shiftx = l_box
                    # Crossed y
                    if np.abs(dy) > h_box:
                        # Crossed top
                        if dy > 0:
                            y += l_box
                            shifty = -l_box
                        # Crossed bottom
                        elif dy < 0:
                            y -= l_box
                            shifty = l_box
                    # Add the wrapped segment to curPos
                    tauSegs.append([(tauPos[-2][0], tauPos[-2][1]), (x, y)])
                    tauSegs.append([(tauPos[-2][0] + shiftx, tauPos[-2][1] + shifty),
                                   (tauPos[-1][0], tauPos[-1][1])])
                
                # Did not cross an edge
                else:
                    tauSegs.append([(tauPos[-2][0], tauPos[-2][1]),
                                    (tauPos[-1][0], tauPos[-1][1])])
            
        tau_segments = mc.LineCollection(tauSegs, lw=1.0, ls="--", colors='r', zorder=1)
        
        # This list holds all unadjusted positions
        dPos.append([pos[0][0], pos[0][1]])
        
        if j > start:
            # Conditional that tests if dPos has gone through an edge
            dx = (dPos[-2][0] - dPos[-1][0])
            dy = (dPos[-2][1] - dPos[-1][1])
            dist = np.sqrt( (dx**2) +  (dy**2) )
            # Jumped the box!
            if dist > h_box:
                # Store the coordinates
                x = dPos[-1][0]
                y = dPos[-1][1]
                shiftx = 0.
                shifty = 0.
                # Crossed x
                if np.abs(dx) > h_box:
                    # Crossed right
                    if dx > 0:
                        x += l_box
                        shiftx = -l_box
                    # Crossed left
                    elif dx < 0:
                        x -= l_box
                        shiftx = l_box
                    
                # Crossed y
                if np.abs(dy) > h_box:
                    # Crossed top
                    if dy > 0:
                        y += l_box
                        shifty = -l_box
                    # Crossed bottom
                    elif dy < 0:
                        y -= l_box
                        shifty = l_box
                # Add the wrapped segment to curPos
                curPos.append([(dPos[-2][0], dPos[-2][1]), (x, y)])
                curPos.append([(dPos[-2][0] + shiftx, dPos[-2][1] + shifty),
                               (dPos[-1][0], dPos[-1][1])])
                my_cols.append(plt.cm.jet((len(dPos) - 1.) / dumps))
                my_cols.append(plt.cm.jet((len(dPos) - 1.) / dumps))
                
            # Did not cross an edge
            else:
                curPos.append([(dPos[-2][0], dPos[-2][1]), (pos[0][0], pos[0][1])])
                my_cols.append(plt.cm.jet((len(dPos) - 1.) / dumps))
                
        # Create line segment collection
        line_segments = mc.LineCollection(curPos, lw=1.0, colors=my_cols, zorder=0)
        # Create frame pad for images
        pad = str(j).zfill(4)
        
        # Plot the figure
        fig, ax = plt.subplots(1, 1)
        ax.add_collection(line_segments)
        ax.add_collection(tau_segments)
        taux = list(list(zip(*tauPos))[0])
        tauy = list(list(zip(*tauPos))[1])
        ax.scatter(taux, tauy, c='r', zorder=2)
        ax.scatter(pos[:,0], pos[:,1], zorder=3)
        ax.text(0.95, 0.025, s=r'$\tau_{r}=$' + '{:0.1f}'.format(tst*3.),
                horizontalalignment='right', verticalalignment='bottom',
                transform=ax.transAxes,
                fontsize=18,
                bbox=dict(facecolor=(1,1,1,0.5), edgecolor=(0,0,0,1), boxstyle='round, pad=0.1'))
        line = lines.Line2D([x1, x2], [y1, y2], ls='--', lw=1.5, c='r', transform=ax.transAxes)
        ax.add_line(line)
        line.set_clip_on(False)
        ax.text(1.05, 0.9, s=r'$l_{p}$',
                transform=ax.transAxes,
                fontsize=18,)
        ax.set_xlim(-h_box, h_box)
        ax.set_ylim(-h_box, h_box)
        ax.axes.set_xticks([])
        ax.axes.set_yticks([])
        ax.axes.set_xticklabels([])
        ax.axes.set_yticks([])
        ax.set_aspect('equal')
        plt.subplots_adjust(0,0,1,1)
        plt.savefig("test_fm" + pad + ".png", dpi=1000)
        plt.close()
