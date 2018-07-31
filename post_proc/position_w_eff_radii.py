'''
#                           This is an 80 character line                       #
This file is intended to look at the average distance between particles (delta).
This should directly result from the forces applied on particles (that's how
potentials work). We will then look at how this changes (for a specific epsilon)
as a function of activity.

What this file does:
    1. Load in source .gsd files (in a list)
    2. Compute center-to-center MODE distance (of last 10 frames)
    3. Index this with the corresponding activity ratio
    4. Plot mode center-to-center distance as a function of activity
'''

# Imports and loading the .gsd file
import sys

import numpy as np
from scipy.interpolate import griddata
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
import math

sigma=1.0
# The old model
epsilon = 1.0

# The new model
def computeEps(activity):
    "Given particle activity, output repulsion well depth"
    epsilon = activity * sigma / 24.0
    return epsilon

# Function that'll grab my parameters from the filenames
def getFromTxt(fname, first, last):
    """Takes a string, text before and after desired text, outs text between"""
    start = fname.index( first ) + len( first )
    end = fname.index( last, start )
    myTxt = fname[start:end]
    return float(myTxt)

# Computes distance
def getDistance(point1, point2x, point2y):
    """"Find the distance between two points"""
    distance = np.sqrt((point2x - point1[0])**2 + (point2y - point1[1])**2)
    return distance

# Computes force given distance with LJ potential
def computeLJForce(r, eps_in):
    """Given a distance, computes the force"""
    forceLJ = 4*eps_in*((12*(sigma**12)*(r**-13))-(6*(sigma**12)*(r**-7)))
    return forceLJ

# Grab the command line arguments
gsdPath = sys.argv[3]
gsdFiles = sys.argv[4:]     # pass starting at 4 to avoid script path and hoomd
sys.path.append(gsdPath)
import gsd
from gsd import hoomd

# Parse out the activities and fractions from the filenames
peA = np.zeros_like(gsdFiles, dtype=np.int)
peB = np.zeros_like(gsdFiles, dtype=np.int)
xA = np.zeros_like(gsdFiles, dtype=np.float64)
ep = np.zeros_like(gsdFiles, dtype=np.int)

# This grabs the parameters of each text file
# for i in range(0, len(gsdFiles)):
#    peA[i] = getFromTxt(gsdFiles[i], "pa", "_pb")
#    peB[i] = getFromTxt(gsdFiles[i], "pb", "_xa")
#    xA[i] = getFromTxt(gsdFiles[i], "xa", ".gsd")

# Only if epsilon is in the filename
for i in range(0, len(gsdFiles)):
    peA[i] = getFromTxt(gsdFiles[i], "pa", "_pb")
    peB[i] = getFromTxt(gsdFiles[i], "pb", "_xa")
    xA[i] = getFromTxt(gsdFiles[i], "xa", "_ep")
    ep[i] = getFromTxt(gsdFiles[i], "ep", ".gsd")

peR = peA.astype(float) / peB.astype(float)         # Compute activity ratio

epsilonA = computeEps(peA)
epsilonB = computeEps(peB)
epsHS = np.zeros(len(peA), dtype=np.float64)
for i in range(0, len(peA)): epsHS[i] = epsilonA[i] if epsilonA[i] > epsilonB[i] else epsilonB[i]
# epsHS[:] = ep[:] # only if you've set epsilon explicitly
epsHS[:] = 1

partFracA = xA/100.0    # Particle fraction
mode = []               # List to store modes in

lowest = 1.122 # get the smallest analyzed diameter

for i in range(0, len(gsdFiles)):
    print('Computing mode center-to-center distance for data in file: {}'.format(gsdFiles[i]))
    # Load in the data (read, binary)
    f = hoomd.open(name=gsdFiles[i], mode='rb') # open gsd file with hoomd
    dumps = f.__len__()                         # get number of frames
    # Start and stop frames
    start = dumps - 1   # gives first frame to read
    end = dumps         # gives last frame to read
    # Instantiate necessary arrays
    positions = np.zeros((end), dtype=np.ndarray)   # array of positions
    types = np.zeros((end), dtype=np.ndarray)       # particle types
    box_data = np.zeros((1), dtype=np.ndarray)      # box dimensions

    # Get relevant data from .gsd file
    with hoomd.open(name=gsdFiles[i], mode='rb') as t:
        # Get the box dimensions
        snap = t[0]
        box_data = snap.configuration.box
        # Load data for each timestep
        for j in range(start, end):
            snap = t[j]                             # snapshot of frame
            types[j] = snap.particles.typeid        # get types
            positions[j] = snap.particles.position  # get positions

    # Get number of each type of particle
    partNum = len(types[start])
    partA = int(partNum * partFracA[i])
    partB = partNum - partA

    # Get the relevant box values
    l_box = box_data[0]
    h_box = l_box / 2.0
    a_box = l_box * l_box

    # Make the mesh
    r_cut = 1.122
    sizeBin = r_cut
    nBins = int(l_box / sizeBin)
    nBins += 1  # account for integer rounding

    # Instantiate lists here to sum temporally
    ALL = []
    AA = []
    BB = []
    AB = []

    for j in range(start, end):

        # Mesh array
        binParts = [[[] for b in range(nBins)] for a in range(nBins)]

        # Easier accessors
        pos = positions[j]
        pos = np.delete(pos, 2, 1)
        typ = types[j]

        effRad = np.zeros((partNum), dtype=np.float64)

        # Put particles in their respective bins
        for k in range(0, partNum):
            # Get mesh indices
            tmp_posX = pos[k][0] + h_box
            tmp_posY = pos[k][1] + h_box
            x_ind = int(tmp_posX / sizeBin)
            y_ind = int(tmp_posY / sizeBin)
            # Append particle id to appropriate bin
            binParts[x_ind][y_ind].append(k)

        # Compute distance, each pair will be counted twice
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
            if r_bin == nBins:
                r_bin -= nBins  # adjust if wrapped
            if t_bin == nBins:
                t_bin -= nBins  # adjust if wrapped
            h_list = [l_bin, x_ind, r_bin]  # list of horizontal bin indices
            v_list = [b_bin, y_ind, t_bin]  # list of vertical bin indices


            neighCount = 0
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
                        r = getDistance(pos[k],
                                        pos[ref][0] + wrapX,
                                        pos[ref][1] + wrapY)
                        r = round(r, 4)  # round value to 4 decimal places

                        # If LJ potential is on, store into a list (omit self)
                        if 0.1 < r <= r_cut:
                            neighCount += 1
                            effRad[k] += r
            if neighCount != 0:
                effRad[k] /= neighCount
            else:
                effRad[k] = 1.122
            # Get lowest value (for heatmap)
            if effRad[k] < lowest:
                lowest = effRad[k]

        # print("Closest center-to-center distance: {}").format(lowest)
        # Plot radius overlaid on position
        x = np.zeros((partNum), dtype=np.float64)
        y = np.zeros((partNum), dtype=np.float64)
        for w in range(0, partNum):
            x[w] = pos[w][0]
            y[w] = pos[w][1]
        plt.scatter(x, y, c=effRad, s=0.25, edgecolors='none')
        plt.xlim(-h_box, h_box)
        plt.ylim(-h_box, h_box)
        plt.xticks([])
        plt.yticks([])
        plt.clim(lowest, 1.0)
        plt.colorbar()
        plt.tight_layout()
        plt.axes().set_aspect('equal')
        plt.savefig('effSigma' +\
                    '_pa' + str(peA[i]) +\
                    '_pb' + str(peB[i]) +\
                    '_xa' + str(int(xA[i])) +\
                    '_ep' + str(int(epsHS[i])) +\
                    '_fm' + str(j) +\
                    '.png', dpi=1000)
        plt.close()
