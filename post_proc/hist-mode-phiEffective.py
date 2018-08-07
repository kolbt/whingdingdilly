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
for i in range(0, len(gsdFiles)):
   peA[i] = getFromTxt(gsdFiles[i], "pa", "_pb")
   peB[i] = getFromTxt(gsdFiles[i], "pb", "_xa")
   xA[i] = getFromTxt(gsdFiles[i], "xa", ".gsd")

# Only if epsilon is in the filename
# for i in range(0, len(gsdFiles)):
#     peA[i] = getFromTxt(gsdFiles[i], "pa", "_pb")
#     peB[i] = getFromTxt(gsdFiles[i], "pb", "_xa")
#     xA[i] = getFromTxt(gsdFiles[i], "xa", "_ep")
#     ep[i] = getFromTxt(gsdFiles[i], "ep", ".gsd")

peR = peA.astype(float) / peB.astype(float)         # Compute activity ratio

epsilonA = computeEps(peA)
epsilonB = computeEps(peB)
epsHS = np.ones(len(peA), dtype=np.float64)
# for i in range(0, len(peA)): epsHS[i] = epsilonA[i] if epsilonA[i] > epsilonB[i] else epsilonB[i]
# epsHS[:] = ep[:] # only if you've set epsilon explicitly

partFracA = xA/100.0    # Particle fraction
mode = []               # List to store modes in

hist_file = 'effective_particle_radius.txt'
f = open(hist_file, 'w') # write file headings
f.write('PeA'.center(10) + ' ' +\
        'PeB'.center(10) + ' ' +\
        'PeRatio'.center(10) + ' ' +\
        'xSlow'.center(10) + ' ' +\
        'Epsilon'.center(10) + ' ' +\
        'ModeAll'.center(10) + ' ' +\
        'ModeAA'.center(10) + ' ' +\
        'ModeAB'.center(10) + ' ' +\
        'ModeBB'.center(10) + ' ' +\
        'ForceAll'.center(10) + ' ' +\
        'ForceAA'.center(10) + ' ' +\
        'ForceAB'.center(10) + ' ' +\
        'ForceBB'.center(10) + ' ' +\
        'PhiEff'.center(10) + '\n')
f.close()

for i in range(0, len(gsdFiles)):
    print('Computing mode center-to-center distance for data in file: {}'.format(gsdFiles[i]))
    # Load in the data (read, binary)
    f = hoomd.open(name=gsdFiles[i], mode='rb') # open gsd file with hoomd
    dumps = f.__len__()                         # get number of frames
    # Start and stop frames
    start = dumps - 10  # gives first frame to read
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
                            ALL.append(r)  # All particles
                            if typ[k] == 0 and typ[ref] == 0:  # AA distance
                                AA.append(r)
                            elif typ[k] == 1 and typ[ref] == 1:  # BB distance
                                BB.append(r)
                            else:  # AB distance
                                AB.append(r)

    # Compute quantities for specific simulation, write to text file

    # ALL
    modeALL = stats.mode(ALL)
    modeALL = round(modeALL[0][0], 4)
    fALL = computeLJForce(modeALL, epsHS[i])
    phiEff = modeALL * 0.6                      # Effective area fraction phi=0.6

    # AA
    try:
        modeAA = stats.mode(AA)
        modeAA = round(modeAA[0][0], 4)
        fAA = computeLJForce(modeAA, epsHS[i])
    except:
        modeAA = 0.0
        fAA = 0.0

    # AB
    try:
        modeAB = stats.mode(AB)
        modeAB = round(modeAB[0][0], 4)
        fAB = computeLJForce(modeAB, epsHS[i])
    except:
        modeAB = 0.0
        fAB = 0.0

    # BB
    try:
        modeBB = stats.mode(BB)
        modeBB = round(modeBB[0][0], 4)
        fBB = computeLJForce(modeBB, epsHS[i])
    except:
        modeBB = 0.0
        fBB = 0.0

    # Monodisperse? Done with exceptions above
    # modeAB = 0
    # modeBB = 0
    # fAB = 0
    # fBB = 0

    # Write simulation data to textfile
    f = open(hist_file, 'a')
    f.write(str(peA[i]).center(10) + ' ' + \
            str(peB[i]).center(10) + ' ' + \
            str(peR[i]).center(10) + ' ' + \
            str(xA[i]).center(10) + ' ' + \
            '{0:.2f}'.format(epsHS[i]).center(10) + ' ' + \
            '{0:.3f}'.format(modeALL).center(10) + ' ' + \
            '{0:.3f}'.format(modeAA).center(10) + ' ' + \
            '{0:.3f}'.format(modeAB).center(10) + ' ' + \
            '{0:.3f}'.format(modeBB).center(10) + ' ' + \
            '{0:.0f}'.format(fALL).center(10) + ' ' + \
            '{0:.0f}'.format(fAA).center(10) + ' ' + \
            '{0:.0f}'.format(fAB).center(10) + ' ' + \
            '{0:.0f}'.format(fBB).center(10) + ' ' + \
            '{0:.3f}'.format(phiEff).center(10) + \
            '\n')
    f.close()

# # Compute LJ values from mode values to get experienced force
# ljForce = np.zeros(len(mode), dtype=np.float64)
# phiEff = np.zeros(len(mode), dtype=np.float64)
# for i in range(0, len(mode)):
#     ljForce[i] = computeLJForce(mode[i], epsHS[i])  # switch to 'epsilon' if old method
#     #ljForce[i] = computeLJForce(mode[i], epsilon)   # this is epsilon = 1 case
#     phiEff[i] = mode[i] * 0.6   # Compute effective are fraction from the mode
#
# # If activity A is non-zero
# if (any(peA)):
#
#     # Need a plot with two axes
#     fig = plt.figure(facecolor='w', edgecolor='k', frameon=True)
#     ax1 = fig.add_subplot(111)
#     ax2 = ax1.twinx()
#     ax3 = ax1.twinx()
#     # Label y-axes and format colors
#     ax1.set_ylabel(r'Effective Diameter $(\sigma)$')
#     ax1.yaxis.label.set_color('b')
#     ax2.set_ylabel(r'Corresponding $F_{LJ}$')
#     ax2.yaxis.label.set_color('k')
#     ax3.set_ylabel(r'$\phi_{Effective}$')
#     ax3.yaxis.label.set_color('b')
#     # Data for second axis
#     ax2.plot(peA, ljForce, 'k.')
#     # Invert the axis
#     ylims = ax2.get_ylim()
#     ax2.set_ylim([ylims[1], ylims[0]])
#     # Set ticks
#     ax1.tick_params('y', colors='b')
#     ax2.tick_params('y', colors='k')
#     ax3.tick_params('y', colors='b')
#     # Turn off gridlines
#     plt.setp(ax1.get_yticklabels(), visible=True)
#     plt.setp(ax2.get_yticklabels(), visible=True)
#
#     # Data for first axis
#     ax1.plot(peA, mode, 'b.')
#     ax3.plot(peA, phiEff, 'b.')
#     # ax3.set_yticks(np.linspace(ax3.get_yticks()[0], ax3.get_yticks()[-1], len(ax1.get_yticks())))
#     ax3.spines['left'].set_position(('axes', -0.2))
#     ax3.yaxis.set_label_position('left')
#     ax3.yaxis.set_ticks_position('left')
#     plt.xlabel('Activity A')
#     plt.xlim(min(peA), max(peA))
#     plt.savefig('peA_vs_sigma.png', facecolor='w', edgecolor='k', frameon=True, bbox_inches='tight', dpi=1000)
#     plt.close()
#
# # If activity B is non-zero, plot ratio
# if (any(peB)):
#
#     # Need a plot with two axes
#     fig = plt.figure(facecolor='w', edgecolor='k', frameon=True)
#     ax1 = fig.add_subplot(111)
#     ax2 = ax1.twinx()
#     ax3 = ax1.twinx()
#     # Label y-axes and format colors
#     ax1.set_ylabel(r'Effective Diameter $(\sigma)$')
#     ax1.yaxis.label.set_color('b')
#     ax2.set_ylabel(r'Corresponding $F_{LJ}$')
#     ax2.yaxis.label.set_color('k')
#     ax3.set_ylabel(r'$\phi_{Effective}$')
#     ax3.yaxis.label.set_color('b')
#     # Data for second axis
#     peR = peA.astype(float) / peB.astype(float)
#     ax2.plot(peR, ljForce, 'k.')
#     # Invert second y axis
#     ylims = ax2.get_ylim()
#     ax2.set_ylim([ylims[1], ylims[0]])
#     # Set ticks
#     ax1.tick_params('y', colors='b')
#     ax2.tick_params('y', colors='k')
#     ax3.tick_params('y', colors='b')
#     # Turn off gridlines
#     plt.setp(ax1.get_yticklabels(), visible=True)
#     plt.setp(ax2.get_yticklabels(), visible=True)
#
#     # Data for first axis
#     ax1.plot(peR, mode, 'b.')
#     ax3.plot(peR, phiEff, 'b.')
#     # Move the axis to the left
#     ax3.spines['left'].set_position(('axes', -0.2))
#     ax3.yaxis.set_label_position('left')
#     ax3.yaxis.set_ticks_position('left')
#     plt.xlabel('Activity Ratio')
#     plt.xlim(min(peR), max(peR))
#     plt.savefig('peRatio_vs_sigma.png', facecolor='w', edgecolor='k', frameon=True, bbox_inches='tight', dpi=1000)
#     plt.close()
#
# # If particle fraction is varied
# if (any(xA - xA[0])):
#     # Make plot with 3 axes
#     fig = plt.figure(facecolor='w', edgecolor='k', frameon=True)
#     ax1 = fig.add_subplot(111)
#     ax2 = ax1.twinx()
#     ax3 = ax1.twinx()
#     # Label y-axes and format colors
#     ax1.set_ylabel(r'Effective Diameter $(\sigma)$')
#     ax1.yaxis.label.set_color('b')
#     ax2.set_ylabel(r'Corresponding $F_{LJ}$')
#     ax2.yaxis.label.set_color('k')
#     ax3.set_ylabel(r'$\phi_{Effective}$')
#     ax3.yaxis.label.set_color('b')
#     # Data for second axis
#     ax2.plot(xA, ljForce, 'k.')
#     # Invert second y axis
#     ylims = ax2.get_ylim()
#     ax2.set_ylim([ylims[1], ylims[0]])
#     # Set ticks
#     ax1.tick_params('y', colors='b')
#     ax2.tick_params('y', colors='k')
#     ax3.tick_params('y', colors='b')
#     # Turn off gridlines
#     plt.setp(ax1.get_yticklabels(), visible=True)
#     plt.setp(ax2.get_yticklabels(), visible=True)
#
#     # Data for first axis
#     ax1.plot(xA, mode, 'b.')
#     ax3.plot(xA, phiEff, 'b.')
#     # Move the axis to the left
#     ax3.spines['left'].set_position(('axes', -0.2))
#     ax3.yaxis.set_label_position('left')
#     ax3.yaxis.set_ticks_position('left')
#     plt.xlabel('Particle Fraction')
#     plt.xlim(min(xA), max(xA))
#     plt.savefig('xA_vs_sigma.png', facecolor='w', edgecolor='k', frameon=True, bbox_inches='tight', dpi=1000)
#     plt.close()
