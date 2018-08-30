'''
#                           This is an 80 character line                       #

This file:
    1.  Reads in data for diameter from textfiles
    2.  Computes the LJ force for given distances
    3.  Plots this data
'''

# Imports and loading the .gsd file
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# Function that'll grab my parameters from the filenames
def getFromTxt(fname, first, last):
    '''Takes a string, text before and after desired text, outs text between'''
    start = fname.index( first ) + len( first )
    end = fname.index( last, start )
    myTxt = fname[start:end]
    return float(myTxt)

def ljForce(r):
    '''Takes distance gives force from Lennard-Jones potential'''
    forceLJ = 4 * epsilon * ((12 * (sigma ** 12) * (r ** -13)) - (6 * (sigma ** 12) * (r ** -7)))
    return forceLJ

sigma = 1.0
# The old model
epsilon = 1.0

# Grab the command line arguments
txtFiles = sys.argv[2:]     # start at 2 to avoid early arguments

# Parse out the activities and fractions from the filenames
peA = np.zeros_like(txtFiles, dtype=np.int)
peB = np.zeros_like(txtFiles, dtype=np.int)
xA = np.zeros_like(txtFiles, dtype=np.float64)
ep = np.zeros_like(txtFiles, dtype=np.int)

try:
    for i in range(0, len(txtFiles)):
        peA[i] = getFromTxt(txtFiles[i], "pa", "_pb")
        peB[i] = getFromTxt(txtFiles[i], "pb", "_xa")
        xA[i] = getFromTxt(txtFiles[i], "xa", "_ep")
        ep[i] = getFromTxt(txtFiles[i], "ep", ".txt")
except:
    for i in range(0, len(txtFiles)):
        peA[i] = getFromTxt(txtFiles[i], "pa", "_pb")
        peB[i] = getFromTxt(txtFiles[i], "pb", "_xa")
        xA[i] = getFromTxt(txtFiles[i], "xa", ".txt")
        ep[i] = 1

try:
    peR = peA.astype(float) / peB.astype(float)         # Compute activity ratio
except:
    peR = np.zeros(len(txtFiles))

# Instantiate arrays I'd like to plot
phaseSep = np.zeros(len(txtFiles), dtype=np.int)
ALL = np.zeros(len(txtFiles), dtype=np.float64)
AA = np.zeros(len(txtFiles), dtype=np.float64)
AB = np.zeros(len(txtFiles), dtype=np.float64)
BB = np.zeros(len(txtFiles), dtype=np.float64)
phiAvg = np.zeros(len(txtFiles), dtype=np.float64)

# Loop through each data series
for i in range(0, len(txtFiles)):

    # Import data into arrays
    tst, \
    gasA, \
    gasB, \
    gasTot, \
    denA, \
    denB, \
    denTot, \
    lgClust, \
    MCS, \
    sigALL, \
    sigAA, \
    sigAB, \
    sigBB, \
    phiEff, \
    lC_Area, \
    totC_Area, \
    lC_density, \
    denDen, \
    gasDen = np.loadtxt(txtFiles[i], skiprows=1, unpack=True)

    # Requirement to be consider phase separated
    partNum = gasTot[0]     # everything starts in a gas
    frames = len(tst)
    sizeMin = partNum * 0.25  # 40% of particles in single cluster
    timeMin = frames * 0.50  # cluster present for half of all frames

    count = 0

    # Get last 10% of simulation
    numAvg = (0.10 * len(lgClust))
    avgTime = len(lgClust) - numAvg

    for j in range(0, len(lgClust)):
        # Average over last
        if j >= avgTime:
            ALL[i] += sigALL[j]
            AA[i] += sigAA[j]
            AB[i] += sigAB[j]
            BB[i] += sigBB[j]
            phiAvg[i] += phiEff[j]
        if lgClust[j] >= sizeMin:
            count += 1

    # Average diameter values
    ALL[i] /= numAvg
    AA[i] /= numAvg
    AB[i] /= numAvg
    BB[i] /= numAvg
    phiAvg[i] /= numAvg

    if count >= timeMin:
        phaseSep = 1

# Now everything is in an array, sort them (for lines)
for i in range(0, len(txtFiles)):
    for j in range(0, len(txtFiles)):
        # Values need to be swapped
        if peA[i] > peA[j] and i < j:
            # Swap A activity
            tmp = peA[j]
            peA[j] = peA[i]
            peA[i] = tmp
            # Swap total diameter
            tmp = ALL[j]
            ALL[j] = ALL[i]
            ALL[i] = tmp
            # Swap AA diameter
            tmp = AA[j]
            AA[j] = AA[i]
            AA[i] = tmp
            # Swap AB diameter
            tmp = AB[j]
            AB[j] = AB[i]
            AB[i] = tmp
            # Swap BB diameter
            tmp = BB[j]
            BB[j] = BB[i]
            BB[i] = tmp
            # Swap phi
            tmp = phiAvg[j]
            phiAvg[j] = phiAvg[i]
            phiAvg[i] = tmp

# Plot the data
plt.plot(peA, ALL, marker='o', c='k', label='Emergent Diameter')

# Get theory on fine r-scale
ys = np.arange(min(ALL), max(ALL), 0.001)
xs = np.zeros_like(ys)
for i in range(0, len(xs)):
    xs[i] = ljForce(ys[i])

# Plot theory
plt.plot(xs, ys, c='g', label=r'$Pe=F_{LJ}$')
plt.plot(0.5 * xs, ys, c='r', label=r'$2Pe=F_{LJ}$')
plt.plot(0.3 * xs, ys, c='b', label=r'$3Pe=F_{LJ}$')

# Axes limits
plt.xlim(min(peA), max(peA))
plt.ylim(min(ALL), max(ALL))
# Labels
plt.xlabel(r'Activity $(Pe)$')
plt.ylabel(r'Center-to-center Distance $(\sigma_{Eff}$)')
# Get information for legend
plt.legend()
# Plot :)
plt.savefig('data_LJ_overlay_monodisperse.png', bbox_inches='tight', dpi=1000)
plt.close()

# # This is an example that works for using multiple axes
# # Instantiate figure
# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# ax2 = ax1.twiny()
# fig.subplots_adjust(bottom=0.2)
#
# # Plot the data
# data = ax1.plot(peA, ALL, c='k', label='Emergent Diameter')
#
# # Get theory on fine r-scale
# ys = np.arange(min(ALL), max(ALL), 0.001)
# xs = np.zeros_like(ys)
# for i in range(0, len(xs)):
#     xs[i] = ljForce(ys[i])
#
# # Plot theory
# first = ax1.plot(xs, ys, c='g', label=r'$Pe=F_{LJ}$')
# second = ax2.plot(xs, ys, c='r', label=r'$2Pe=F_{LJ}$')
# third = ax1.plot(0.45 * xs, ys, c='b', label=r'$Pe2=F_{LJ}$')
#
# # Additional plot restrictions
#
# # Axes limits
# ax1.set_xlim(min(peA), max(peA))
# ax2.set_xlim(2 * min(peA), 2 * max(peA))
# plt.ylim(min(ALL), max(ALL))
# # Labels
# ax1.set_xlabel(r'Activity $(Pe)$')
# ax2.set_xlabel(r'Twice Activity $(2Pe)$')
# plt.ylabel(r'Center-to-center Distance $(\sigma_{Eff}$)')
# # Move second axis to bottom
# ax2.xaxis.set_ticks_position("bottom")
# ax2.xaxis.set_label_position("bottom")
# ax2.spines["bottom"].set_position(("axes", -0.15))
# # Get information for legend
# lns = data + first + second + third
# labs = [l.get_label() for l in lns]
# plt.legend(lns, labs)
# # Plot :)
# plt.show()










