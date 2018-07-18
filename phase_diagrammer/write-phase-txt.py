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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math

# Function that'll grab my parameters from the filenames
def getFromTxt(fname, first, last):
    """Takes a string, text before and after desired text, outs text between"""
    start = fname.index( first ) + len( first )
    end = fname.index( last, start )
    myTxt = fname[start:end]
    return float(myTxt)

sigma = 1.0
# The old model
epsilon = 1.0

# Grab the command line arguments
txtFiles = sys.argv[3:]     # pass starting at 4 to avoid script path and hoomd

# Parse out the activities and fractions from the filenames
peA = np.zeros_like(txtFiles, dtype=np.int)
peB = np.zeros_like(txtFiles, dtype=np.int)
xA = np.zeros_like(txtFiles, dtype=np.float64)
ep = np.zeros_like(txtFiles, dtype=np.int)

# Only if epsilon is in the filename
for i in range(0, len(txtFiles)):
    peA[i] = getFromTxt(txtFiles[i], "pa", "_pb")
    peB[i] = getFromTxt(txtFiles[i], "pb", "_xa")
    xA[i] = getFromTxt(txtFiles[i], "xa", "_ep")
    ep[i] = getFromTxt(txtFiles[i], "ep", ".txt")

try:
    peR = peA.astype(float) / peB.astype(float)         # Compute activity ratio
except:
    peR = np.zeros(len(txtFiles))

# Pairwise potential data (overwritten in this case)
# epsilonA = computeEps(peA)
# epsilonB = computeEps(peB)
epsHS = np.zeros(len(peA), dtype=np.float64)
# for i in range(0, len(peA)): epsHS[i] = epsilonA[i] if epsilonA[i] > epsilonB[i] else epsilonB[i]
epsHS[:] = ep[:] # only if you've set epsilon explicitly
partFracA = xA/100.0    # Particle fraction

# Requirement to be consider phase separated
partNum = 20000
frames = 4000
sizeMin = partNum * 0.4     # 40% of particles in single cluster
timeMin = frames * 0.5      # cluster present for half of all frames

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
    denDen, \
    gasDen, \
    clustAr, \
    gasAr, \
    MCS = np.loadtxt(txtFiles[i], skiprows=1, unpack=True)

    # See if data satisfies phase separated requirements
    count = 0
    phaseSep = 0
    for j in range(0, len(lgClust)):
        if lgClust >= sizeMin:
            count += 1

    if count >= timeMin:
        phasSep = 1
        print("{}: Y").format(txtFiles[i])
    else:
        print("{}: N").format(txtFiles[i])

    # Save to phase separation array









# Write phase seperation array to text file

