'''
#                           This is an 80 character line                       #
This file reads in processed data from text files and outputs whether or not the
given system undergoes phase separation.  This is currently computed by the size
and lifetime of the largest cluster, but, could also feasibly arise from the
percentage of clustered particles in the system (not necessarily one cluster).

This file:
    1.  Reads in post-processed text data
    2.  Examines the largest cluster size
    3.  Writes to a phase separation file the values: peA, peB, xA, ps
'''

# Imports and loading the .gsd file
import sys
import numpy as np
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
txtFiles = sys.argv[2:]     # start at 2 to avoid early arguments

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
sizeMin = partNum * 0.25     # 25% of particles in single cluster
timeMin = frames * 0.5      # cluster present for half of all frames
# Name to write to
phase_file = "phase-separation-data.txt"
f = open(phase_file, 'w')
f.write(('Act_A').center(10) + ' ' + \
        ('Act_B').center(10) + ' ' + \
        ('Frac_A').center(10) + ' ' + \
        ('Phase_Sep').center(10) + \
        '\n')
f.close()

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
        if lgClust[j] >= sizeMin:
            count += 1


    if count >= timeMin:
        phaseSep = 1
        f = open(phase_file, 'a')
        f.write(str(peA[i]).center(10) + ' ' + \
                str(peB[i]).center(10) + ' ' + \
                str(xA[i]).center(10) + ' ' + \
                str(phaseSep).center(10) + \
                '\n')
    else:
        f = open(phase_file, 'a')
        f.write(str(peA[i]).center(10) + ' ' + \
                str(peB[i]).center(10) + ' ' + \
                str(xA[i]).center(10) + ' ' + \
                str(phaseSep).center(10) + \
                '\n')
    f.close()
