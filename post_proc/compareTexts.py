'''
#                           This is an 80 character line                       #
'''

import sys
import numpy as np
import matplotlib.pyplot as plt

# Define what functions you'll need here
def getFromTxt(fname, first, last):
    "Takes a string, text before and after desired text, outs text between"
    start = fname.index( first ) + len( first )
    end = fname.index( last, start )
    myTxt = fname[start:end]
    return float(myTxt)
# Above function kindly provided by user "cji" on stackoverflow
#https://stackoverflow.com/questions/3368969/find-string-between-two-substrings

# Grab the textfiles to run on
txtFiles = sys.argv[2:] # pass starting at 2 to avoid script path

# Parse out the activities and fractions from the filenames
peA = np.zeros_like(txtFiles, dtype=np.float64)
peB = np.zeros_like(txtFiles, dtype=np.float64)
xA = np.zeros_like(txtFiles, dtype=np.float64)

# This grabs the parameters of each text file
for i in range(0, len(txtFiles)):
    peA[i] = getFromTxt(txtFiles[i], "pa", "_pb")
    peB[i] = getFromTxt(txtFiles[i], "pb", "_xa")
    xA[i] = getFromTxt(txtFiles[i], "xa", ".txt")

partFracA = xA/100.0

# Initialize the arrays so that all text files fit in each one
tsteps = np.zeros(len(txtFiles), dtype=np.ndarray)
gasA = np.zeros(len(txtFiles), dtype=np.ndarray)
gasB = np.zeros(len(txtFiles), dtype=np.ndarray)
gasAll = np.zeros(len(txtFiles), dtype=np.ndarray)
denseA = np.zeros(len(txtFiles), dtype=np.ndarray)
denseB = np.zeros(len(txtFiles), dtype=np.ndarray)
denseAll = np.zeros(len(txtFiles), dtype=np.ndarray)
lgClust = np.zeros(len(txtFiles), dtype=np.ndarray)
dpDensity = np.zeros(len(txtFiles), dtype=np.ndarray)
gpDensity = np.zeros(len(txtFiles), dtype=np.ndarray)
dpArea = np.zeros(len(txtFiles), dtype=np.ndarray)
gpArea = np.zeros(len(txtFiles), dtype=np.ndarray)
# It is absolutley archaic that I have to initialize these like this ^

# Pull data from text files
for i in range(0, len(txtFiles)):
    tsteps[i],\
    gasA[i],\
    gasB[i],\
    gasAll[i],\
    denseA[i],\
    denseB[i],\
    denseAll[i],\
    lgClust[i],\
    dpDensity[i],\
    gpDensity[i],\
    dpArea[i],\
    gpArea[i] = np.loadtxt(txtFiles[i], skiprows=1, unpack=True)

    partNum = gasAll[i][0]
#    tsteps *= dt
    gasA[i] /= partNum * partFracA[i]
    gasB[i] /= partNum * (1-partFracA[i])
    gasAll[i] /= partNum
    denseA[i] /= partNum * partFracA[i]
    denseB[i] /= partNum * (1-partFracA[i])
    denseAll[i] /= partNum
    lgClust[i] /= partNum

# Compute Brownian time
#time = np.arange(0, 600, 0.5)

# Gas phase
for i in range(0, len(txtFiles)):
    time = np.arange(0, len(gasA[i])/2.0, 0.5)
    plt.plot(time, gasA[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Fraction in Gas Phase')
plt.ylim(0,1)
plt.show()
#plt.savefig('gasPhase_' + out + '.png', dpi=1000)
#plt.close()






