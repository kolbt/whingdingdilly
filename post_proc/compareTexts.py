'''
#                           This is an 80 character line                       #
Purpose: read data from multiple text files and compare (plot) in one convenient
image :) Wow, go you
'''

import sys
import numpy as np
import matplotlib.pyplot as plt

# Use this to back-calculate time in units of Brownian tau
kT = 1.0                        # temperature
threeEtaPiSigma = 1.0           # drag coefficient
sigma = 1.0                     # particle diameter
D_t = kT / threeEtaPiSigma      # translational diffusion constant
D_r = (3.0 * D_t) / (sigma**2)  # rotational diffusion constant
tauBrown = (sigma**2) / D_t     # brownian time scale (invariant)

# Define what functions you'll need here
def getFromTxt(fname, first, last):
    "Takes a string, text before and after desired text, outs text between"
    start = fname.index( first ) + len( first )
    end = fname.index( last, start )
    myTxt = fname[start:end]
    return float(myTxt)
# Above function kindly provided by user "cji" on stackoverflow
#https://stackoverflow.com/questions/3368969/find-string-between-two-substrings

def computeVel(activity):
    "Given particle activity, output intrinsic swim speed"
    velocity = (activity * sigma) / (3 * (1/D_r))
    return velocity

def computeActiveForce(velocity):
    "Given particle activity, output repulsion well depth"
    activeForce = velocity * threeEtaPiSigma
    return activeForce

def computeEps(activeForce):
    "Given particle activity, output repulsion well depth"
    epsilon = activeForce * sigma / 24.0
    return epsilon

def computeTauLJ(epsilon):
    "Given epsilon, compute lennard-jones time unit"
    tauLJ = ((sigma**2) * threeEtaPiSigma) / epsilon
    return tauLJ

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
dt = np.zeros(len(txtFiles), dtype=np.float64)
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
    gasA[i] /= partNum * partFracA[i]
    gasB[i] /= partNum * (1-partFracA[i])
    gasAll[i] /= partNum
    denseA[i] /= partNum * partFracA[i]
    denseB[i] /= partNum * (1-partFracA[i])
    denseAll[i] /= partNum
    lgClust[i] /= partNum

# Compute Brownian time
for i in range(0, len(txtFiles)):
    # Compute parameters from activities
    if peA[i] != 0:                     # A particles are NOT Brownian
        vA = computeVel(peA[i])
        FpA = computeActiveForce(vA)
        epsA = computeEps(FpA)
        tauA = computeTauLJ(epsA)
    else:                               # A particles are Brownian
        vA = 0.0
        FpA = 0.0
        epsA = kT
        tauA = computeTauLJ(epsA)

    if peB[i] != 0:                     # B particles are NOT Brownian
        vB = computeVel(peB[i])
        FpB = computeActiveForce(vB)
        epsB = computeEps(FpB)
        tauB = computeTauLJ(epsB)
    else:                               # B particles are Brownian
        vB = 0.0
        FpB = 0.0
        epsB = kT
        tauB = computeTauLJ(epsB)

    # Get adjusted time units
    tauLJ = (tauA if (tauA <= tauB) else tauB)  # use the smaller tauLJ
    ratio = tauLJ / tauBrown                    # get rid of LJ units
    dt[i] = ratio * 0.00001                     # tstep size
    tsteps[i] *= dt[i]

# Now plot everything

# Gas phase
for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], gasA[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Fraction in Gas Phase')
plt.ylim(0,1)
plt.savefig('gasPhaseA.png', dpi=1000)
plt.close()

for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], gasB[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Fraction in Gas Phase')
plt.ylim(0,1)
plt.savefig('gasPhaseB.png', dpi=1000)
plt.close()

for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], gasAll[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Fraction in Gas Phase')
plt.ylim(0,1)
plt.savefig('gasPhaseAll.png', dpi=1000)
plt.close()

# Dense phase
for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], denseA[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Fraction in Dense Phase')
plt.ylim(0,1)
plt.savefig('densePhaseA.png', dpi=1000)
plt.close()

for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], denseB[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Fraction in Dense Phase')
plt.ylim(0,1)
plt.savefig('densePhaseB.png', dpi=1000)
plt.close()

for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], denseAll[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Fraction in Dense Phase')
plt.ylim(0,1)
plt.savefig('densePhaseAll.png', dpi=1000)
plt.close()

# Largest Cluster
for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], lgClust[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Fraction in Largest Cluster')
plt.ylim(0,1)
plt.savefig('largestCluster.png', dpi=1000)
plt.close()

# Cluster Density
for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], gpDensity[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Gas Phase Density')
#plt.ylim(0,1)
plt.savefig('gpDensity.png', dpi=1000)
plt.close()

for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], dpDensity[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Dense Phase Density')
#plt.ylim(0,1)
plt.savefig('dpDensity.png', dpi=1000)
plt.close()

# Cluster Area
for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], gpArea[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Gas Phase Area')
#plt.ylim(0,1)
plt.savefig('gpArea.png', dpi=1000)
plt.close()

for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], dpArea[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Dense Phase Area')
#plt.ylim(0,1)
plt.savefig('dpArea.png', dpi=1000)
plt.close()

# Steady state values








