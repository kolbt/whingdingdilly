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
MCS = np.zeros(len(txtFiles), dtype=np.ndarray)
sigAll = np.zeros(len(txtFiles), dtype=np.ndarray)
sigAA = np.zeros(len(txtFiles), dtype=np.ndarray)
sigAB = np.zeros(len(txtFiles), dtype=np.ndarray)
sigBB = np.zeros(len(txtFiles), dtype=np.ndarray)
phiEff = np.zeros(len(txtFiles), dtype=np.ndarray)
lgClustA = np.zeros(len(txtFiles), dtype=np.ndarray)
totClustA = np.zeros(len(txtFiles), dtype=np.ndarray)
lgClustD = np.zeros(len(txtFiles), dtype=np.ndarray)
dpDensity = np.zeros(len(txtFiles), dtype=np.ndarray)
gpDensity = np.zeros(len(txtFiles), dtype=np.ndarray)
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
    MCS[i],\
    sigALL[i],\
    sigAA[i],\
    sigAB[i],\
    sigBB[i],\
    phiEff[i],\
    lgClustA[i],\
    totClustA[i],\
    lgClustD[i],\
    dpDensity[i],\
    gpDensity[i] = np.loadtxt(txtFiles[i], skiprows=1, unpack=True)

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

# MCS
for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], MCS[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Mean Cluster Size')
plt.savefig('MCS.png', dpi=1000)
plt.close()

# Largest Cluster Area
for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], lgClustA[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Largest Cluster Area')
plt.savefig('largest_cluster_area.png', dpi=1000)
plt.close()

# Total Cluster Area
for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], totClustA[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Total Cluster Area')
plt.savefig('total_cluster_area.png', dpi=1000)
plt.close()

# Largest Cluster Density
for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], lgClustD[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Largest Cluster Density')
plt.savefig('largest_cluster_density.png', dpi=1000)
plt.close()

# Dense phase density
for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], dpDensity[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Dense Phase Density')
plt.savefig('dpDensity.png', dpi=1000)
plt.close()

# Gas phase density
for i in range(0, len(txtFiles)):
    plt.plot(tsteps[i], gpDensity[i], label=str(peA[i]))
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Gas Phase Density')
plt.savefig('gpDensity.png', dpi=1000)
plt.close()

################################################################################
# Steady state values
################################################################################
ssGasA = np.zeros((len(txtFiles)), dtype=np.float64)
ssGasB = np.zeros((len(txtFiles)), dtype=np.float64)
ssGasAll = np.zeros((len(txtFiles)), dtype=np.float64)
ssDenseA = np.zeros((len(txtFiles)), dtype=np.float64)
ssDenseB = np.zeros((len(txtFiles)), dtype=np.float64)
ssDenseAll = np.zeros((len(txtFiles)), dtype=np.float64)
ssLgClust = np.zeros((len(txtFiles)), dtype=np.float64)
ssDPDensity = np.zeros((len(txtFiles)), dtype=np.float64)
ssGPDensity = np.zeros((len(txtFiles)), dtype=np.float64)
ssDPArea = np.zeros((len(txtFiles)), dtype=np.float64)
ssGPArea = np.zeros((len(txtFiles)), dtype=np.float64)

# Get steady state values (last 100 dumps)
for i in range(0, len(txtFiles)):
    ssGasA[i] = np.mean(gasA[i][-100:-1])
    ssGasB[i] = np.mean(gasB[i][-100:-1])
    ssGasAll[i] = np.mean(gasAll[i][-100:-1])
    ssDenseA[i] = np.mean(denseA[i][-100:-1])
    ssDenseB[i] = np.mean(denseB[i][-100:-1])
    ssDenseAll[i] = np.mean(denseAll[i][-100:-1])
    ssLgClust[i] = np.mean(lgClust[i][-100:-1])
    ssDPDensity[i] = np.mean(dpDensity[i][-100:-1])
    ssGPDensity[i] = np.mean(gpDensity[i][-100:-1])
    ssDPArea[i] = np.mean(dpArea[i][-100:-1])
    ssGPArea[i] = np.mean(gpArea[i][-100:-1])

# Gas Phase #####################

# What is varying? (xA or PeA)
#if xA[0] == xA[1] == 100:
#    plt.scatter(peA, ssGasA, c=peA)
#    plt.xlabel(r'Activity')
#elif xA[0] == xA[1] == 0:
#    plt.scatter(peB, ssGasA, c=peB)
#    plt.xlabel(r'Activity')
#else:
#    ratio = np.zeros_like(peA)
#    ratio = peA / peB
#    plt.scatter(ratio, ssGasA, label=str(ratio), c=ratio)
#    plt.xlabel(r'Activity Ratio')
#
#plt.ylabel('Steady-State Fraction of A-Particles in Gas')
#plt.legend()
##plt.ylim(0,1)
#plt.savefig('SteadyState_gasA.png', dpi=1000)
#plt.close()

for i in range(0, len(txtFiles)):
    # All monodisperse B, use activity
    if xA[i] == xA[i-1] == 0:
        plt.scatter(peB[i], ssGasA[i], c='k')
        plt.xlabel(r'Activity')
    # All monodisperse A, use activity
    elif xA[i] == xA[i-1] == 100:
        plt.scatter(peA[i], ssGasA[i], c='k')
        plt.xlabel(r'Activity')
    # Binary, varying xA
    elif xA[i] != xA[i-1]:
        plt.scatter(xA[i], ssGasA[i], c='k')
        plt.xlabel(r'Particle Fraction $x_{A}$')
        plt.xlim(0,1)
    # Binary, varying peA
    else:
        ratio = float(peA[i] / peB[i])
        plt.scatter(ratio, ssGasA[i], c='k')
        plt.xlabel(r'Activity Ratio')

plt.ylabel('Steady-State Fraction of A-Particles in Gas')
plt.savefig('SteadyState_gasA.png', dpi=1000)
plt.close()

for i in range(0, len(txtFiles)):
    # Monodisperse B, use activity
    if xA[i] == 0:
        plt.scatter(peB[i], ssGasB[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Monodisperse A, use activity
    elif xA[i] == 100:
        plt.scatter(peA[i], ssGasB[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Binary, varying xA
    elif xA[i] != xA[i-1]:
        plt.scatter(xA[i], ssGasB[i], c='k')
        plt.xlabel(r'Particle Fraction $x_{A}$')
        plt.xlim(0,1)
    # Binary, use activity ratio
    else:
        ratio = float(peA[i] / peB[i])
        plt.scatter(ratio, ssGasB[i], label=str(ratio), c='k')
        plt.xlabel(r'Activity Ratio')

plt.ylabel('Steady-State Fraction of B-Particles in Gas')
plt.savefig('SteadyState_gasB.png', dpi=1000)
plt.close()

for i in range(0, len(txtFiles)):
    # Monodisperse B, use activity
    if xA[i] == 0:
        plt.scatter(peB[i], ssGasAll[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Monodisperse A, use activity
    elif xA[i] == 100:
        plt.scatter(peA[i], ssGasAll[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Binary, varying xA
    elif xA[i] != xA[i-1]:
        plt.scatter(xA[i], ssGasAll[i], c='k')
        plt.xlabel(r'Particle Fraction $x_{A}$')
        plt.xlim(0,1)
    # Binary, use activity ratio
    else:
        ratio = float(peA[i] / peB[i])
        plt.scatter(ratio, ssGasAll[i], label=str(ratio), c='k')
        plt.xlabel(r'Activity Ratio')

plt.ylabel('Steady-State Fraction of All Particles in Gas')
plt.savefig('SteadyState_gasAll.png', dpi=1000)
plt.close()

# Dense Phase ###################
for i in range(0, len(txtFiles)):
    # Monodisperse B, use activity
    if xA[i] == 0:
        plt.scatter(peB[i], ssDenseA[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Monodisperse A, use activity
    elif xA[i] == 100:
        plt.scatter(peA[i], ssDenseA[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Binary, varying xA
    elif xA[i] != xA[i-1]:
        plt.scatter(xA[i], ssDenseA[i], c='k')
        plt.xlabel(r'Particle Fraction $x_{A}$')
        plt.xlim(0,1)
    # Binary, use activity ratio
    else:
        ratio = float(peA[i] / peB[i])
        plt.scatter(ratio, ssDenseA[i], label=str(ratio), c='k')
        plt.xlabel(r'Activity Ratio')

plt.ylabel('Steady-State Fraction of A-Particles in Dense Phase')
plt.savefig('SteadyState_denseA.png', dpi=1000)
plt.close()

for i in range(0, len(txtFiles)):
    # Monodisperse B, use activity
    if xA[i] == 0:
        plt.scatter(peB[i], ssDenseB[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Monodisperse A, use activity
    elif xA[i] == 100:
        plt.scatter(peA[i], ssDenseB[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Binary, varying xA
    elif xA[i] != xA[i-1]:
        plt.scatter(xA[i], ssDenseB[i], c='k')
        plt.xlabel(r'Particle Fraction $x_{A}$')
        plt.xlim(0,1)
    # Binary, use activity ratio
    else:
        ratio = float(peA[i] / peB[i])
        plt.scatter(ratio, ssDenseB[i], label=str(ratio), c='k')
        plt.xlabel(r'Activity Ratio')

plt.ylabel('Steady-State Fraction of B-Particles in Dense Phase')
plt.savefig('SteadyState_denseB.png', dpi=1000)
plt.close()

for i in range(0, len(txtFiles)):
    # Monodisperse B, use activity
    if xA[i] == 0:
        plt.scatter(peB[i], ssDenseAll[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Monodisperse A, use activity
    elif xA[i] == 100:
        plt.scatter(peA[i], ssDenseAll[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Binary, varying xA
    elif xA[i] != xA[i-1]:
        plt.scatter(xA[i], ssDenseAll[i], c='k')
        plt.xlabel(r'Particle Fraction $x_{A}$')
        plt.xlim(0,1)
    # Binary, use activity ratio
    else:
        ratio = float(peA[i] / peB[i])
        plt.scatter(ratio, ssDenseAll[i], label=str(ratio), c='k')
        plt.xlabel(r'Activity Ratio')

plt.ylabel('Steady-State Fraction of All Particles in Dense Phase')
plt.savefig('SteadyState_denseAll.png', dpi=1000)
plt.close()

# Largest Cluster ###############
for i in range(0, len(txtFiles)):
    # Monodisperse B, use activity
    if xA[i] == 0:
        plt.scatter(peB[i], ssLgClust[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Monodisperse A, use activity
    elif xA[i] == 100:
        plt.scatter(peA[i], ssLgClust[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Binary, varying xA
    elif xA[i] != xA[i-1]:
        plt.scatter(xA[i], ssLgClust[i], c='k')
        plt.xlabel(r'Particle Fraction $x_{A}$')
        plt.xlim(0,1)
    # Binary, use activity ratio
    else:
        ratio = float(peA[i] / peB[i])
        plt.scatter(ratio, ssLgClust[i], label=str(ratio), c='k')
        plt.xlabel(r'Activity Ratio')

plt.ylabel('Steady-State Largest Cluster Size')
plt.savefig('SteadyState_dpDensity.png', dpi=1000)
plt.close()

# Cluster Density ###############
for i in range(0, len(txtFiles)):
    # Monodisperse B, use activity
    if xA[i] == 0:
        plt.scatter(peB[i], ssDPDensity[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Monodisperse A, use activity
    elif xA[i] == 100:
        plt.scatter(peA[i], ssDPDensity[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Binary, varying xA
    elif xA[i] != xA[i-1]:
        plt.scatter(xA[i], ssDPDensity[i], c='k')
        plt.xlabel(r'Particle Fraction $x_{A}$')
        plt.xlim(0,1)
    # Binary, use activity ratio
    else:
        ratio = float(peA[i] / peB[i])
        plt.scatter(ratio, ssDPDensity[i], label=str(ratio), c='k')
        plt.xlabel(r'Activity Ratio')

plt.ylabel('Steady-State Dense Phase Density')
plt.savefig('SteadyState_dpDensity.png', dpi=1000)
plt.close()

for i in range(0, len(txtFiles)):
    # Monodisperse B, use activity
    if xA[i] == 0:
        plt.scatter(peB[i], ssGPDensity[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Monodisperse A, use activity
    elif xA[i] == 100:
        plt.scatter(peA[i], ssGPDensity[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Binary, varying xA
    elif xA[i] != xA[i-1]:
        plt.scatter(xA[i], ssGPDensity[i], c='k')
        plt.xlabel(r'Particle Fraction $x_{A}$')
        plt.xlim(0,1)
    # Binary, use activity ratio
    else:
        ratio = float(peA[i] / peB[i])
        plt.scatter(ratio, ssGPDensity[i], label=str(ratio), c='k')
        plt.xlabel(r'Activity Ratio')

plt.ylabel('Steady-State Gas Phase Density')
plt.savefig('SteadyState_gpDensity.png', dpi=1000)
plt.close()

# Cluster Area ##################
for i in range(0, len(txtFiles)):
    # Monodisperse B, use activity
    if xA[i] == 0:
        plt.scatter(peB[i], ssDPArea[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Monodisperse A, use activity
    elif xA[i] == 100:
        plt.scatter(peA[i], ssDPArea[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Binary, varying xA
    elif xA[i] != xA[i-1]:
        plt.scatter(xA[i], ssDPArea[i], c='k')
        plt.xlabel(r'Particle Fraction $x_{A}$')
        plt.xlim(0,1)
    # Binary, use activity ratio
    else:
        ratio = float(peA[i] / peB[i])
        plt.scatter(ratio, ssDPArea[i], label=str(ratio), c='k')
        plt.xlabel(r'Activity Ratio')

plt.ylabel('Steady-State Dense Phase Area')
plt.savefig('SteadyState_dpArea.png', dpi=1000)
plt.close()

for i in range(0, len(txtFiles)):
    # Monodisperse B, use activity
    if xA[i] == 0:
        plt.scatter(peB[i], ssGPArea[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Monodisperse A, use activity
    elif xA[i] == 100:
        plt.scatter(peA[i], ssGPArea[i], label=str(peA[i]), c='k')
        plt.xlabel(r'Activity')
    # Binary, varying xA
    elif xA[i] != xA[i-1]:
        plt.scatter(xA[i], ssGPArea[i], c='k')
        plt.xlabel(r'Particle Fraction $x_{A}$')
        plt.xlim(0,1)
    # Binary, use activity ratio
    else:
        ratio = float(peA[i] / peB[i])
        plt.scatter(ratio, ssGPArea[i], label=str(ratio), c='k')
        plt.xlabel(r'Activity Ratio')

plt.ylabel('Steady-State Gas Phase Area')
plt.savefig('SteadyState_gpArea.png', dpi=1000)
plt.close()

