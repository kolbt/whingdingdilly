'''
#                           This is an 80 character line                       #
INTENT: take text file data and plot it.  Saving both the text files and the
        plots allows for simple replotting without recalculating.
        
Note: I'll need to generate another file which compares between simulations

Column:   1   2-4(A,B,Tot)  5-7(A,B,Tot)      8     9-10(DP,GP)  11-12(DP,GP)
Data:   Time     gasPop       densePop    Lg_clust    Density        Area

    ( 1.) Read data into arrays
    ( 2.) Make plots of each
    ( 3.) Celebrate with a cup of coffee :)
'''

# Import the goods
import sys
import matplotlib.pyplot as plt
import numpy as np

# Get parameters passed in from filename
peA = int(sys.argv[1])                      # activity A
peB = int(sys.argv[2])                      # activity B
partPercA = int(sys.argv[3])                # percentage A particles
partFracA = float(partPercA) / 100.0        # fraction A particles

# File to pull data from
inTxt = "all_pa"+str(peA)+\
"_pb"+str(peB)+\
"_xa"+str(partPercA)+\
".txt"

# Base name to write to
out = "pa"+str(peA)+\
"_pb"+str(peB)+\
"_xa"+str(partPercA)

# Use this to back-calculate time in units of Brownian tau
kT = 1.0                        # temperature
threeEtaPiSigma = 1.0           # drag coefficient
sigma = 1.0                     # particle diameter
D_t = kT / threeEtaPiSigma      # translational diffusion constant
D_r = (3.0 * D_t) / (sigma**2)  # rotational diffusion constant
tauBrown = (sigma**2) / D_t     # brownian time scale (invariant)

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

# Compute parameters from activities
if peA != 0:                        # A particles are NOT Brownian
    vA = computeVel(peA)
    FpA = computeActiveForce(vA)
    epsA = computeEps(FpA)
    tauA = computeTauLJ(epsA)
else:                               # A particles are Brownian
    vA = 0.0
    FpA = 0.0
    epsA = kT
    tauA = computeTauLJ(epsA)

if peB != 0:                        # B particles are NOT Brownian
    vB = computeVel(peB)
    FpB = computeActiveForce(vB)
    epsB = computeEps(FpB)
    tauB = computeTauLJ(epsB)
else:                               # B particles are Brownian
    vB = 0.0
    FpB = 0.0
    epsB = kT
    tauB = computeTauLJ(epsB)

# Import data into arrays
tsteps,\
gasA,\
gasB,\
gasAll,\
denseA,\
denseB,\
denseAll,\
lgClust,\
dpDensity,\
gpDensity,\
dpArea,\
gpArea = np.loadtxt(inTxt, skiprows=1, unpack=True)

# Get adjusted time units
tauLJ = (tauA if (tauA <= tauB) else tauB)  # use the smaller tauLJ
ratio = tauLJ / tauBrown                    # get rid of LJ units
dt = ratio * 0.00001                        # now you have size of a timestep

# Normalize data (terms of tauBrown and percent)
partNum = gasAll[0]
tsteps *= dt
gasA /= partNum * partFracA
gasB /= partNum * (1-partFracA)
gasAll /= partNum
denseA /= partNum * partFracA
denseB /= partNum * (1-partFracA)
denseAll /= partNum
lgClust /= partNum

# Now you have your data, plot it starting with...

# Gas phase
plt.plot(tsteps, gasA, label='Slow-type')
plt.plot(tsteps, gasB, label='Fast-type')
plt.plot(tsteps, gasAll, label='All')
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Fraction in Gas Phase')
plt.ylim(0,1)
plt.savefig('gasPhase_' + out + '.png', dpi=1000)
plt.close()

# Dense phase
plt.plot(tsteps, denseA, label='Slow-type')
plt.plot(tsteps, denseB, label='Fast-type')
plt.plot(tsteps, denseAll, label='All')
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Fraction in Dense Phase')
plt.ylim(0,1)
plt.savefig('densePhase_' + out + '.png', dpi=1000)
plt.close()

# Largest Cluster
plt.plot(tsteps, lgClust)
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Fraction in Largest Cluster')
plt.ylim(0,1)
plt.savefig('lgCluster_' + out + '.png', dpi=1000)
plt.close()

# Densities
plt.plot(tsteps, dpDensity, label='Dense Phase')
plt.plot(tsteps, gpDensity, label='Gas Phase')
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Density')
plt.savefig('phaseDensity_' + out + '.png', dpi=1000)
plt.close()

# Areas
plt.plot(tsteps, dpArea, label='Dense Phase')
plt.plot(tsteps, gpArea, label='Gas Phase')
plt.legend()
plt.xlabel(r'Brownian Time $(\tau_{Brownian})$')
plt.ylabel('Area')
plt.savefig('phaseArea_' + out + '.png', dpi=1000)
plt.close()
