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

# File to pull data from
inTxt = "effective_particle_radius.txt"

# Base name to write to
out = "pa"+str(peA)+\
"_pb"+str(peB)+\
"_xa"+str(partPercA)

# Import data into arrays
peA,\
peB,\
peR,\
xA,\
eps,\
modeALL,\
modeAA,\
modeAB,\
modeBB,\
fALL,\
fAA,\
fAB,\
fBB,\
phiEff = np.loadtxt(inTxt, skiprows=1, unpack=True)

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
