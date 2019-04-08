'''
Let's plot Brady's theory for our variables...
'''

import sys
import os
import numpy as np
import math
import matplotlib.pyplot as plt

def computeActivity(VP, TAUR, SIG):
    return (3.0 * VP * TAUR / SIG)

tau_R = 1.0 / 3.0
sigma = 1.0
particle_Area = np.pi * ((sigma / 2.0)**2)
# 2D
#phi_CP = np.pi / (2.0 * np.sqrt(3.0))

# 3D
phi_CP = np.pi / (3.0 * np.sqrt(2.0))

# Random close packing in 3D
phi_CP = 0.64

def computeBradyMono(PHI, PE, APART, PHICP):
    '''Computes pressure, straight from the Brady paper'''
    ndense = APART / PHI
    return ndense * ((PE**2)/6.0) * ( 1 - PHI - PHI**2 + PHI*(9.0/(2.0 * PE)) * (PHICP/(1-PHI)) )

def computeMonoSpinodal(PHI, PHICP):
    top = (3 * (PHI**2)) + (2 * PHI) - 1
    term = 1.0 - (PHI/PHICP)
    bottom1 = 2.0 * PHI * ((term)**-1)
    bottom2 = ((PHI**2)/PHICP) * ((term)**-2)
    return top / (3 * (bottom1 + bottom2))

inPhis = np.arange(0.2, phi_CP, 0.001)
PeRs = np.zeros_like(inPhis)
Pes = np.zeros_like(inPhis)
for i in xrange(len(inPhis)):
    PeRs[i] = computeMonoSpinodal(inPhis[i], phi_CP)
    Pes[i] = (PeRs[i])**-1 * (3.0/2.0)

plt.semilogy(inPhis, PeRs)
plt.ylim(10**-3, 10**-1)
plt.xlabel(r'$\phi$')
plt.ylabel(r'$Pe_{R}$')
plt.show()

plt.plot(inPhis, Pes)
plt.xlim(0.336, 0.7)
plt.ylim(0, 500)
plt.xlabel(r'$\phi$')
plt.ylabel(r'$Pe$')
plt.show()
