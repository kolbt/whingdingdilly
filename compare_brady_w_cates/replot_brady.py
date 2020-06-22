'''
Let's plot Brady's theory for our variables...
'''

import sys
import os
import numpy as np
import math
import matplotlib.pyplot as plt

tau_R = 1.0 / 3.0
sigma = 1.0
particle_Area = np.pi * ((sigma / 2.0)**2)

def compute2DMonoSpinodal(PHI, PHICP):
    top = (3 * 0.2 * (PHI**2)) + (2 * PHI) - 1
    term = 1.0 - (PHI/PHICP)
    bottom1 = 2.0 * PHI * ((term)**-1)
    bottom2 = ((PHI**2)/PHICP) * ((term)**-2)
    return top / ((4.0 / np.pi) * (bottom1 + bottom2))

def compute3DMonoSpinodal(PHI, PHICP):
    top = (3 * (PHI**2)) + (2 * PHI) - 1
    term = 1.0 - (PHI/PHICP)
    bottom1 = 2.0 * PHI * ((term)**-1)
    bottom2 = ((PHI**2)/PHICP) * ((term)**-2)
    return top / (3.0 * (bottom1 + bottom2))

inPhis = np.arange(0.2, 0.9, 0.001)
#inPhis /= particle_Area
PeRs2D = np.zeros_like(inPhis)
PeRs3D = np.zeros_like(inPhis)
Pes2D = np.zeros_like(inPhis)
Pes3D = np.zeros_like(inPhis)
phi02D = 0.9 # this is what Brady uses
#phi02D /= particle_Area
phi03D = 0.64 # this is what Brady uses
for i in xrange(len(inPhis)):
    PeRs2D[i] = compute2DMonoSpinodal(inPhis[i], phi02D)
    PeRs3D[i] = compute3DMonoSpinodal(inPhis[i], phi03D)
    # This fixes issues with plotting against conventional Pe
    if inPhis[i] <= 0.444:
        Pes2D[i] = 1000
    else:
        Pes2D[i] = (PeRs2D[i])**-1 * (3.0/2.0)

    if inPhis[i] <= 0.336 or inPhis[i] >= 0.6:
        Pes3D[i] = 1000
    else:
        Pes3D[i] = (PeRs3D[i])**-1 * (3.0/2.0)
    # Fix PeR weirdness
    if inPhis[i] >= 0.62:
        PeRs3D[i] = 10**-4
        
print((0.028)**-1 * (3.0/2.0))
        
# Let's plot the 2D expression
inPhis /= particle_Area
plt.plot(Pes2D, inPhis)
plt.xlim(0, 120)
plt.ylim(0.1, 1.3)
plt.xlabel(r'Activity ($Pe$)')
plt.ylabel(r'$\rho/\rho_{cp}$')
plt.show()


