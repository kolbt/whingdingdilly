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

def computeBradyMono(PHI, PE, APART, PHICP):
    '''Computes pressure, straight from the Brady paper'''
    ndense = APART / PHI
    return ndense * ((PE**2)/6.0) * ( 1 - PHI - PHI**2 + PHI*(9.0/(2.0 * PE)) * (PHICP/(1-PHI)) )

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
PeRs2D = np.zeros_like(inPhis)
PeRs3D = np.zeros_like(inPhis)
Pes2D = np.zeros_like(inPhis)
Pes3D = np.zeros_like(inPhis)
#phi02D = np.pi / (2.0 * np.sqrt(3.0))
phi02D = 0.9 # this is what Brady uses
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

# Plotting individually:
#plt.semilogy(inPhis, PeRs2D)
#plt.xlim(0.2, 0.9)
#plt.ylim(2*10**-3, 4*10**-1)
#plt.xlabel(r'$\phi$')
#plt.ylabel(r'$Pe_{R}$')
#plt.title(r'2D Spinodal')
#plt.show()
#
#plt.plot(inPhis, Pes2D)
#plt.xlim(0.2, 0.9)
#plt.ylim(0, 500)
#plt.xlabel(r'$\phi$')
#plt.ylabel(r'$Pe$')
#plt.title(r'2D Spinodal')
#plt.show()
#
#plt.plot(inPhis, Pes3D)
#plt.xlim(0.2, 0.7)
#plt.ylim(0, 500)
#plt.xlabel(r'$\phi$')
#plt.ylabel(r'$Pe$')
#plt.title(r'3D Spinodal')
#plt.show()

## Make a figure
#fig, ax = plt.subplots(1, 4, figsize=(14,4))
#ax[0].semilogy(inPhis, PeRs2D)
#ax[1].semilogy(inPhis, PeRs3D)
#ax[2].plot(inPhis, Pes2D)
#ax[3].plot(inPhis, Pes3D)
## Limits
#ax[0].set_xlim(0.2, 0.9)
#ax[1].set_xlim(0.2, 0.9)
#ax[2].set_xlim(0.2, 0.9)
#ax[3].set_xlim(0.2, 0.9)
#ax[0].set_ylim(10**-3, 10**-1)
#ax[1].set_ylim(10**-3, 10**-1)
#ax[2].set_ylim(0, 500)
#ax[3].set_ylim(0, 500)
## Labels
#ax[0].set_title(r'2D')
#ax[1].set_title(r'3D')
#ax[2].set_title(r'2D')
#ax[3].set_title(r'3D')
#ax[0].set_ylabel(r'$Pe_{R}$')
#ax[2].set_ylabel(r'$Pe$')
#ax[0].set_xlabel(r'$\phi$')
#ax[1].set_xlabel(r'$\phi$')
#ax[2].set_xlabel(r'$\phi$')
#ax[3].set_xlabel(r'$\phi$')
#ax[0].text(0.75, 1.1, 'Spinodal (vs. $Pe_{R}$)', size=16, transform=ax[0].transAxes)
#ax[2].text(0.75, 1.1, 'Spinodal (vs. $Pe$)', size=16, transform=ax[2].transAxes)
## Ticks
#ax[1].set_yticklabels([])
#ax[3].set_yticklabels([])
#for i in xrange(4):
#    ax[i].tick_params(direction='in', which='both')
#plt.tight_layout(pad=2.0, w_pad=0.5, h_pad=1.0)
#plt.savefig('Brady_spinodals.png', dpi=1000)
#plt.close()




# We need to do this with symbolic python
import sympy as sp

# 2D Parameters
beta = 4.0 / np.pi
xi = 0.2
phi0 = 0.90
# 3D Parameters
beta = 3.0
xi = 1.0
phi0 = 0.64

PHI, PER = sp.symbols("PHI PER")
# Each term gets its own line
sp.plot_implicit(sp.Eq(
                       ((PHI / (1.0 - PHI)) *
                        ((1.0/PHI) +
                         (((beta * PER * phi0) / phi0) * ((1.0-(PHI/phi0))**-1)) +
                         (3.0*xi*(PHI**2)) -
                         (PHI) +
                         (((beta * PER * phi0) / phi0) * (1.0 - phi0) * ((1.0-(PHI/phi0))**-2)) -
                         (3.0) +
                         ((beta * PER * phi0)))) -
                       (1.0) +
                       (2.0*PHI) +
                       (xi * (PHI**2)) -
                       (((beta * PER * phi0) / phi0) *
                        ((2.0*PHI)*((1.0-(PHI/phi0))**-1)) +
                        (((PHI**2)/phi0)*((1.0-(PHI/phi0))**-2)))),
                 (PHI, 0.2, 0.7),
                 (PER, 10**-3, 10**1),
                 yscale='log')

# Modify this plot with the following link
# https://stackoverflow.com/questions/40747474/sympy-and-plotting
