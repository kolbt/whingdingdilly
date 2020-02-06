'''
We want to plot different coexistance curves (pe vs. phi) for a few different epsilon
'''

import numpy as np
import matplotlib.pyplot as plt

# Functions
angle = np.pi
# Lennard-Jones potential
def ljPotential(r, eps, sigma=1.):
    div = (sigma/r)
    U = ( 4. * eps * ((div)**12 - (div)**6) ) + eps
    return U
# Lennard-Jones force
def ljForce(r, eps, sigma=1.):
    div = (sigma/r)
    dU = (24. * eps / r) * ((2*(div**12)) - (div)**6)
    return dU
# Force of a collision at some angle and activity
def collisionForce(pe, angle):
    return pe - (pe * np.cos(angle))
# From area fraction, get lattice spacing
def phiToLat(phiIn):
    '''Read in phi, output the lattice spacing'''
    phiCP = np.pi / (2. * np.sqrt(3.))
    latCP = 1.
    return np.sqrt(phiCP / phiIn)
# From lattice spacing, get area fraction
def latToPhi(latIn):
    '''Read in lattice spacing, output phi'''
    phiCP = np.pi / (2. * np.sqrt(3.))
    latCP = 1.
    return phiCP / (latIn**2)
# Get effective separation/density from a particular collision type
def convergeR(pe, eps, angle):
    x = []
    y = []
    zi = []
    zii = []
    mati = []
    matii = []
    count = 0
    for i in pe:
        mati.append([])
        matii.append([])
        for j in eps:
            r = 1.
            # Is the LJ-Force greater than the collision force?
            while ljForce(r, j) < collisionForce(i, angle):
                # Decrement distance
                r -= 0.001
            x.append(i)
            y.append(j)
            zi.append(r)
            zii.append(latToPhi(r))
            mati[count].append(r)
            matii[count].append(zii[-1])
        count += 1
    return x, y, zi, zii, np.asarray(mati).T, np.asarray(matii).T
    
# Let's do this at constant epsilon
def convergeConstEps(pe, eps, angle):
    out = []
    for i in pe:
        r = 1.
        while ljForce(r, eps) < collisionForce(i, angle):
            r -= 0.0001
        out.append(latToPhi(r))
    return out

epsRange = np.arange(0.1, 1.1, 0.1)
peRange = np.arange(0., 500., 1.)

# Get phi vs pe at constant epsilon
for i in epsRange:
    plt.plot(peRange, convergeConstEps(peRange, i, angle), c=plt.cm.jet(i/max(epsRange)), label="{0:.1f}".format(i))
plt.legend(title=r'$\epsilon$')
plt.xlim(0, 500)
plt.xlabel(r'Activity $(Pe)$')
plt.ylabel(r'Liquid phase area fraction $(\phi_{l})$')
plt.savefig('dense_phase_coexistence.png', dpi=1000, bbox_inches='tight', pad_inches=0)
plt.close()
