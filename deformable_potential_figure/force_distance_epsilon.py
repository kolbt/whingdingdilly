'''
We want to know the relationship between Pe, epsilon and sigma...
'''

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
import pandas as pd
import seaborn as sns

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

# This gives me the relationship for a head-on collision: pe vs eps vs sigma
peRange = np.arange(0., 500., 1.)
epsRange = np.arange(0.01, 1., 0.01)
x, y, zi, zii, mati, matii = convergeR(peRange, epsRange, angle)

fig = plt.figure(figsize=(10, 8))
ax = []
# 3D Plot for: Pe, epsilon, sigma
print("Plotting 3D surfaces...")
ax.append(fig.add_subplot(221, projection='3d'))
ax[0].plot_trisurf(x, y, zi)
ax[0].set_xlabel(r'Activity (Pe)')
ax[0].set_ylabel(r'Repulsive strength $(\epsilon)$')
ax[0].set_zlabel(r'Equilibrium separation $(\sigma_{eff})$')
ax[0].zaxis._axinfo['juggled'] = (1,2,0)
    
# 3D Plot for: Pe, epsilon, phi
ax.append(fig.add_subplot(222, projection='3d'))
ax[1].plot_trisurf(x, y, zii)
#ax[1].view_init(25, 325)
ax[1].set_xlabel(r'Activity (Pe)')
ax[1].set_ylabel(r'Repulsive strength $(\epsilon)$')
ax[1].set_zlabel(r'Liquid phase area fraction $(\phi_{l})$')
ax[1].zaxis._axinfo['juggled'] = (1,2,0)

# Heatmap for: Pe, epsilon, sigma
print("Plotting heatmaps...")
ax.append(fig.add_subplot(223))
imi = ax[2].imshow(mati, extent=[min(peRange), max(peRange), min(epsRange), max(epsRange)], aspect='auto', origin='lower')
ax[2].set_xlabel('Activity (Pe)')
ax[2].set_ylabel(r'Repulsive strength $(\epsilon)$')
cbari = fig.colorbar(imi, ax=ax[2])
cbari.ax.get_yaxis().labelpad = 15
cbari.ax.set_ylabel(r'$\sigma$', rotation=270)

# Heatmap for: Pe, epsilon, phi
ax.append(fig.add_subplot(224))
imii = ax[3].imshow(matii, extent=[min(peRange), max(peRange), min(epsRange), max(epsRange)], aspect='auto', origin='lower')
ax[3].set_xlabel('Activity (Pe)')
ax[3].set_ylabel(r'Repulsive strength $(\epsilon)$')
cbarii = fig.colorbar(imii, ax=ax[3])
cbarii.ax.get_yaxis().labelpad = 15
cbarii.ax.set_ylabel(r'$\phi_{l}$', rotation=270)

# Now save the figure as an image
fig.tight_layout(w_pad=1.0)
plt.savefig("activity_epsilon_r_phi.png", dpi=1000, bbox_inches='tight', pad_inches=0.2)
plt.close()
