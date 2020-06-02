'''
A schematic figure to illustrate how Pe (and F_act) sets deformability:
    -LJ potential
    -Overlay 2 Forces (Strong under week)
    -Corresponds to collision angle
'''

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
from PIL import Image
import matplotlib.patches as patches
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
from collections import OrderedDict
from matplotlib.lines import Line2D
fsize = 9
plt.rcParams.update({'font.size': fsize})
params = {'legend.fontsize': fsize,
          'axes.labelsize': fsize,
          'axes.titlesize': fsize,
          'xtick.labelsize': fsize,
          'ytick.labelsize': fsize}
plt.rcParams.update(params)

# Get lattice spacing for particle size
def ljForce(r, eps, sigma=1.):
    div = (sigma/r)
    dU = (24. * eps / r) * ((2*(div**12)) - (div)**6)
    return dU
    
def avgCollisionForce(pe, power=1.):
    '''Computed from the integral of possible angles'''
    magnitude = np.sqrt(28.)
    return (magnitude * (pe)) / (np.pi)
    
def conForRClust(pe, eps):
    out = []
    r = 1.112
    skip = [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001]
    for j in skip:
        while ljForce(r, eps) < avgCollisionForce(pe):
            r -= j
        r += j
    out = r
    return out
    
def latToPhi(latIn):
    '''Read in lattice spacing, output phi'''
    phiCP = np.pi / (2. * np.sqrt(3.))
    return phiCP / (latIn**2)
    
def fc(phi, pe, a):
    sig = 1.
    kap = 4.05
    # This should decrease with decreasing lattice spacing
    num = (4.*phi*pe*(a)) - (3.*(np.pi**2)*kap*(sig**2))
    # This should increase with increasing lattice spacing
    den = (4.*phi*pe*(a)) - ((latToPhi(a)**-1)*(3.*(np.pi**2)*kap*phi*(sig**2.)))
    return num / den
    
def compPhiG(pe, a, kap=4.05, sig=1.):
    num = 3. * (np.pi**2) * kap * sig
    den = 4. * pe * a
    return num / den
    
def clustFrac(phi, phiG, a, sig=1.):
    phiL = latToPhi(a)
    ApL = np.pi * (sig**2) / 4.
    Ap = np.pi * (sig**2) / 4.
    num = (phiL*phiG) - (phiL*phi)
    den = ((ApL/Ap)*phi*phiG) - (phi*phiL)
    ans = num / den
    return ans
    
# Lennard-Jones pressure
def ljPress(r, eps, sigma=1.):
    div = (sigma/r)
    dU = (24. * eps / r) * ((2.*(div**12.)) - (div)**6.)
    return dU * np.sqrt(28.) * r
        
# Let's create a mesh and evaluate this
pes = np.arange(35., 500., 1.)
phi = 0.65
eps = [0.0001, 0.001, 0.01, 0.1, 1.]

# Let's use analytical theory and kinetic theory to get cluster radius
epsRange = [0.0001, 0.001, 0.01, 0.1, 1.]
epsRange = [1., 0.0001]
peRange = np.arange(35., 500., 1.)
phi = 0.65
phiRange = np.arange(0.45, 0.75, 0.1)
phiLines = ['-', ':', '--']
N = 100000.

# Now we can plot the surface tension
fig, ax = plt.subplots(1, 3, figsize=(15, 5))

phiCP = np.pi / (2. * np.sqrt(3))
for m in range(0, len(phiRange)):
    lat = []
    pColl = []
    pLJ = []
    cfs = []
    Rls = []
    peCrit = []
    for i in range(0, len(epsRange)):
        lat.append([])
        pColl.append([])
        pLJ.append([])
        cfs.append([])
        Rls.append([])
        for j in range(0, len(peRange)):
            # Compute lattice spacing
            lat[i].append(conForRClust(peRange[j], epsRange[i]))
            
            # Compute pressure
            curPLJ = ljPress(lat[i][-1], epsRange[i]) / (np.pi * (lat[i][-1]**2) * 0.25 * phiCP)
            # Append to list
            pLJ[i].append(curPLJ/(10.**4))
            
            # Compute cluster fraction
            phiG = compPhiG(peRange[j], lat[i][-1])
            cf = clustFrac(phiRange[m], phiG, lat[i][-1])
            if cf < 0. or cf > 1.:
                cf = 0.
            cfs[i].append(cf)
            
            # Get the critical activity
            if j > 0:
                if cfs[i][-2] == 0. and cfs[i][-1] > 0.:
                    peCrit.append(peRange[j])
            
            # Get the radius (for some N)
            Nl = cfs[i][-1] * N
            Al = Nl * (np.pi * (lat[i][-1]**2) * 0.25)
            Rl = np.sqrt(Al / np.pi)
            Rls[i].append(Rl)

    # Plot the analytical component
    for i in range(0, len(epsRange)):
        st = [a * b for a, b in zip(pLJ[i], Rls[i])]
        divPe = [ a / (b - peCrit[i]) for a, b in zip(st, peRange)]
        ax[0].plot(peRange, st, c=plt.cm.jet(float(i)/(len(epsRange)-1)),
                 lw=1.5, zorder=0, label=epsRange[i], ls=phiLines[m])
        ax[1].plot(peRange - peCrit[i], st, c=plt.cm.jet(float(i)/(len(epsRange)-1)),
                 lw=1.5, zorder=0, label=epsRange[i], ls=phiLines[m])
        ax[2].plot(peRange - peCrit[i], divPe, c=plt.cm.jet(float(i)/(len(epsRange)-1)),
                 lw=1.5, zorder=0, label=epsRange[i], ls=phiLines[m])

# Make a simulation vs theory legend
phi_leg = []
for i in range(0, len(phiRange)):
    phi_leg.append(Line2D([0], [0], lw=1.5, ls=phiLines[i], c='k', markeredgecolor='none',
                          label="{:.2f}".format(phiRange[i]), markerfacecolor='none',
                          markersize=10))
add_leg = ax[2].legend(handles=phi_leg, frameon=False, handletextpad=0.1, title=r'$\phi$',
                    bbox_transform=ax[2].transAxes, bbox_to_anchor=[0.75, 1.0])
ax[2].add_artist(add_leg)


handles, labels = ax[2].get_legend_handles_labels()
for i in range(0, len(labels)):
    for j in range(0, len(labels)):
        if labels[j] < labels[i] and j > i:
            labels[i], labels[j] = labels[j], labels[i]
            handles[i], handles[j] = handles[j], handles[i]
by_label = OrderedDict(zip(labels, handles))
ax[2].legend(by_label.values(), by_label.keys(), title=r'Softness $(\epsilon)$', frameon=False)
    
# Set x limits
ax[0].set_xlim(0, 500)
ax[1].set_xlim(0, 500)
ax[2].set_xlim(0, 500)
# Set y limits
ax[0].set_ylim(0,)
ax[1].set_ylim(0,)
ax[2].set_ylim(0,)

# Set x labels
ax[0].set_xlabel(r'Activity $(Pe)$')
ax[1].set_xlabel(r'Activity $(Pe-Pe_{c})$')
ax[2].set_xlabel(r'Activity $(Pe-Pe_{c})$')

# Set y labels
ax[0].set_ylabel(r'$\gamma \ = \ \left[\Pi^{p}_{bulk} \cdot R_{dense}\right]/10^{4}$')
ax[1].set_ylabel(r'$\gamma$')
ax[2].set_ylabel(r'$\gamma / (Pe-Pe_{c}) $')

plt.tight_layout()
plt.savefig('eps_phi_analytical_effects.pdf', dpi=500)
plt.close()
        
