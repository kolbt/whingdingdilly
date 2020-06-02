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

    

#def phig(phi, pe, a):
#    sig = 1.
#    kap = 4.5
#    num = kap * (np.pi**2) * (a**2)
#    den = 4. * pe * (sig**2)
#    return num / den
    
#def phig(phi, pe, a):
#    sig = 1.
#    kap = 4.5
#    num = kap * (np.pi**2)
#    den = 4. * pe * (a**1)
#    return num / den
        
# Let's create a mesh and evaluate this
pes = np.arange(35., 500., 1.)
phi = 0.65
#aa = np.arange(0.1, 10., 0.1)

eps = [0.0001, 0.001, 0.01, 0.1, 1.]
for i in eps:
    aa = []
    fcs = []
    phigs = []
    outCF = []
    for j in pes:
        aa.append(conForRClust(j, i))
        instPhiG = compPhiG(j, aa[-1])
        outCF.append(clustFrac(phi, instPhiG, aa[-1]))
#        fcs.append(fc(phi, j, aa[-1]))
#        phigs.append(phig(0.6, j, aa[-1]))
#    plt.plot(pes, fcs, label=i, c=plt.cm.jet(eps.index(i)/float(len(eps))), ls='--')
#    plt.plot(pes, phigs, label=i, c=plt.cm.jet(eps.index(i)/float(len(eps))))
    plt.plot(pes, outCF, label=i, c=plt.cm.jet(eps.index(i)/(float(len(eps))-1) ))
    
ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
for i in range(0, len(labels)):
    for j in range(0, len(labels)):
        if labels[j] < labels[i] and j > i:
            labels[i], labels[j] = labels[j], labels[i]
            handles[i], handles[j] = handles[j], handles[i]
by_label = OrderedDict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), title=r'Softness $(\epsilon)$')
    
plt.xlim(0, 500)
plt.ylim(0., 1.)
plt.xlabel('Activity')
plt.ylabel('Cluster fraction')
plt.legend()
plt.savefig('soft_kinetic_theory.pdf', dpi=1000)
plt.close()
