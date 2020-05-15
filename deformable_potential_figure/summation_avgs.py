'''
Compute the summation of interparticle distance to get average,
how does this overlay with the force average?
'''

import numpy as np
import matplotlib.pyplot as plt
import math

def ljPotential(r, eps, sigma=1.):
    div = (sigma/r)
    U = ( 4. * eps * ((div)**12 - (div)**6) ) + eps
    return U

def ljForce(r, eps, sigma=1.):
    div = (sigma/r)
    dU = (24. * eps / r) * ((2*(div**12)) - (div)**6)
    return dU
    
def convergeConstPeEps(pe, eps, theta):
    r = 2.**(1./6.)
    skip = [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001]
    for j in skip:
        while ljForce(r, eps) < (pe * np.cos(theta)):
            r -= j
        r += j
    return r
    
def convergeAvgF(eps, favg):
    r = 2.**(1./6.)
    skip = [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001]
    for j in skip:
        while ljForce(r, eps) < favg:
            r -= j
        r += j
    return r

# Softness to use
epsilons = [1., 0.1, 0.01, 0.001, 0.0001]
cols = [plt.cm.jet(0.0), plt.cm.jet(0.25), plt.cm.jet(0.5), plt.cm.jet(0.75), plt.cm.jet(1.0)]
# Width of my rieman sum
width = 0.01
# Limits of my integration
#a = -np.pi / 2.
a = 0.
b = np.pi / 2.
#print(np.cos(a))
#print(np.cos(b))
# Angles to integrate over
thetas = np.arange(a, b, width)

# Do this for range of pes
pes = np.arange(0., 500., 1.)

for k in range(0, len(epsilons)):
    # Lists to hold values
    listAvgR = []
    listOutR = []
    for m in range(0, len(pes)):
        # Hold the avg distance and force
        avgR = 0.
        avgF = 0.
        for i in range(0, len(thetas)):
            # Interparticle separation for this theta
            instR = convergeConstPeEps(pes[m], epsilons[k], thetas[i])
            avgR += (instR)
        #    avgR += (instR * width)
            # Pair force for this theta (must be > 0)
            instF = pes[m] * np.cos(thetas[i])
            if instF >= 0.:
                avgF += instF
            else:
                print(pes[m])
                print(thetas[i])
                print(np.cos(thetas[i]))

        # Compute the average value of the interparticle distance
        avgR /= len(thetas)
        listAvgR.append(avgR)
        #avgR /= np.pi

        # Now compute the interparticle separation from the average force
        avgF /= len(thetas)
        outR = convergeAvgF(epsilons[k], avgF)
        listOutR.append(outR)

    plt.plot(pes+50., listAvgR, c=cols[k], ls='-', label='Average distance')
    plt.plot(pes+50., listOutR, c=cols[k], ls='--', label='Average force')
    
from matplotlib.lines import Line2D
leg = [Line2D([0], [0], lw=1.5, c='k', ls='-', markeredgecolor='none', label=r'Average distance', markerfacecolor='none', markersize=10),
       Line2D([0], [0], lw=1.5, c='k', ls='--', markeredgecolor='none', label=r'Average force', markerfacecolor='none', markersize=10)]

plt.legend(handles=leg, loc='best', frameon=False)
plt.xlim(0., 500.)
plt.ylim(0.4, 2.**(1./6.))
plt.xlabel('Activity')
plt.ylabel('Interparticle distance')
plt.savefig('averaging_test.pdf', dpi=1000)
