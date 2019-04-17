import sys
import matplotlib
#matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd

## Here are my rc parameters for matplotlib
mpl.rc('font', serif='Helvetica Neue')
mpl.rcParams.update({'font.size': 18})
#mpl.rcParams['figure.figsize'] = 3.2, 2.8
mpl.rcParams['figure.dpi'] = 1000
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['lines.linewidth'] = 2.0

mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['xtick.minor.width'] = 1.5
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['ytick.minor.width'] = 1.5

xminor = MultipleLocator(10)
xmajor = MultipleLocator(50)
yminor = MultipleLocator(0.1)
ymajor = MultipleLocator(0.5)

# Theory fits 99.3% of the data
kappa = 3.8
#kappa = 4.05
phi_min = 0.6

def solvePartFrac(PeA, PeB):
    xA = ((3 * (np.pi**2) * kappa) - (4 * phi_min * PeB)) / ((4 * phi_min) * (PeA - PeB))
    if xA <= 0.0:
        xA = 5
    return xA

phi0=0.90
xi=0.2
beta=4.0/np.pi
def computeOmega(PE, PHI):
    if PE == 0:
        return 0
    '''This will make the spinodal function a bit cleaner'''
    omega = (((3.0/2.0)*(PE))**2) *\
        (1-(2*PHI)-(3*xi*(PHI**2)) +\
         ((beta*((3.0/2.0)*(PE**-1))) *\
          ((2*PHI)*((1-(PHI/phi0))**-1) +\
           (PHI/phi0)*((1-(PHI/phi0))**-2))))
    return omega

def bradySpinodal(PeA, PeB, phi):
    '''
        This method uses the spinodal derived from setting the entire expression
        for the derivative of the active pressure to 0 (not each component).
    '''
    try:
        xA = computeOmega(PeB, phi) / (computeOmega(PeB, phi)-computeOmega(PeA, phi))
    except ZeroDivisionError:
        xA = None
    if xA > 1.0 or xA < 0.0:
        xA = None
    return xA

x = np.arange(0.0, 160.0, 0.001)
y = np.zeros_like(x)
spin = np.zeros_like(x)

file = 'hard_sphere_phase_separation.csv'
ps_data = pd.read_csv(file)
ps_data['xA'] /= 100.0

distinctPePlanes = []
for i in xrange(len(ps_data['PeB'])):
    if ps_data['PeB'][i] not in distinctPePlanes:
        distinctPePlanes.append(ps_data['PeB'][i])

count = 1
for i in xrange(len(distinctPePlanes)):
    for j in xrange(len(ps_data['PeB'])):
        # If data belongs to this plane, plot it
        if ps_data['PeB'][j] == distinctPePlanes[i]:
            if ps_data['PS'][j] == 1:
                plt.scatter(ps_data['PeA'][j], ps_data['xA'][j], c='k')
            else:
                plt.scatter(ps_data['PeA'][j], ps_data['xA'][j], facecolor = 'none', edgecolor = 'k')
    for j in xrange(len(y)):
        y[j] = solvePartFrac(x[j], distinctPePlanes[i])
        spin[j] = bradySpinodal(float(x[j]), float(distinctPePlanes[i]), 0.6)
    plt.plot(x, y, c='k')
    plt.plot(x, spin, c='b')
    plt.title(r'$Pe_{B}=$' + str(distinctPePlanes[i]), fontsize=30, y=1.02)
    plt.xlabel(r'$Pe_{A}$')
    plt.ylabel('$x_{A}$')
    plt.xlim((-5, 155.0))
    plt.ylim((-0.05, 1.05))

    plt.axes().xaxis.set_minor_locator(xminor)
    plt.axes().yaxis.set_minor_locator(yminor)
    plt.axes().xaxis.set_major_locator(xmajor)
    plt.axes().yaxis.set_major_locator(ymajor)

    plt.tight_layout()
    plt.savefig('brady_HS_peB_' + str(distinctPePlanes[i]) + '.png' ,dpi=1000)
    plt.savefig('brady_HS_peB_' + str(count) + '.png' ,dpi=1000)
    plt.close()

    count += 1


