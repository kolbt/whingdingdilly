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
kappa = 3.6475626
#kappa = 4.05
phi_min = 0.6

def solvePartFrac(PeA, PeB):
    xA = ((3 * (np.pi**2) * kappa) - (4 * phi_min * PeB)) / ((4 * phi_min) * (PeA - PeB))
    if xA <= 0.0:
        xA = 5
    return xA
    
def solvePeNetCrit(PeA, PeB, PeNCrit=45.):
    '''Pe_net = (PeA * xA) + (PeB - (PeB * xA))'''
    xA = (PeNCrit - PeB) / (PeA - PeB)
    return xA

x = np.arange(0.0, 160.0, 0.001)
y = np.zeros_like(x)
y2 = np.zeros_like(x)

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
        y2[j] = solvePeNetCrit(x[j], distinctPePlanes[i])
    plt.plot(x, y, c='k', linestyle='--')
#    plt.plot(x, y2, c='b', linestyle='--')
    plt.title(r'$Pe_{B}=$' + str(distinctPePlanes[i]), fontsize=30, y=1.02)
    plt.xlabel(r'$Pe_{A}$')
    plt.ylabel('$x_{A}$')
    plt.xlim((-5, 155.0))
    plt.ylim((-0.05, 1.05))

    plt.axes().xaxis.set_minor_locator(xminor)
    plt.axes().yaxis.set_minor_locator(yminor)
    plt.axes().xaxis.set_major_locator(xmajor)
    plt.axes().yaxis.set_major_locator(ymajor)

    ax = plt.gca()
    ratio = 0.685
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)

    plt.tight_layout()
    plt.savefig('HS_peB_' + str(distinctPePlanes[i]) + '.png' ,dpi=1000)
#    plt.savefig('HS_peB_0' + str(count) + '.png' ,dpi=1000)
    plt.close()

    count += 1


