import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
from scipy.misc import imread

# Import the PeB plane value
inPeB = float(sys.argv[1])
print(inPeB)

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
#kappa = 3.8
kappa = 4.05
phi_min = 0.6

def solvePartFrac(PeA, PeB):
    xA = ((3 * (np.pi**2) * kappa) - (4 * phi_min * PeB)) / ((4 * phi_min) * (PeA - PeB))
    if xA <= 0.0:
        xA = 5
    return xA

x = np.arange(0.0, 160.0, 0.001)
y = np.zeros_like(x)

fig = plt.figure()
ax = fig.add_subplot(111)

for j in xrange(len(y)):
    y[j] = solvePartFrac(x[j], inPeB)
ax.plot(x, y, c='k', linestyle='--', zorder=1)

# Plot the phase plane as background
img = imread('pb' + str(int(inPeB)) + '.png')
ax.imshow(img, zorder=0, extent=[-5.0, 155, -0.05, 1.05])
ax.set_title(r'$Pe_{B}=$' + str(int(inPeB)), fontsize=30, y=1.02)
ax.set_xlabel(r'$Pe_{A}$')
ax.set_ylabel(r'$x_{A}$')
ax.set_xlim((-5, 155.0))
ax.set_ylim((-0.05, 1.05))

plt.axes().xaxis.set_minor_locator(xminor)
plt.axes().yaxis.set_minor_locator(yminor)
plt.axes().xaxis.set_major_locator(xmajor)
plt.axes().yaxis.set_major_locator(ymajor)

# This should fix aspect ratio issues
ratio = (float(img.shape[0])/img.shape[1])
print(ratio)
xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)

#plt.tight_layout()
plt.savefig('HS_peB_' + str(int(inPeB)) + '.png' ,dpi=1000)
plt.close()
