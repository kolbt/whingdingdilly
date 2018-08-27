import sys
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
# 3D axes
from mpl_toolkits.mplot3d import Axes3D

kappa = 3.8
phi_min = 0.6

def solvePartFrac(PeA, PeB):
    xA = ((3 * (np.pi**2) * kappa) - (4 * phi_min * PeB)) / ((4 * phi_min) * (PeA - PeB))
    return xA

x = np.arange(0.0, 150.0, 0.001)
y = np.zeros_like(x)

# Plot the experimental data
# see file format in write-phase-txt.py
phase_file = "phase-separation-data.txt"

peA, \
peB, \
xA, \
ps = np.loadtxt(phase_file, skiprows=1, unpack=True)

xA /= 100.0

# Make a new array that has symmetry incorporated
peA_all = np.zeros((len(peA) * 2), dtype=np.float32)
peB_all = np.zeros((len(peB) * 2), dtype=np.float32)
xA_all = np.zeros((len(xA) * 2), dtype=np.float32)
ps_all = np.zeros((len(ps) * 2), dtype=np.int)
# Make list of distinct PeB values
distinctPeBs = []

sym = len(ps)
for i in range(0, len(ps)):
    # Copy in the original point
    peA_all[i] = peA[i]
    peB_all[i] = peB[i]
    xA_all[i] = xA[i]
    ps_all[i] = ps[i]
    # Take advantage of symmetry
    peA_all[i + sym] = peB[i]
    peB_all[i + sym] = peA[i]
    xA_all[i + sym] = 1.0 - xA[i]
    ps_all[i + sym] = ps[i]
    # Get all distinct PeB values
    if peB[i] not in distinctPeBs:
        distinctPeBs.append(peB[i])

for h in range(0, len(distinctPeBs)):
    for i in range(0, len(ps_all)):
        if peB_all[i] == distinctPeBs[h]:
            if ps_all[i] == 1:
                plt.scatter(peA_all[i], xA_all[i], c='k')
            else:
                plt.scatter(peA_all[i], xA_all[i], facecolor = 'none', edgecolor = 'k')
    # Overlay theory on raw data
    for j in range(0, len(x)):
        y[j] = solvePartFrac(x[j], distinctPeBs[h])
    plt.plot(x, y, 'k')
    # Axes labels, limits, etc.
    plt.title("PeB={}".format(distinctPeBs[h]))
    plt.xlim(0, 150)
    plt.ylim(0, 1.0)
    plt.savefig('plane' + str(int(distinctPeBs[h] / 10)) + '.png', dpi=1000)
    plt.close()

