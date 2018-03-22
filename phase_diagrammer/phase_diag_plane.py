import sys
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Your files are on Hagrid: /Volumes/Hagrid/clust_1000/pe#B/*/pe#B_sep.txt
    
# start with 1D array
file = str(sys.argv[1])
PeB = int(sys.argv[2])
data = np.loadtxt(file, dtype=np.int8)

# change this array so that it is 11(rows) x 16(columns)
phase_diag = np.zeros((11,16), dtype=np.int8)

# for loop the fuck outta this
r=0
c=0

for i in range(0,len(data)):
    phase_diag[r][c] = data[i]
    r += 1
    if r == 11:
        r = 0
        c += 1

# Plot the experimental data
for lll in range(0, 11):
    for mmm in range(0,16):
        xA = float(lll) / 10.0
        PeA = float(mmm) * 10.0
        if phase_diag[lll][mmm] == 1:
            plt.scatter(PeA, xA, c='k')
        else:
            plt.scatter(PeA, xA, facecolor='none', edgecolor='k')

# Instead of using imshow, let's make a scatterplot and overlay theory
#kappa = 2.22
#phi_min = 0.40

kappa = 3.8
phi_min = 0.6

def solvePartFrac(PeA, PeB):
    xA = ((3 * (np.pi**2) * kappa) - (4 * phi_min * PeB)) / ((4 * phi_min) * (PeA - PeB))
    return xA

x = np.arange(0.0, 150.0, 0.001)
y = np.zeros_like(x)

for iii in range(0, len(x)):
    y[iii] = solvePartFrac(x[iii], PeB)
plt.plot(x, y, 'k')
plt.xlim([0.0, 150.0])
plt.ylim([0.0, 1.0])
plt.xlabel('Activity A')
#plt.xlabel(r'Time ($\tau$)')
plt.ylabel('Particle Fraction')
plt.savefig('phase_diagram_pb' + str(PeB) + '.png', dpi=1000)
plt.close()

def findSpinodal()


