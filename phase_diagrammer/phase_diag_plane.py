import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

## start with 1D array
#file = str(sys.argv[1])
#pb = float(sys.argv[2])
#data = np.loadtxt(file, dtype=np.int8)
#
## change this array so that it is 11(rows) x 16(columns)
#phase_diag = np.zeros((11,16), dtype=np.int8)
#
## for loop the fuck outta this
#r=10
#c=0
##for i in range(0,len(data)):
##    phase_diag[r][c] = data[i]
##    r += 1
##    if r == 11:
##        r = 0
##        c += 1
#
#for i in range(0,len(data)):
#    phase_diag[r][c] = data[i]
#    r -= 1
#    if r == -1:
#        r = 10
#        c += 1

# WORD. Now plot that shit
# Instead of using imshow, let's make a scatterplot and overlay theory
kappa = 1.2
phi = 0.6

def theory(xA, PeA, PeB):
    out = (3 * (np.pi**2) * kappa) / (4 * phi * (xA * (PeA - PeB) + PeB))
    return out

def solvePartFrac(PeA, PeB):
    xA = ((3 * (np.pi**2) * kappa) / (4 * phi)) - (PeB / (PeA - PeB))
    return xA

x = np.arange(0, 160, 10)
y = np.zeros_like(x)
for iii in range(0, 15):
    y[iii]=solvePartFrac(x[iii], 150)


plt.plot(x, y)
plt.show()

