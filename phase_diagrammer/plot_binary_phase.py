import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# start with 1D array
file = str(sys.argv[1])
data = np.loadtxt(file, dtype=np.int8)

# change this array so that it is 11(rows) x 16(columns)
phase_diag = np.zeros((11,16), dtype=np.int8)

# for loop the fuck outta this
r=0
c=0
for i in range(0,len(data))
    phase_diag[r][c] = data[i]
    r += 1
    if r == 12:
        r = 0
        c += 1

# WORD. Now plot that shit
plt.imshow(phase_diag, cmap='hot', interpolation='nearest')
plt.show()


