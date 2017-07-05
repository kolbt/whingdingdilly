import sys
import matplotlib
#matplotlib.use('Agg')
matplotlib.use('macosx')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# start with 1D array
file = str(sys.argv[1])
data = np.loadtxt(file, dtype=np.int8)

phase_3D = np.zeros((11,16,16), dtype=np.int8)

for j in range(0,16):
    # change this array so that it is 11(rows) x 16(columns)
    #phase_diag = np.zeros((11,16), dtype=np.int8)

    # for loop the fuck outta this
    r=0
    c=0
    for i in range(0,len(data)):
        phase_3D[r][c][j] = data[i][j]
        r += 1
        if r == 11:
            r = 0
            c += 1

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for idx in range(0,11):
    for idy in range(0,16):
        for idz in range(0,16):
            if phase_3D[idx,idy,idz] != 0:
                ax.scatter(idx, idy, idz, c='b')
#            else:
#                ax.scatter(idx, idy, idz, c='w')

ax.set_xlabel('Particle Fraction of A')
ax.set_xlim3d(0,10)
ax.set_ylabel('Pe_A')
ax.set_ylim3d(0,15)
ax.set_zlabel('Pe_B')
ax.set_zlim3d(0,15)
plt.show()

# WORD. Now plot that shit
#plt.plot(phase_3D, cmap='Blues')
fig.savefig('phase_3D.png', dpi=1000)
