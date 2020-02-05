'''
We want to know the relationship between Pe, epsilon and sigma...
'''

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d

# Functions
angle = np.pi
def ljPotential(r, eps, sigma=1.):
    div = (sigma/r)
    U = ( 4. * eps * ((div)**12 - (div)**6) ) + eps
    return U

def ljForce(r, eps, sigma=1.):
    div = (sigma/r)
    dU = (24. * eps / r) * ((2*(div**12)) - (div)**6)
    return dU
    
def collisionForce(pe, angle):
    return pe - (pe * np.cos(angle))
            
def convergeR(pe, eps, angle):
    x = []
    y = []
    z = []
    for i in pe:
        for j in eps:
            r = 1.
            # Is the LJ-Force greater than the collision force?
            while ljForce(r, j) < collisionForce(i, angle):
                # Decrement distance
                r -= 0.001
            x.append(i)
            y.append(j)
            z.append(r)
    return x, y, z

peRange = np.arange(0., 500., 0.1)
epsRange = np.arange(0.1, 1., 0.1)
x, y, z = convergeR(peRange, epsRange, angle)

fig = plt.figure()
ax = []
ax.append(fig.add_subplot(111, projection='3d'))
ax[0].plot_trisurf(x, y, z)
plt.show()


