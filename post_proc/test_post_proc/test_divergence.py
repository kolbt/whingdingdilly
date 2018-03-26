import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import math

def computeDivergence(field):
    "return the divergence of a n-D field"
    return np.sum(np.gradient(field),axis=0)

#interp = 20
#xx, yy = np.mgrid[-5:5:interp*1j, -5:5:interp*1j]

#T = np.arctan2(yy, xx)

eval = np.zeros((10, 10), dtype=np.float32)
#for i in range(0, interp):
#    for j in range(0, interp):
#        eval[i][j] = np.arctan2(yy[i], xx[j])
x = np.arange(-5, 5, 1)
y = np.arange(-5, 5, 1)
for iii in range(0, 10):
    for jjj in range(0, 10):
        eval[iii][jjj] = math.atan2(y[jjj], x[iii])


dx, dy = np.gradient(eval)

plt.subplot(221)
plt.imshow(dx.T, extent=(-5,5,-5,5), origin='lower')
plt.xlim(-5, 5)
plt.xticks(())
plt.ylim(-5, 5)
plt.yticks(())
#plt.colorbar()

plt.subplot(222)
plt.imshow(dy.T, extent=(-5,5,-5,5), origin='lower')
plt.xlim(-5, 5)
plt.xticks(())
plt.ylim(-5, 5)
plt.yticks(())
#plt.colorbar()

plt.subplot(223)
#plt.quiver(xx, yy, np.cos(T), np.sin(T), edgecolor='k', facecolor='None', linewidth=.5)
plt.quiver(x, y, eval[:][0], eval[0][:], edgecolor='k', facecolor='None', linewidth=.5)
plt.xlim(-5, 5)
plt.xticks(())
plt.ylim(-5, 5)
plt.yticks(())

div = computeDivergence(eval)

plt.subplot(224)
plt.imshow(div.T, extent=(-5,5,-5,5), origin='lower')
plt.xlim(-5, 5)
plt.xticks(())
plt.ylim(-5, 5)
plt.yticks(())
plt.show()

