'''
A schematic figure to illustrate how Pe (and F_act) sets deformability:
    -LJ potential
    -Overlay 2 Forces (Strong under week)
    -Corresponds to collision angle
'''

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d

eps = 0.01
sigma = 1.
def ljPotential(r, eps=0.1, sigma=1.):
    div = (sigma/r)
    U = ( 4. * eps * ((div)**12 - (div)**6) ) + eps
    return U

def ljForce(r, eps=0.1, sigma=1.):
    div = (sigma/r)
    dU = (24. * eps / r) * ((2*(div**12)) - (div)**6)
    return dU
    
def collisionForce(pe, angle):
    return pe - (pe * np.cos(angle))
    
from sympy.solvers import solve
from sympy import solveset, S
from sympy.abc import x
from sympy import Symbol
def forceToDist(f, eps=0.1, sigma=1.):
    r = Symbol('r', real=True, positive=True)
    solution = solve(( (f/(24. * eps)) * (r**13.) ) + ( 2. * (sigma**6.) * (r**6.) ) - ( 2 * (sigma**12) ), r, numerical=True)
#    solution = solve([r>=0.5, r<=1.0, ( (f/(24. * eps)) * (r**13.) ) + ( 2. * (sigma**6.) * (r**6.) ) - ( 2 * (sigma**12) )], r)
#    solution = solveset(( (f/(24. * eps)) * (x**13.) ) + ( 2. * (sigma**6.) * (x**6.) ) - ( 2 * (sigma**12) ), x, domain=S.Reals)
    return solution[0]

# Angle of collision (left particle points toward right particle)
angle = 0.      # right particle pointing away (no deformation)
angle = np.pi   # head on collision
# Compute the weak and strong collision force
peWeak = 20.
peStrong = 500.
# Assume that the left particle always points toward the right
fWeak = collisionForce(peWeak, angle)
fStrong = collisionForce(peStrong, angle)

# Compute the distance that corresponds to the force
#rWeak = forceToDist(fWeak, eps=eps)
#print("Weak force: r={}").format(rWeak)
#rWeak = 0.776736185849486 # pe = 50
#rWeak = 0.825094041592472 # pe = 20
rWeak = 0.704511014939217

#rStrong = forceToDist(fStrong, eps=eps)
#print("Strong force: r={}").format(rStrong)
#rStrong = 0.719245773085951 # pe = 150
#rStrong = 0.658845113101655 # pe = 500
rStrong = 0.554278202533698


fig = plt.figure()
ax = []
ax.append(fig.add_subplot(221))                     #top left
ax.append(fig.add_subplot(222, projection='3d'))    #top right
ax.append(fig.add_subplot(223))                     #bottom left
ax.append(fig.add_subplot(224, projection='3d'))    #bottom right
# Plot the LJ potential and Force:
# Base plot
dist = np.arange(0.0001, ((2.**(1./6.))*sigma)*2., 0.001)
ax[0].plot(dist, ljPotential(dist, eps=eps), c='k', lw=5., label='LJ-Potential')
ax[2].plot(dist, ljForce(dist, eps=eps), c='k', lw=5., label='LJ-Force')
# Plot for PeStrong
strongRange = np.arange(rStrong, ((2.**(1./6.))*sigma)*2., 0.001)
ax[0].plot(strongRange, ljPotential(strongRange, eps=eps), c='g', lw=2.5, label='Strong')
ax[2].plot(strongRange, ljForce(strongRange, eps=eps), c='g', lw=2.5, label='Strong')
# Plot for PeWeak
weakRange = np.arange(rWeak, ((2.**(1./6.))*sigma)*2., 0.001)
ax[0].plot(weakRange, ljPotential(weakRange, eps=eps), c='b', lw=1.25, label='Weak')
ax[2].plot(weakRange, ljForce(weakRange, eps=eps), c='b', lw=1.25, label='Weak')

# Limits
ax[0].set_xlim(0.65, (2.**(1./6.))*sigma)
ax[2].set_xlim(0.65, (2.**(1./6.))*sigma)
ax[0].set_ylim(0., 300.*eps)
ax[2].set_ylim(0., 3.*peStrong)
ax[0].legend()
ax[2].legend()

# Plot the overlap of spheres
# For wire mesh
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
# For solid sphere
uS, vS = np.mgrid[0:2*np.pi:1000j, 0:np.pi:500j]
xS = np.cos(u)*np.sin(v)
yS = np.sin(u)*np.sin(v)
zS = np.cos(v)
ax[1].plot_wireframe(x - rWeak, y, z, color="#808080")
ax[1].plot_wireframe(x + rWeak, y, z, color="#808080")
ax[1].plot_surface((xS*rWeak) - rWeak, yS*rWeak, zS*rWeak, color="b")
ax[1].plot_surface((xS*rWeak) + rWeak, yS*rWeak, zS*rWeak, color="b")
ax[1].set_axis_off()
ax[1].view_init(0, 90)
#ax[1].set_xlim(-2.*rWeak, 2.*rWeak)
#ax[1].set_ylim(-1.5*rWeak, 1.5*rWeak)
#ax[1].set_zlim(-1.5*rWeak, 1.5*rWeak)
ax[1].set_xlim(-2., 2.)
ax[1].set_ylim(-1.5, 1.5)
ax[1].set_zlim(-1.5, 1.5)
ax[1].dist = 5.

ax[3].plot_wireframe(x - rStrong, y, z, color="#808080")
ax[3].plot_wireframe(x + rStrong, y, z, color="#808080")
ax[3].plot_surface((xS*rStrong) - rStrong, yS*rStrong, zS*rStrong, color="g")
ax[3].plot_surface((xS*rStrong) + rStrong, yS*rStrong, zS*rStrong, color="g")
ax[3].set_axis_off()
ax[3].view_init(0, 90)
#ax[3].set_xlim(-2.*rStrong, 2.*rStrong)
#ax[3].set_ylim(-1.5*rStrong, 1.5*rStrong)
#ax[3].set_zlim(-1.5*rStrong, 1.5*rStrong)
ax[3].set_xlim(-2., 2.)
ax[3].set_ylim(-1.5, 1.5)
ax[3].set_zlim(-1.5, 1.5)
ax[3].dist = 5.

plt.savefig("particle_deformation_eps" + str(eps) + ".png", dpi=1000, bbox_inches='tight', pad_inches=0)
plt.close()
