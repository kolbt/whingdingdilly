'''
#                           This is an 80 character line                       #
We want to simulate soft spheres pushing through a stationary HCP wall
'''
# Initial imports
import sys
import os

sys.path.append('/Users/kolbt/Desktop/compiled/hoomd-blue/build')

import hoomd
from hoomd import md
from hoomd import dem
from hoomd import deprecated
from hoomd import data
import numpy as np

# Some parameters:
nCol = 100                  # number of lattice rows
nRow = 20                   # number of lattice columns
lat = 1.0                   # lattice spacing
ver = np.sqrt(0.75) * lat   # vertical shift between lattice rows
hor = lat / 2.0             # horizontal shift between lattice rows

# Grab the box coordinates
xBox = (nCol * lat) - hor
xBox *= 1.5
yBox = (nRow * ver) - ver

# Make the mesh find positions for the HCP particles
pos = []
y = 0
for i in xrange(nRow):
    # Every other row needs to be shifted in the x direction
    if i % 2 == 0:
        x = 0.0
    else:
        x = hor

    for j in xrange(nCol):
        pos.append([x, y, 0.5]) # place particle
        x += lat                # move a lattice spacing to the right

    # Every row increases by the same y-value
    y += ver

nHex = len(pos)     # number of lattice particles
partNum = nHex + 1  # particle we are testing

# Plot to check and make sure things are working
#import matplotlib.pyplot as plt
#xPos, yPos, zPos = zip(*pos)
#plt.scatter(xPos, yPos)
#plt.xlim(0, xBox)
#plt.ylim(0, yBox)
#plt.show()

# Now we make the system in hoomd
hoomd.context.initialize()
# A small shift to help with the periodic box
per = ver
snap = hoomd.data.make_snapshot(N = partNum,
                                box = hoomd.data.boxdim(Lx=xBox + per,
                                                        Ly=yBox + per,
                                                        dimensions=2),
                                particle_types = ['A', 'B'])

# Get position for subject particle
xSub = ((nCol * lat) - hor) * 1.05
ySub = yBox / 2.0
subject = [xSub, ySub, 0.5]
pos.append(subject)

# You have to shift the positions
xShift = xBox / 2.0
yShift = yBox / 2.0
for i in xrange(len(pos)):
    pos[i][0] -= xShift
    pos[i][1] -= yShift

# Set positions/types for all particles
snap.particles.position[:] = pos[:]
snap.particles.typeid[:] = 0
snap.particles.types[:] = 'A'
snap.particles.typeid[-1] = 1
snap.particles.types[-1] = 'B'

# Initialize the system
system = hoomd.init.read_snapshot(snap)
all = hoomd.group.all()
hcp = hoomd.group.type(type='A')
active = hoomd.group.type(type='B')

# Set particle potentials
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=2**(1/6), nlist=nl)
lj.set_params(mode='shift')
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('A', 'B', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0)

# Add a wall potential that acts on the HCP particles
wallShift = 5.0 * lat
wallShift = 0.0
xyzLWall = (((xBox + per) / -2.0), 0, 0)
xyzRWall = ((nCol * lat) - xShift + wallShift, 0, 0)
leftWall = hoomd.md.wall.group(hoomd.md.wall.plane(origin=xyzLWall,
                                                   normal=(1,0,0),
                                                   inside=True))
rightWall = hoomd.md.wall.group(hoomd.md.wall.plane(origin=xyzRWall,
                                                    normal=(1,0,0),
                                                    inside=False))
leftLJWall = hoomd.md.wall.slj(leftWall, r_cut=1.112)
rightLJWall = hoomd.md.wall.slj(rightWall, r_cut=1.112)
leftLJWall.force_coeff.set('A', epsilon=10.0, sigma=1.0)
leftLJWall.force_coeff.set('B', epsilon=0.0, sigma=1.0)
rightLJWall.force_coeff.set('A', epsilon=10.0, sigma=1.0)
rightLJWall.force_coeff.set('B', epsilon=0.0, sigma=1.0)

# Ballistic motion towards the HCP phase
force = (-500, 0, 0)
hoomd.md.force.constant(fvec=force, group=active)

# Alternatively we can make it a true active particle

# brownian integration
hoomd.md.integrate.mode_standard(dt=0.00001)
bd = hoomd.md.integrate.brownian(group=all, kT=1.0, seed=123)

#write dump
hoomd.dump.gsd("test_hex.gsd",
               period=1000,
               group=all,
               overwrite=True,
               phase=-1,
               dynamic=['attribute', 'property', 'momentum'])

#run
hoomd.run(100000)
