'''
#                        This is an 80 character line                          #
We want to simulate soft spheres getting sorted by a density gradient...

IN:
    -list of densities
    -layer thickness
    -particle activity
'''
# Initial imports
import sys
import os

sys.path.append('/Users/kolbt/Desktop/compiled/hoomd-blue/build')
#sys.path.append('/nas/longleaf/home/kolbt/programs/cpu-hoomd/hoomd-blue/build')
sys.path.append('/nas/longleaf/home/kolbt/programs/hoomd_2.2.1/hoomd-blue/build')

import hoomd
from hoomd import md
from hoomd import dem
from hoomd import deprecated
from hoomd import data
import numpy as np

import matplotlib.pyplot as plt

# Set some constants
kT = 1.0                        # temperature
threeEtaPiSigma = 1.0           # drag coefficient
sigma = 1.0                     # particle diameter
D_t = kT / threeEtaPiSigma      # translational diffusion constant
D_r = (3.0 * D_t) / (sigma**2)  # rotational diffusion constant
tauBrown = (sigma**2) / D_t     # brownian time scale (invariant)

# Repulsion and timestep
def computeTauLJ(epsilon):
    "Given epsilon, compute lennard-jones time unit"
    tauLJ = ((sigma**2) * threeEtaPiSigma) / epsilon
    return tauLJ
    
eps = kT                                # repulsive depth
tauLJ = computeTauLJ(eps)               # LJ time unit
dt = 0.000001 * tauLJ                   # timestep
dumpPerBrownian = 30.                   # number of dumps per 1 tauB
simLength = 20. * tauBrown              # how long to run (in tauBrown)
totTsteps = int(simLength / dt)         # how many tsteps to run
numDumps = simLength * dumpPerBrownian  # total number of frames dumped
dumpFreq = totTsteps / numDumps         # normalized dump frequency
dumpFreq = int(dumpFreq)                # ensure this is an integer
seed = 71996                            # a random seed
seed2= 2394                             # activity seed

# Parameters I'll need for this
densities = [1., 0.98, 0.96, 0.94, 0.92, 0.90]
# Thickness/height is in number of particles
thickness = 10.
height = 50
activity = 500.

def vertShift(lat):
    '''This will change depending on the density of the layer'''
    return np.sqrt(0.75) * lat
    
def horzShift(lat):
    '''Lattice-spacing dependent'''
    return lat / 2.

# Start at (0,0) for ease, slide all positions over once placed
pos = []
id = []
layer = 0
lshift = 0.
wallx = []
for i in densities:
    vert = vertShift(i)
    horz = horzShift(i)
    # Shift particles over according to previous lattices
    if i != densities[0]:
        lshift += ((thickness + 1.5) * i)
    wallx.append(lshift - 0.1)
    # Reset y layer
    ny = 0
    while ny <= height:
        # Should the layer get shifted
        if ny % 2 != 0:
            xshift = horz
        else:
            xshift = 0.
            
        # Reset x layer
        nx = 0
        while nx <= thickness:
            xpos = (nx * i) + xshift + lshift
            ypos = (ny * vert)
            zpos = 0.5
            pos.append([xpos, ypos, zpos])
            id.append(layer)
            # Increment x layer
            nx += 1
        # Increment y layer
        ny += 1
    layer += 1

# Let's give this a haircut (so that height is uniform)
xmin = -10.
xmax = pos[-1][0]
ymin = 0
ymax = pos[-1][1] + 0.5
for i in range(len(pos) - 1, -1, -1):
    if pos[i][0] > xmax:
        xmax = pos[i][0]
    if pos[i][1] > ymax:
        del pos[i]
        del id[i]
     
## Plotting script to check particle/wall placement
#xyz = zip(*pos)
#x = list(xyz[0])
#y = list(xyz[1])
#z = list(xyz[2])
#plt.xlim(xmin, xmax)
#plt.ylim(ymin, ymax)
#plt.scatter(x, y, c=id, s=5)
#for i in wallx:
#    plt.axvline(x=i)
#ax = plt.gca()
#ax.set_aspect('equal')
#plt.show()

xBox = xmax - xmin
yBox = ymax - ymin
for i in pos:
    i[0] -= ((xBox + (2. * xmin)) / 2.)
    i[1] -= (yBox / 2.)
for i in xrange(len(wallx)):
    wallx[i] -= ((xBox + (2. * xmin)) / 2.)
    
## Plotting script to check box adjustment
#xyz = zip(*pos)
#x = list(xyz[0])
#y = list(xyz[1])
#z = list(xyz[2])
#plt.scatter(x, y, c=id, s=5)
#for i in wallx:
#    plt.axvline(x=i)
#plt.xlim(-(xBox/2.), (xBox/2.))
#plt.ylim(-(yBox/2.), (yBox/2.))
#ax = plt.gca()
#ax.set_aspect('equal')
#plt.show()

nHex = len(pos)     # number of lattice particles
partNum = nHex + 1  # particle we are testing

# Now we make the system in hoomd
hoomd.context.initialize()
# A small shift to help with the periodic box
snap = hoomd.data.make_snapshot(N = partNum,
                                box = hoomd.data.boxdim(Lx=xBox,
                                                        Ly=yBox,
                                                        dimensions=2),
                                particle_types = ['A', 'B'])

# Get position for subject particle
xSub = ((nCol * lat) - hor) + 1.
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
lj.pair_coeff.set('A', 'A', epsilon=eps, sigma=sigma)
lj.pair_coeff.set('A', 'B', epsilon=eps, sigma=sigma)
lj.pair_coeff.set('B', 'B', epsilon=eps, sigma=sigma)

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

# Active motion initially oriented towards the HCP phase
activity = []
tuple = (-swimForce, 0, 0)
activity.append(tuple)
hoomd.md.force.active(group=active,
                      seed=seed2,
                      f_lst=activity,
                      rotation_diff=D_r,
                      orientation_link=False,
                      orientation_reverse_link=True)

# brownian integration
hoomd.md.integrate.mode_standard(dt=dt)
bd = hoomd.md.integrate.brownian(group=all, kT=kT, seed=seed)

out = "active_pe" + str(swimForce) + "_lattice" + str(lat) + ".gsd"

#write dump
hoomd.dump.gsd(out,
               period=dumpFreq,
               group=all,
               overwrite=True,
               phase=-1,
               dynamic=['attribute', 'property', 'momentum'])

#run
hoomd.run(totTsteps)
