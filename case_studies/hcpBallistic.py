'''
#                        This is an 80 character line                          #
We want to simulate soft spheres pushing through a stationary HCP wall...

    Reads in: 1.) # of columns in HCP phase
              2.) # of rows in HCP phase
              3.) lattice spacing
              4.) the swim force (for the ballistic particle)
              
    and runs a simulation accordingly :) Use this to sweep through
    lattice spacing vs force. Plot to understand penetration of faster
    particles into soft sphere dense phase.
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
dumpPerBrownian = 25.                   # number of dumps per 1 tauB
simLength = 20. * tauBrown              # how long to run (in tauBrown)
totTsteps = int(simLength / dt)         # how many tsteps to run
numDumps = simLength * dumpPerBrownian  # total number of frames dumped
dumpFreq = totTsteps / numDumps         # normalized dump frequency
dumpFreq = int(dumpFreq)                # ensure this is an integer
seed = 71996                            # a random seed

## Some parameters (from command line):
#nCol = int(sys.argv[1])         # number of lattice rows
#nRow = int(sys.argv[2])         # number of lattice columns
#lat = float(sys.argv[3])        # lattice spacing
#lat /= 100.                     # bash only works with integers
#swimForce = float(sys.argv[4])  # swim force of active particle

# Replace to have permanent input file
nCol = ${nCol}              # number of lattice rows
nCol = int(nCol)
nRow = ${nRow}              # number of lattice columns
nRow = int(nRow)
lat = ${latCount}           # lattice spacing
lat = float(lat)
lat /= 100.                 # bash only works with integers
swimForce = ${swimCount}    # swim force of active particle
swimForce = float(swimForce)

ver = np.sqrt(0.75) * lat       # vertical shift between lattice rows
hor = lat / 2.0                 # horizontal shift between lattice rows

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

# Ballistic motion towards the HCP phase
force = (-swimForce, 0, 0)
hoomd.md.force.constant(fvec=force, group=active)

# Alternatively we can make it a true active particle

# brownian integration
hoomd.md.integrate.mode_standard(dt=dt)
bd = hoomd.md.integrate.brownian(group=all, kT=kT, seed=seed)

out = "pe" + str(swimForce) + "_lattice" + str(lat) + ".gsd"

#write dump
hoomd.dump.gsd(out,
               period=dumpFreq,
               group=all,
               overwrite=True,
               phase=-1,
               dynamic=['attribute', 'property', 'momentum'])

#run
hoomd.run(totTsteps)
