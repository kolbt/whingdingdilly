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
seed2= 2394                             # activity seed

# Some parameters (from command line):
pe1 = float(sys.argv[1])
pe2 = float(sys.argv[2])
pe3 = float(sys.argv[3])
r1 = float(sys.argv[4])
w2 = float(sys.argv[5])
w3 = float(sys.argv[6])

def computeLat(activity):
    "Get lattice spacing based on effective diameter"
    m = -0.0957707250171
    b = 0.306774161185
    lattice = (activity**m) * (np.exp(b))
    return lattice
    
def computeDistance(x, y):
    return np.sqrt((x**2) + (y**2))
    
# List of activities
peList = [ pe1, pe2, pe3 ]
# List of ring radii
rList = [ 0, r1, r1 + w2, r1 + w2 + w3 ]
# List to store particle positions and types
pos = []
typ = []

for i in xrange(len(peList)):
    rMin = rList[i]             # starting distance for particle placement
    rMax = rList[i + 1]         # maximum distance for particle placement
    lat = computeLat(pe)        # activity-dependent lattice spacing
    ver = np.sqrt(0.75) * lat   # vertical shift between lattice rows
    hor = lat / 2.0             # horizontal shift between lattice rows
    
    x = rMin
    y = rMin
    shift = 0
    while y < rMax:
        r = computeDistance(x, y)
        # Check if x-position is large enough
        if r < rMin:
            x += lat
            continue
            
        # Check if x-position is too large
        if r >= rMax:
            y += ver
            shift += 1
            if shift % 2:
                x = hor
            else:
                x = 0
            continue
        
        # If the loop makes it this far, append
        pos.append((x, y))
        if x != 0 and y != 0:
            pos.append((-x, y))
            pos.append((-x, -y))
            pos.append((x, -y))
        # y must be zero
        elif x != 0:
            pos.append((x, -y))
        # x must be zero
        elif y!= 0:
            pos.append((-x, y))
        
            
        
        
        
    
    

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
