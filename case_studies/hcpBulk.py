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
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

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
seed = 71996                            # brownian integration seed
seed2 = 9376                            # orientation seed
seed3 = 2394                            # activity seed

## Some parameters (from command line):
#nCol = int(sys.argv[1])         # number of lattice rows
#nRow = nCol                     # keep square box
#lat = float(sys.argv[2])        # lattice spacing
#lat /= 100.                     # bash only works with integers
#swimForce = float(sys.argv[3])  # swim force of active particle

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
xBox = (nCol * lat)
yBox = (nRow * ver)

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
partNum = nHex      # change on lattice particle to the active particle

# Now we make the system in hoomd
hoomd.context.initialize()
# A small shift to help with the periodic box
per = ver
snap = hoomd.data.make_snapshot(N = partNum,
                                box = hoomd.data.boxdim(Lx=xBox,
                                                        Ly=yBox,
                                                        dimensions=2),
                                particle_types = ['A', 'B'])

actInd = (nHex / 2) + (nRow / 2) - 1
# Incorporate ID and type into all particle array
myID = []
myTyp = []
for i in xrange(partNum):
    if i != actInd:
        myID.append(0)
        myTyp.append('A')
    else:
        myID.append(1)
        myTyp.append('B')
        
# You have to shift the positions
xShift = (xBox / 2.0)
yShift = (yBox / 2.0)
for i in xrange(len(pos)):
    pos[i][0] -= (xShift - (hor / 2.))
    pos[i][1] -= (yShift - (ver / 2.))
    
## Plot to check box and positions
#pltx, plty, pltz = zip(*pos)
#plt.scatter(pltx, plty, c=myID)
#ax = plt.gca()
#ax.add_patch(Rectangle((-xShift, -yShift), xBox, yBox, ec='k', fc='none', lw=1.5))
#ax.set_aspect('equal')
#plt.xlim(-xShift - 5, xShift + 5)
#plt.ylim(-yShift - 5, yShift + 5)
#plt.show()

# Set positions/types for all particles
snap.particles.position[:] = pos[:]
snap.particles.typeid[:] = myID[:]
snap.particles.types[:] = myTyp[:]

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

# Active motion (with random initial direction)
activity = []
np.random.seed(seed2)
randomOrient = np.random.rand(1) * 2 * np.pi
x = (np.cos(randomOrient)) * swimForce
y = (np.sin(randomOrient)) * swimForce
z = 0.
tuple = (x, y, z)
activity.append(tuple)
hoomd.md.force.active(group=active,
                      seed=seed3,
                      f_lst=activity,
                      rotation_diff=D_r,
                      orientation_link=False,
                      orientation_reverse_link=True)

# Brownian integration
hoomd.md.integrate.mode_standard(dt=dt)
bd = hoomd.md.integrate.brownian(group=all, kT=kT, seed=seed)

out = "active_in_bulk_pe" + str(swimForce) + "_lattice" + str(lat) + ".gsd"
outBulk = "bulkHCP_pe" + str(swimForce) + "_lattice" + str(lat) + ".gsd"

# Write dump
hoomd.dump.gsd(out,
               period=1,
               group=active,
               overwrite=True,
               phase=-1,
               dynamic=['attribute', 'property', 'momentum'])

# Run
hoomd.run(totTsteps)
