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

#import matplotlib.pyplot as plt

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
# Must be multiple of 3 (to account for tau_r)
dumpPerBrownian = 12.                   # number of dumps per 1 tauB
simLength = 100. * tauBrown             # how long to run (in tauBrown)
totTsteps = int(simLength / dt)         # how many tsteps to run
numDumps = simLength * dumpPerBrownian  # total number of frames dumped
dumpFreq = totTsteps / numDumps         # normalized dump frequency
dumpFreq = int(dumpFreq)                # ensure this is an integer
seed = 71996                            # a random seed
seed2= ${inSeed}                        # activity seed

# Parameters I'll need for this
densities = [1., 0.98, 0.96, 0.94, 0.92, 0.90]
# Thickness/height is in number of particles
thickness = 20.
height = 100
swimForce = ${swimCount}

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

# Get the unique type IDs
uniqueID = []
uniqueChar = []
for i in id:
    if i not in uniqueID:
        uniqueID.append(i)
        uniqueChar.append(chr(ord('@') + i+1))
# Don't forget to add the active particle
uniqueID.append(max(uniqueID) + 1)
uniqueChar.append( chr(ord('@') + uniqueID[-1]+1) )

# Convert all IDs to chars
chars = []
for i in id:
    chars.append(chr(ord('@') + i+1))

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
relax = 0.5
for i in pos:
    i[0] -= ((xBox + (2. * xmin)) / 2.)
    i[1] -= (yBox / 2.)
for i in xrange(len(wallx)):
    wallx[i] -= ((xBox + (2. * xmin)) / 2.)
    wallx[i] -= 0.1
wallx.append((xBox / 2.) + (relax / 2.))
    
#print(len(wallx))
## Plotting script to check box adjustment
#xyz = zip(*pos)
#x = list(xyz[0])
#y = list(xyz[1])
#z = list(xyz[2])
#plt.scatter(x, y, c=id, s=5)
#for i in wallx:
#    plt.axvline(x=i)
#plt.xlim(-(xBox/2.) - 5, (xBox/2.) + 5)
#plt.ylim(-(yBox/2.) - 5, (yBox/2.) + 5)
#ax = plt.gca()
#ax.set_aspect('equal')
#plt.show()

nHex = len(pos)     # number of lattice particles
partNum = nHex + 1  # particle we are testing

# Now we make the system in hoomd
hoomd.context.initialize()
# A small shift to help with the periodic box
snap = hoomd.data.make_snapshot(N = partNum,
                                box = hoomd.data.boxdim(Lx=xBox + relax,
                                                        Ly=yBox + relax,
                                                        dimensions=2),
                                particle_types = uniqueChar)

# Get position/id/char for subject particle
xSub = -(xBox/2.) + 5.
ySub = 0.
subject = [xSub, ySub, 0.5]
pos.append(subject)
id.append(uniqueID[-1])
chars.append(uniqueChar[-1])

# Set positions/types for all particles
snap.particles.position[:] = pos[:]
snap.particles.typeid[:] = id[:]
snap.particles.types[:] = chars[:]

# Initialize the system
system = hoomd.init.read_snapshot(snap)
all = hoomd.group.all()
active = hoomd.group.type(type=uniqueChar[-1])

# Set particle potentials
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=2**(1/6), nlist=nl)
lj.set_params(mode='shift')
for i in xrange(len(uniqueChar)):
    for j in range(i, len(uniqueChar)):
        lj.pair_coeff.set(uniqueChar[i],
                          uniqueChar[j],
                          epsilon=eps, sigma=sigma)

# Add wall potentials that maintain density gradient
walls = []
wallPots = []
for i in xrange(len(wallx)):
    if i != (len(wallx) - 1):
        walls.append(hoomd.md.wall.group(hoomd.md.wall.plane(origin=(wallx[i], 0, 0),
                                                             normal=(1,0,0),
                                                             inside=True)))
    else:
        walls.append(hoomd.md.wall.group(hoomd.md.wall.plane(origin=(wallx[i], 0, 0),
                                                             normal=(1,0,0),
                                                             inside=False)))
    wallPots.append(hoomd.md.wall.slj(walls[i], r_cut=(2**(1./6.)) ))
    for j in xrange(len(uniqueChar)):
        # Final wall needs to bound the particles on the right
        if  i == (len(wallx) - 1) and j == (len(uniqueChar) - 2):
            wallPots[i].force_coeff.set(uniqueChar[j], epsilon=10.0, sigma=1.0)
        else:
            # Each wall will only interact with one particle type (not the active particle)
            if i == j and j != (len(uniqueChar) - 1):
                wallPots[i].force_coeff.set(uniqueChar[j], epsilon=10.0, sigma=1.0)
            # Other interactions are set to zero
            else:
                wallPots[i].force_coeff.set(uniqueChar[j], epsilon=0.0, sigma=1.0)

# Brownian integration
brownEquil = 10000
hoomd.md.integrate.mode_standard(dt=dt)
bd = hoomd.md.integrate.brownian(group=all, kT=kT, seed=seed)
hoomd.run(brownEquil)

# Active motion initially oriented towards the HCP phase
activity = []
tuple = (swimForce, 0, 0)
activity.append(tuple)
hoomd.md.force.active(group=active,
                      seed=seed2,
                      f_lst=activity,
                      rotation_diff=D_r,
                      orientation_link=False,
                      orientation_reverse_link=True)

out = "gradient_density_pe" + str(swimForce) + ".gsd"

# Write output
hoomd.dump.gsd(out,
               period=dumpFreq,
               group=all,
               overwrite=True,
               phase=-1,
               dynamic=['attribute', 'property', 'momentum'])

# Run
hoomd.run(totTsteps)
