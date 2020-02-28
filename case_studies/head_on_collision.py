'''
#                           This is an 80 character line                       #
We want to answer the question of what the equilibrium separation of two active
particles oriented at one another should be...
'''

import sys

# Run locally
sys.path.append('/Users/kolbt/Desktop/compiled/hoomd-blue/build')

import hoomd
from hoomd import md
from hoomd import dem
from hoomd import deprecated
from hoomd import data

import matplotlib
#matplotlib.use('Agg')
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
    tauLJ = ((sigma**2) * threeEtaPiSigma)
    return tauLJ
    
eps = kT                                # repulsive depth
tauLJ = computeTauLJ(eps)               # LJ time unit
dt = 0.000001 * tauLJ                   # timestep
dumpPerBrownian = 1.                   # number of dumps per 1 tauB
simLength = 1.0 * tauBrown            # how long to run (in tauBrown)
totTsteps = int(simLength / dt)         # how many tsteps to run
numDumps = simLength * dumpPerBrownian  # total number of frames dumped
dumpFreq = totTsteps / numDumps         # normalized dump frequency
dumpFreq = int(dumpFreq)                # ensure this is an integer
seed = 71996                            # a random seed
seed2 = 2394                            # orientation seed
seed3 = 183                             # activity seed

# Some parameters:

# Number of particles
N = 3.
# Particle positions
pos = [[-1., 0., 0.5], [0., 0., 0.5], [1., 0., 0.5]]
# Types
typ = [0, 0, 0]


# Now we make the system in hoomd
hoomd.context.initialize()
# A small shift to help with the periodic box
snap = hoomd.data.make_snapshot(N = 3,
                                box = hoomd.data.boxdim(Lx=100.,
                                                        Ly=100.,
                                                        dimensions=2),
                                particle_types = ['A'])

# Set positions/types for all particles
snap.particles.position[:] = pos[:]
snap.particles.typeid[:] = typ[:]

# Initialize the system
system = hoomd.init.read_snapshot(snap)
all = hoomd.group.all()

# Set particle potentials
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=2**(1/6), nlist=nl)
lj.set_params(mode='shift')
lj.pair_coeff.set('A', 'A', epsilon=eps, sigma=1.)

# Brownian integration
hoomd.md.integrate.mode_standard(dt=dt)
bd = hoomd.md.integrate.brownian(group=all, kT=0., seed=seed)

# Set body force of each particle explicity
pe = 100.
Fa_act = pe
Fb_act = 0
Fc_act = pe
aForce = (Fa_act, 0., 0.)
bForce = (-Fb_act, 0., 0.)
cForce = (-Fc_act, 0., 0.)
# Implement the activities in hoomd
agroup = hoomd.group.tag_list(name="a", tags=[0])
bgroup = hoomd.group.tag_list(name="b", tags=[1])
cgroup = hoomd.group.tag_list(name="c", tags=[2])
hoomd.md.force.constant(fvec=aForce, group=agroup)
hoomd.md.force.constant(fvec=bForce, group=bgroup)
hoomd.md.force.constant(fvec=cForce, group=cgroup)

# We need to put a hard wall on the left
leftWall = hoomd.md.wall.group(hoomd.md.wall.plane(origin=(-10, 0, 0),
                                                   normal=(1,0,0),
                                                   inside=True))
leftLJWall = hoomd.md.wall.slj(leftWall, r_cut=1.112)
leftLJWall.force_coeff.set('A', epsilon=10.0, sigma=1.0)

# Write dump
hoomd.dump.gsd('head_on_collision.gsd',
               period=20000,
               group=all,
               overwrite=True,
               phase=-1,
               dynamic=['attribute', 'property', 'momentum'])

# Run
hoomd.run(totTsteps)


