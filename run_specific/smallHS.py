'''
#                           This is an 80 character line                       #
Purpose: run MONODISPERSE hard spheres to approximate the reentrant effect
observed for soft particles. Note that the actual activity of these particles
is much higher than is being reported (we maintain units in terms of the larger
diameter, sigma=1.0)
'''
# Initial imports
import sys
import os
import psutil

# Read in bash arguments
hoomdPath = "${hoomd_path}"     # path to hoomd-blue
gsdPath = "${gsd_path}"         # path to gsd
runFor = ${runfor}              # simulation length (in tauLJ)
dumpPerBrownian = ${dump_freq}  # how often to dump data
pe = ${pe}                      # activity of A particles
partNum = ${part_num}           # total number of particles
intPhi = ${phi}                 # system area fraction
phi = float(intPhi)/100.0

seed1 = ${seed1}                # seed for position
seed2 = ${seed2}                # seed for bd equilibration
seed3 = ${seed3}                # seed for initial orientations
seed4 = ${seed4}                # seed for A activity

# Remaining imports
sys.path.append(hoomdPath)
import hoomd
from hoomd import md
from hoomd import deprecated
import numpy as np

# Set some constants
kT = 1.0                        # temperature
threeEtaPiSigma = 1.0           # drag coefficient
sigma = 1.0                     # particle diameter
D_t = kT / threeEtaPiSigma      # translational diffusion constant
D_r = (3.0 * D_t) / (sigma**2)  # rotational diffusion constant
tauBrown = (sigma**2) / D_t     # brownian time scale (invariant)

def computeVel(activity):
    "Given particle activity, output intrinsic swim speed"
    velocity = (activity * sigma) / (3 * (1/D_r))
    return velocity

def computeActiveForce(velocity):
    "Given particle activity, output repulsion well depth"
    activeForce = velocity * threeEtaPiSigma
    return activeForce

def computeNetEps(activeForce):
    "Given particle activity, output repulsion well depth"
    # The minimum overall epsilon should be Brownian
    epsBrown = 1.0
    epsilon = ((4.0 * activeForce) / 24.0) + epsBrown
    return epsilon

def computeTauLJ(epsilon):
    "Given epsilon, compute lennard-jones time unit"
    tauLJ = ((sigma**2) * threeEtaPiSigma) / epsilon
    return tauLJ

def computeSigma(activity):
    "We fitted a power law relation to mode separation vs activity"
    m = -0.0958
    b = 0.3068
    effectiveSigma = (activity**m) * (np.exp(b))
    return effectiveSigma

# Compute parameters from activities
if pe != 0:                        # A particles are NOT Brownian
    v = computeVel(pe)
    Fp = computeActiveForce(v)
    eps = computeNetEps(Fp)
    effSig = computeSigma(pe)
else:                               # A particles are Brownian
    v = 0.0
    Fp = 0.0
    eps = 10.0
    effSig = 1.0

tauLJ = computeTauLJ(eps)                       # get LJ time unit
cut = (2**(1./6.)) * effSig                     # the cutoff for the LJ potential

dt = 0.00001 * tauLJ                            # timestep size
simLength = runFor * tauBrown                   # how long to run (in tauBrown)
simTauLJ = simLength / tauLJ                    # how long to run (in tauLJ)
totTsteps = int(simLength / dt)                 # how many tsteps to run
numDumps = float(simLength * dumpPerBrownian)   # frames in 1 tauB
dumpFreq = float(totTsteps / numDumps)          # normalized dump frequency
dumpFreq = int(dumpFreq)                        # ensure this is an integer

print "Brownian tau in use:", tauBrown
print "Lennard-Jones tau in use:", tauLJ
print "Timestep in use:", dt
print "Epsilon in use:", eps
print "Total number of timesteps:", totTsteps
print "Total number of output frames:", numDumps
print "Dumped snapshots per 1 tauB:", dumpPerBrownian
print "Brownian run time:", simLength
print "Activity:", pe
print "Effective diameter:", effSig

# Initialize system
hoomd.context.initialize()
# We can still use phi_p as input, the radius is assumed to be 0.5
system = hoomd.deprecated.init.create_random(N = partNum,
                                             phi_p = phi,
                                             name = 'A',
                                             min_dist = 0.70,
                                             seed = seed1,
                                             dimensions = 2)

# Assigning groups and lengths to particles
all = hoomd.group.all()
N = len(all)

# Define potential between pairs
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=cut, nlist=nl)
lj.set_params(mode='shift')
lj.pair_coeff.set('A', 'A', epsilon=eps, sigma=effSig)

# General integration parameters
brownEquil = 100000
hoomd.md.integrate.mode_standard(dt=dt)
hoomd.md.integrate.brownian(group=all, kT=kT, seed=seed2)
hoomd.run(brownEquil)

#set the activity of each type
np.random.seed(seed3)                           # seed for random orientations
angle = np.random.rand(partNum) * 2 * np.pi     # random particle orientation

activity = []
for i in range(0,partNum):
    x = (np.cos(angle[i])) * pe
    y = (np.sin(angle[i])) * pe
    z = 0
    tuple = (x, y, z)
    activity.append(tuple)
hoomd.md.force.active(group=all,
                      seed=seed4,
                      f_lst=activity,
                      rotation_diff=D_r,
                      orientation_link=False,
                      orientation_reverse_link=True)

# Get filenames for various file types
name = "pe" + str(pe) +\
"_ep" + str(int(eps)) +\
"_phi" + str(intPhi)

gsdName = name + ".gsd"

hoomd.dump.gsd(gsdName,
               period=dumpFreq,
               group=all,
               overwrite=False,
               phase=-1,
               dynamic=['attribute', 'property', 'momentum'])

hoomd.run(totTsteps)
