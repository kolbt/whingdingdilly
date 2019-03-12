'''
#                           This is an 80 character line                       #
    This is intended to investigate how hard our particles need to be. We want
    to maintain a ratio of active force to LJ well-depth:
                    epsilon = alpha * F_A * sigma / 24.0
    This code will investigate alpha, in order to find the smallest value that
    will maintain a "hard-sphere" potential. (This will optimize computation
    time while keeping our model realistic)
'''
# Initial imports
import sys
import os
import psutil

# Read in bash arguments
hoomdPath = "${hoomd_path}"     # path to hoomd-blue
gsdPath = "${gsd_path}"         # path to gsd
runFor = ${runfor}              # simulation length (in tauLJ)
dumpFreq = ${dump_freq}         # how often to dump data
partPercA = ${part_frac_a}      # percentage of A particles
partFracA = float(partPercA) / 100.0
peA = ${pe_a}                   # activity of A particles
peB = ${pe_b}                   # activity of B particles
partNum = ${part_num}           # total number of particles
intPhi = ${phi}                 # system area fraction
phi = float(intPhi) / 100.0

seed1 = ${seed1}                # seed for position
seed2 = ${seed2}                # seed for bd equilibration
seed3 = ${seed3}                # seed for initial orientations
seed4 = ${seed4}                # seed for A activity
seed5 = ${seed5}                # seed for B activity

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
    # This gives:
    # v_0 = Pe * sigma / tau_B = Pe * sigma / 3 * tau_R
    velocity = (activity * sigma) / (3 * (1/D_r))
    return velocity

def computeActiveForce(velocity):
    "Given particle activity, output repulsion well depth"
    # This is multiplied by Brownian time and gives:
    #          Pe = 3 * v_0 * tau_R / sigma
    # the conventional description of the Peclet number
    activeForce = velocity * threeEtaPiSigma
    return activeForce

def computeEps(alpha, activeForce):
    "Given particle activity, output repulsion well depth"
    # Here is where we will be testing the ratio we use (via alpha)
    epsilon = (alpha * activeForce * sigma / 24.0) + 1.0
    # Add 1 because of integer rounding
    epsilon = int(epsilon) + 1
    return epsilon

def computeTauLJ(epsilon):
    "Given epsilon, compute lennard-jones time unit"
    tauLJ = ((sigma**2) * threeEtaPiSigma) / epsilon
    return tauLJ

# Compute parameters from activities
if peA != 0:                        # A particles are NOT Brownian
    vA = computeVel(peA)
    FpA = computeActiveForce(vA)
    epsA = computeEps(FpA)
    tauA = computeTauLJ(epsA)
else:                               # A particles are Brownian
    vA = 0.0
    FpA = 0.0
    epsA = kT
    tauA = computeTauLJ(epsA)

if peB != 0:                        # B particles are NOT Brownian
    vB = computeVel(peB)
    FpB = computeActiveForce(vB)
    epsB = computeEps(FpB)
    tauB = computeTauLJ(epsB)
else:                               # B particles are Brownian
    vB = 0.0
    FpB = 0.0
    epsB = kT
    tauB = computeTauLJ(epsB)

#epsAB = (epsA + epsB) / 2.0                # AB interaction well depth
tauLJ = (tauA if (tauA <= tauB) else tauB)  # use the smaller tauLJ
epsA = (epsA if (epsA >= epsB) else epsB)   # use the larger epsilon
epsB = epsA                                 # make sure all use this
epsAB = epsA                                # make sure all use this
dt = 0.00001 * tauLJ                        # timestep size
simLength = runFor * tauBrown               # how long to run (in tauBrown)
simTauLJ = simLength / tauLJ                # how long to run (in tauLJ)
totTsteps = int(simLength / dt)             # how many tsteps to run
numDumps = float(simLength / 0.005)         # dump data every 0.5 tauBrown
dumpFreq = float(totTsteps / numDumps)      # normalized dump frequency
dumpFreq = int(dumpFreq)                    # ensure this is an integer

print "Brownian tau in use:", tauBrown
print "Lennard-Jones tau in use:", tauLJ
print "Timestep in use:", dt
print "Epsilon in use:", epsAB
print "Total number of timesteps:", totTsteps
print "Total number of output frames:", numDumps
print "File dump frequency:", dumpFreq

# Initialize system
hoomd.context.initialize()
system = hoomd.deprecated.init.create_random(N = partNum,
                                             phi_p = phi,
                                             name = 'A',
                                             min_dist = 0.70,
                                             seed = seed1,
                                             dimensions = 2)
# Add B-type particles
system.particles.types.add('B')
snapshot = system.take_snapshot()
partA = partNum * partFracA         # get the total number of A particles
partA = int(partA)                  # make sure it is an integer
partB = partNum - partA             # get the total number of B particles
partB = int(partB)                  # make sure it is an integer
mid = int(partA)                    # starting index to assign B particles

# Assign particles a type within the snapshot
if partPercA == 0:                      # take care of all b case
    mid = 0
    for i in range(mid, partNum):
        system.particles[i].type = 'B'
elif partPercA != 100:                  # mix of each or all A
    for i in range(mid, partNum):
        system.particles[i].type = 'B'

# Assigning groups and lengths to particles
all = hoomd.group.all()
gA = hoomd.group.type(type = 'A', update=True)
gB = hoomd.group.type(type = 'B', update=True)
N = len(all)
Na = len(gA)
Nb = len(gB)

# Define potential between pairs
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=2**(1/6), nlist=nl)
lj.set_params(mode='shift')
lj.pair_coeff.set('A', 'A', epsilon=epsA, sigma=1.0)
lj.pair_coeff.set('A', 'B', epsilon=epsAB, sigma=1.0)
lj.pair_coeff.set('B', 'B', epsilon=epsB, sigma=1.0)

# General integration parameters
brownEquil = 100000
hoomd.md.integrate.mode_standard(dt=dt)
hoomd.md.integrate.brownian(group=all, kT=kT, seed=seed2)
hoomd.run(brownEquil)

#set the activity of each type
np.random.seed(seed3)                           # seed for random orientations
angle = np.random.rand(partNum) * 2 * np.pi    # random particle orientation

# Case 1: Mixture
if partPercA != 0 and partPercA != 100:
    # First assign A-type active force vectors (w/ peA)
    activity_a = []
    for i in range(0,mid):
        x = (np.cos(angle[i])) * peA    # x active force vector
        y = (np.sin(angle[i])) * peA    # y active force vector
        z = 0                           # z active force vector
        tuple = (x, y, z)               # made into a tuple
        activity_a.append(tuple)        # add to activity A list

    # Now assign B-type active force vectors (w/ peB)
    activity_b = []
    for i in range(mid,partNum):
        x = (np.cos(angle[i])) * peB
        y = (np.sin(angle[i])) * peB
        z = 0
        tuple = (x, y, z)
        activity_b.append(tuple)
    # Set A-type activity in hoomd
    hoomd.md.force.active(group=gA,
                          seed=seed4,
                          f_lst=activity_a,
                          rotation_diff=D_r,
                          orientation_link=False,
                          orientation_reverse_link=True)
    # Set B-type activity in hoomd
    hoomd.md.force.active(group=gB,
                          seed=seed5,
                          f_lst=activity_b,
                          rotation_diff=D_r,
                          orientation_link=False,
                          orientation_reverse_link=True)
else:
    # Case 2: All B system
    if partPercA == 0:
        activity_b = []
        for i in range(0,partNum):
            x = (np.cos(angle[i])) * peB
            y = (np.sin(angle[i])) * peB
            z = 0
            tuple = (x, y, z)
            activity_b.append(tuple)
        hoomd.md.force.active(group=gB,
                              seed=seed5,
                              f_lst=activity_b,
                              rotation_diff=D_r,
                              orientation_link=False,
                              orientation_reverse_link=True)
    # Case 3: All A system
    else:
        activity_a = []
        for i in range(0,partNum):
            x = (np.cos(angle[i])) * peA
            y = (np.sin(angle[i])) * peA
            z = 0
            tuple = (x, y, z)
            activity_a.append(tuple)
        hoomd.md.force.active(group=gA,
                              seed=seed4,
                              f_lst=activity_a,
                              rotation_diff=D_r,
                              orientation_link=False,
                              orientation_reverse_link=True)

# Get filenames for various file types
name = "pa" + str(peA) +\
"_pb" + str(peB) +\
"_xa" + str(partPercA) +\
"_ep" + str(int(epsAB))

gsdName = name + ".gsd"
sqliteName = name + ".sqlite"
logName = "log_" + name + ".gsd"
chkName = os.getcwd() + "/" + logName
lastTen = "lastTen" + name + ".gsd"

### Dump for MSD ###
# Get the early times for MSD output
# Time intervals I want, in tauBrown
# times = [ 0.00001,
#           0.00010,
#           0.00100,
#           0.01000,
#           0.10000,
#           1.00000 ]
times = [ 0.0001,
          0.0010,
          0.0100,
          0.1000,
          1.0000,
          10.0000]
# Instantiate array for dumping timesteps
# logDump = np.zeros((len(times) - 1) * 9)
logDump = np.zeros(((len(times) - 1) * 90) + 1)

# Get number of timesteps in 1 tauBrownian
tStepPerTau = int(float(totTsteps) / (float(runFor) * float(tauBrown)))

# Little loop to give the desired values
count = 0
for i in range(0, len(times) - 1 ):
    # vals = np.arange(times[i], times[i + 1], times[i])
    vals = np.arange(times[i], times[i + 1], times[i] / 10)
    for j in range(0, len(vals)):
        logDump[count] = vals[j]
        count += 1

logDump *= tStepPerTau
logDump += brownEquil

# Remove .gsd files if they exist
try:
    os.remove(logName)
except OSError:
    pass

try:
    os.remove(gsdName)
except OSError:
    pass

def dump_spec(timestep):

    if timestep in logDump.astype(int):

        hoomd.dump.gsd(filename=logName,
                       period=None,
                       group=all,
                       overwrite=False,
                       time_step=timestep,
                       dynamic=['attribute', 'property', 'momentum'])

        # Need to close duplicate logfile processes
        proc = psutil.Process()
        of = proc.open_files()
        for i in xrange(len(of)):
            if of[i][0] == chkName:
                os.close(of[i][1])

    if timestep > logDump[-1]:
        analyzer.disable()

analyzer = hoomd.analyze.callback(callback = dump_spec, period = 1)
####################

hoomd.dump.gsd(gsdName,
               period=dumpFreq,
               group=all,
               overwrite=False,
               phase=-1,
               dynamic=['attribute', 'property', 'momentum'])

hoomd.run(totTsteps)
