'''
#                           This is an 80 character line                       #
This file will restart a simulation at a time when a cluster is nucleating. This
will allow us to take a look at composition of the dense phase at short
timescales.

What this file does:
    1. Load in a gsd file at the nucleation time
    2. Run the simulation with the same parameters
    3. Process the resultant .gsd to get fine timescale cluster formation
'''

# Imports and loading the .gsd file
import sys
import numpy as np

hoomdPath = sys.argv[1]
print(hoomdPath)
gsdPath = sys.argv[2]
peA = int(sys.argv[3])
peB = int(sys.argv[4])
xA = int(sys.argv[5])
ep = float(sys.argv[6])
seed1 = int(sys.argv[7])
seed2 = int(sys.argv[8])
seed3 = int(sys.argv[9])
seed4 = int(sys.argv[10])
seed5 = int(sys.argv[11])
myFrame = int(sys.argv[12])

# sys.path.append(gsdPath)
sys.path.append(hoomdPath)

# import gsd
# from gsd import hoomd
import hoomd
from hoomd import md

# Get necessary parameters
sigma = 1.0
threeEtaPiSigma = 1.0
kT = 1.0
D_t = kT / threeEtaPiSigma      # translational diffusion constant
D_r = (3.0 * D_t) / (sigma**2)  # rotational diffusion constant
tauBrown = (sigma**2) / D_t     # brownian time scale (invariant)
def computeTauLJ(epsilon):
    "Given epsilon, compute lennard-jones time unit"
    tauLJ = ((sigma**2) * threeEtaPiSigma) / epsilon
    return tauLJ
tauLJ = computeTauLJ(ep)
dt = 0.00001 * tauLJ

fname = 'pa' + str(peA) +\
'_pb' + str(peB) +\
'_xa' + str(xA) +\
'.gsd'
nuc = 'nuc_' + fname
nucFreq = 100

# BEFORE READING IN GSD TO RUN
# Read in and put orientation data into array (for myFrame)
# Why isn't this default behavior?

hoomd.context.initialize()
system = hoomd.init.read_gsd(filename=fname, frame=myFrame)

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
lj.pair_coeff.set('A', 'A', epsilon=ep, sigma=1.0)
lj.pair_coeff.set('A', 'B', epsilon=ep, sigma=1.0)
lj.pair_coeff.set('B', 'B', epsilon=ep, sigma=1.0)

activity_a = [(peA,0,0) for i in range(Na)]
activity_b = [(peB,0,0) for i in range(Nb)]

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

# General integration parameters
hoomd.md.integrate.mode_standard(dt=dt)
hoomd.md.integrate.brownian(group=all, kT=kT, seed=seed2)

hoomd.dump.gsd(nuc,
               period=nucFreq,
               group=all,
               overwrite=False,
               phase=-1,
               dynamic=['attribute', 'property', 'momentum'])

# Run the simulation
# hoomd.run(nucFreq[-1])
hoomd.run(10000)
