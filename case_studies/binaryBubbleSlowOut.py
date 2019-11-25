'''
#                        This is an 80 character line                          #

These will only be binary simulations

'''
# Initial imports
import sys
import os

# Run locally
sys.path.append('/Users/kolbt/Desktop/compiled/hoomd-blue/build')
# Run on the cpu
#sys.path.append('/nas/longleaf/home/kolbt/programs/cpu-hoomd/hoomd-blue/build')
# Run on the gpu
sys.path.append('/nas/longleaf/home/kolbt/programs/hoomd_2.2.1/hoomd-blue/build')

import hoomd
from hoomd import md
from hoomd import dem
from hoomd import deprecated
from hoomd import data

import math
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
dumpPerBrownian = 500.                  # number of dumps per 1 tauB
simLength = 1.0 * tauBrown              # how long to run (in tauBrown)
totTsteps = int(simLength / dt)         # how many tsteps to run
numDumps = simLength * dumpPerBrownian  # total number of frames dumped
dumpFreq = totTsteps / numDumps         # normalized dump frequency
dumpFreq = int(dumpFreq)                # ensure this is an integer
seed = 71996                            # a random seed
seed2 = 2394                            # orientation seed
seed3 = 183                             # activity seed

# Some parameters:
N = 100000.                 # number of particles
xF = 0.5                    # fraction of each type
Ns = N * (1. - xF)          # number of slow particles
Nf = N - Ns                 # number of fast particles
pes = 50                    # slow activity
pef = 500                   # fast activity
phig = 0.04                 # area fraction of gas phase
phitot = 0.6                # total area fraciton
Ng = (phig / phitot) * N    # number of gas particles

def computeLat(activity):
    "Get lattice spacing based on effective diameter"
    m = -0.0957707250171
    b = 0.306774161185
    lattice = (activity**m) * (np.exp(b))
    return lattice
    
lats = computeLat(pes)
latf = computeLat(pef)

def areaType(Nx, latx):
    Ax = Nx * np.pi * 0.25 * (latx**2)
    return Ax
    
areas = areaType(Ns - Ng, lats)
areaf = areaType(Nf, latf)
areatot = areaf + areas
rtot = np.sqrt(areatot / np.pi)
# This depends on choice of interior species
rIn = np.sqrt(areaf / np.pi)
    
def computeDistance(x, y):
    return np.sqrt((x**2) + (y**2))
    
def interDist(x1, y1, x2, y2):
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
def orientToOrigin(x, y, act):
    "Using similar triangles to find sides"
    x *= -1
    y *= -1
    hypRatio = act / np.sqrt(x**2 + y**2)
    xAct = hypRatio * x
    yAct = hypRatio * y
    return xAct, yAct
    
# List of activities
peList = [ pef, pes ]
# List of ring radii
rList = [ 0, rIn, rtot ]
# Depth of alignment
rAlign = 3.
# List to store particle positions and types
pos = []
typ = []
rOrient = []
# z-value for simulation initialization
z = 0.5

for i in xrange(len(peList)):
    rMin = rList[i]             # starting distance for particle placement
    rMax = rList[i + 1]         # maximum distance for particle placement
    lat = computeLat(peList[i]) # activity-dependent lattice spacing
    ver = np.sqrt(0.75) * lat   # vertical shift between lattice rows
    hor = lat / 2.0             # horizontal shift between lattice rows
    
    x = 0
    y = 0
    shift = 0
    while y < rMax:
        r = computeDistance(x, y)
        # Check if x-position is large enough
        if r < rMin:
            x += lat
            continue
            
        # Check if x-position is too large
        if r >= (rMax - (lat/2.)):
            y += ver
            shift += 1
            if shift % 2:
                x = hor
            else:
                x = 0
            continue
        
        # Whether or not particle is oriented
        if r > (rMax - rAlign):
            # Aligned
            rOrient.append(1)
        else:
            # Random
            rOrient.append(0)
        
        # If the loop makes it this far, append
        pos.append((x, y, z))
        typ.append(i)
        if x != 0 and y != 0:
            # Mirror positions, alignment and type
            pos.append((-x, y, z))
            pos.append((-x, -y, z))
            pos.append((x, -y, z))
            rOrient.append(rOrient[-1])
            rOrient.append(rOrient[-1])
            rOrient.append(rOrient[-1])
            typ.append(i)
            typ.append(i)
            typ.append(i)
        # y must be zero
        elif x != 0:
            pos.append((-x, y, z))
            rOrient.append(rOrient[-1])
            typ.append(i)
        # x must be zero
        elif y!= 0:
            pos.append((x, -y, z))
            rOrient.append(rOrient[-1])
            typ.append(i)
         
        # Increment counter
        x += lat

# Set this according to phiTotal
areaParts = N * np.pi * (0.25)
abox = (areaParts / 0.6)
lbox = np.sqrt(abox)
hbox = lbox / 2.
tooClose = 0.8

# Make a mesh for random particle placement
def getNBins(length, minSz=(2**(1./6.))):
    "Given box size, return number of bins"
    initGuess = int(length) + 1
    nBins = initGuess
    # This loop only exits on function return
    while True:
        if length / nBins > minSz:
            return nBins
        else:
            nBins -= 1
# Round up size of bins to account for floating point inaccuracy
def roundUp(n, decimals=0):
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier
# Compute mesh
r_cut = 2**(1./6.)
nBins = (getNBins(lbox, r_cut))
sizeBin = roundUp((lbox / nBins), 6)

# Place particles in gas phase
count = 0
gaspos = []
binParts = [[[] for b in range(nBins)] for a in range(nBins)]
while count < Ng:
    place = 1
    # Generate random position
    gasx = (np.random.rand() - 0.5) * lbox
    gasy = (np.random.rand() - 0.5) * lbox
    r = computeDistance(gasx, gasy)
    
    # Is this an HCP bin?
    if r <= (rList[-1] + (tooClose / 2.)):
        continue
    
    # Are any gas particles too close?
    tmpx = gasx + hbox
    tmpy = gasy + hbox
    indx = int(tmpx / sizeBin)
    indy = int(tmpy / sizeBin)
    # Get index of surrounding bins
    lbin = indx - 1  # index of left bins
    rbin = indx + 1  # index of right bins
    bbin = indy - 1  # index of bottom bins
    tbin = indy + 1  # index of top bins
    if rbin == nBins:
        rbin -= nBins  # adjust if wrapped
    if tbin == nBins:
        tbin -= nBins  # adjust if wrapped
    hlist = [lbin, indx, rbin]  # list of horizontal bin indices
    vlist = [bbin, indy, tbin]  # list of vertical bin indices

    # Loop through all bins
    for h in range(0, len(hlist)):
        for v in range(0, len(vlist)):
            # Take care of periodic wrapping for position
            wrapX = 0.0
            wrapY = 0.0
            if h == 0 and hlist[h] == -1:
                wrapX -= lbox
            if h == 2 and hlist[h] == 0:
                wrapX += lbox
            if v == 0 and vlist[v] == -1:
                wrapY -= lbox
            if v == 2 and vlist[v] == 0:
                wrapY += lbox
            # Compute distance between particles
            if binParts[hlist[h]][vlist[v]]:
                for b in range(0, len(binParts[hlist[h]][vlist[v]])):
                    # Get index of nearby particle
                    ref = binParts[hlist[h]][vlist[v]][b]
                    r = interDist(gasx, gasy,
                                  gaspos[ref][0] + wrapX,
                                  gaspos[ref][1] + wrapY)
                    # Round to 4 decimal places
                    r = round(r, 4)
                    # If too close, generate new position
                    if r <= tooClose:
                        place = 0
                        break
            if place == 0:
                break
        if place == 0:
            break
            
    # Is it safe to append the particle?
    if place == 1:
        binParts[indx][indy].append(count)
        gaspos.append((gasx, gasy, z))
        rOrient.append(0)       # not oriented
        typ.append(1)           # final particle type, same as outer ring
        count += 1              # increment count

# Get each coordinate in a list
pos = pos + gaspos
x, y, z = zip(*pos)

## Plot as scatter
#cs = np.divide(typ, float(len(peList)))
##cs = rOrient
#plt.scatter(x, y, s=1., c=cs, cmap='jet', edgecolors='none')
#ax = plt.gca()
#ax.set_aspect('equal')
#plt.show()

partNum = len(pos)
# Get the number of types
uniqueTyp = []
for i in typ:
    if i not in uniqueTyp:
        uniqueTyp.append(i)
# Get the number of each type
particles = [ 0 for x in xrange(len(uniqueTyp)) ]
for i in xrange(len(uniqueTyp)):
    for j in typ:
        if uniqueTyp[i] == j:
            particles[i] += 1
# Convert types to letter values
unique_char_types = []
for i in uniqueTyp:
    unique_char_types.append( chr(ord('@') + i+1) )
char_types = []
for i in typ:
    char_types.append( chr(ord('@') + i+1) )

# Get a list of activities for all particles
pe = []
for i in typ:
    pe.append(peList[i])

# Now we make the system in hoomd
hoomd.context.initialize()
# A small shift to help with the periodic box
snap = hoomd.data.make_snapshot(N = partNum,
                                box = hoomd.data.boxdim(Lx=lbox,
                                                        Ly=lbox,
                                                        dimensions=2),
                                particle_types = unique_char_types)

# Set positions/types for all particles
snap.particles.position[:] = pos[:]
snap.particles.typeid[:] = typ[:]
snap.particles.types[:] = char_types[:]

# Initialize the system
system = hoomd.init.read_snapshot(snap)
all = hoomd.group.all()
groups = []
for i in unique_char_types:
    groups.append(hoomd.group.type(type=i))

# Set particle potentials
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=2**(1/6), nlist=nl)
lj.set_params(mode='shift')
for i in xrange(len(unique_char_types)):
    for j in range(i, len(unique_char_types)):
        lj.pair_coeff.set(unique_char_types[i],
                          unique_char_types[j],
                          epsilon=eps, sigma=sigma)

# Set activity of each group
np.random.seed(seed2)                           # seed for random orientations
angle = np.random.rand(partNum) * 2 * np.pi     # random particle orientation
activity = []
for i in xrange(partNum):
    if rOrient[i] == 0:
        x = (np.cos(angle[i])) * pe[i]
        y = (np.sin(angle[i])) * pe[i]
    else:
        x, y = orientToOrigin(pos[i][0], pos[i][1], pe[i])
    z = 0.
    tuple = (x, y, z)
    activity.append(tuple)
# Implement the activities in hoomd
hoomd.md.force.active(group=all,
                      seed=seed3,
                      f_lst=activity,
                      rotation_diff=D_r,
                      orientation_link=False,
                      orientation_reverse_link=True)

# Brownian integration
hoomd.md.integrate.mode_standard(dt=dt)
bd = hoomd.md.integrate.brownian(group=all, kT=kT, seed=seed)

# Name the file from parameters
out = "outside_slow_pe"
for i in peList:
    out += str(int(i))
    out += "_"
out += "r"
for i in range(1, len(rList)):
    out += str(int(rList[i]))
    out += "_"
out += "rAlign_" + str(rAlign) + ".gsd"

# Write dump
hoomd.dump.gsd(out,
               period=dumpFreq,
               group=all,
               overwrite=True,
               phase=-1,
               dynamic=['attribute', 'property', 'momentum'])

# Run
hoomd.run(totTsteps)
