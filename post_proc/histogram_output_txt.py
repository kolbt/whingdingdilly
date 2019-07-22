'''
#                           This is an 80 character line                       #
This file is intended to look at the average distance between particles in the
dense phase (delta). This should directly result from the forces applied on
particles in the dense phase (that's how potentials work), but, we want to make
sure we are using an appropriate method of implementing the well depth
(epsilon).

What this file does:
    1. Read in your .gsd file
    2. Mesh the system
    3. Compute distance between 'close' particles (use r_cut=1.122)
    4. Store as lists of three types: AA, AB, BB
    5. Plot lists as population histograms
'''

# Imports and loading the .gsd file
import sys

peA = int(sys.argv[1])                     # activity A
peB = int(sys.argv[2])                     # activity B
part_perc_a = int(sys.argv[3])              # percentage A particles
part_frac_a = float(part_perc_a) / 100.0    # fraction A particles
hoomd_path = str(sys.argv[4])               # local path to hoomd-blue
gsd_path = str(sys.argv[5])                 # local path to gsd
try:
    inEps = int(sys.argv[6])                    # epsilon (if given)
except:
    inEps = 1

try:
    phi = int(sys.argv[7])
except:
    phi = 60
sys.path.append(hoomd_path)     # ensure hoomd is in your python path
sys.path.append(gsd_path)       # ensure gsd is in your python path

import hoomd
from hoomd import md
from hoomd import deprecated

import gsd
from gsd import hoomd
from gsd import pygsd

import freud
from freud import parallel
from freud import box
from freud import density
from freud import cluster

import numpy as np
from scipy.interpolate import griddata
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
# sns.set(color_codes=True)

import math

# Set some constants
kT = 1.0                        # temperature
threeEtaPiSigma = 1.0           # drag coefficient
sigma = 1.0                     # particle diameter
D_t = kT / threeEtaPiSigma      # translational diffusion constant
D_r = (3.0 * D_t) / (sigma**2)  # rotational diffusion constant
tauBrown = (sigma**2) / D_t     # brownian time scale (invariant)

# Any functions I'll need go here
def getDistance(point1, point2x, point2y):
    "Find the distance between two points"
    distance = np.sqrt((point2x - point1[0])**2 + (point2y - point1[1])**2)
    return distance

def computeVel(activity):
    "Given particle activity, output intrinsic swim speed"
    velocity = (activity * sigma) / (3 * (1/D_r))
    return velocity

def computeActiveForce(velocity):
    "Given particle activity, output repulsion well depth"
    activeForce = velocity * threeEtaPiSigma
    return activeForce

def computeEps(activeForce):
    "Given particle activity, output repulsion well depth"
    epsilon = activeForce * sigma / 24.0
    return epsilon

def computeTauLJ(epsilon):
    "Given epsilon, compute lennard-jones time unit"
    tauLJ = ((sigma**2) * threeEtaPiSigma) / epsilon
    return tauLJ

def getLJForce(r, epsilon):
    "Given a distance compute the force being experienced"
    ljForce = 48*epsilon*(((sigma**12)/(r**13))-((sigma**6)/(r**7)))
    return ljForce

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

eps = (epsA if (epsA >= epsB) else epsB)    # use the larger epsilon
Fp = (FpA if (FpA >= FpB) else FpB)         # use the larger active force

try:
    in_file = "pa" + str(peA) + \
              "_pb" + str(peB) + \
              "_xa" + str(part_perc_a) + \
              ".gsd"
    # File base
    b_file = "pa" + str(peA) + \
             "_pb" + str(peB) + \
             "_xa" + str(part_perc_a)
    f = hoomd.open(name=in_file, mode='rb') # open gsd file with hoomd
# Use this if phi is not in filename
except:
    in_file = "pa"+str(peA)+\
              "_pb"+str(peB)+\
              "_xa"+str(part_perc_a)+\
              "_ep"+str(inEps)+\
              ".gsd"
    # File base
    b_file = "pa" + str(peA) + \
             "_pb" + str(peB) + \
             "_xa" + str(part_perc_a) + \
             "_ep" + str(inEps)
    f = hoomd.open(name=in_file, mode='rb')  # open gsd file with hoomd
# Use this if phi is in filename
#except:
#    in_file = "pa"+str(peA)+\
#              "_pb"+str(peB)+\
#              "_xa"+str(part_perc_a)+\
#              "_ep"+str(inEps)+\
#              "_phi"+str(phi)+\
#              ".gsd"
#    # File base
#    b_file = "pa" + str(peA) + \
#             "_pb" + str(peB) + \
#             "_xa" + str(part_perc_a) + \
#             "_ep" + str(inEps) + \
#             "_phi" + str(phi)
#    f = hoomd.open(name=in_file, mode='rb')  # open gsd file with hoomd

dumps = f.__len__()                     # get number of timesteps dumped

start = 0       # gives first frame to read
end = dumps     # gives last frame to read

# Run for last 20 tsteps?
start = dumps - 2
#end = 401

positions = np.zeros((end), dtype=np.ndarray)       # array of positions
types = np.zeros((end), dtype=np.ndarray)           # particle types
box_data = np.zeros((1), dtype=np.ndarray)          # box dimensions
timesteps = np.zeros((end), dtype=np.float64)       # timesteps

# Get relevant data from .gsd file
with hoomd.open(name=in_file, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    for iii in range(start, end):
        snap = t[iii]                               # snapshot of frame
        types[iii] = snap.particles.typeid          # get types
        positions[iii] = snap.particles.position    # get positions
        timesteps[iii] = snap.configuration.step    # get timestep

timesteps -= timesteps[0]       # get rid of brownian run time

# Get number of each type of particle
part_num = len(types[start])
part_A = int(part_num * part_frac_a)
part_B = part_num - part_A

# Feed data into freud analysis software
l_box = box_data[0]
h_box = l_box / 2.0
a_box = l_box * l_box
f_box = box.Box(Lx = l_box, Ly = l_box, is2D = True)    # make freud box

# Make the mesh
#nBins = 200
#sizeBin = l_box / nBins
r_cut = 1.122
sizeBin = r_cut
nBins = int(l_box / sizeBin)
nBins += 1                      # account for integer rounding

#print(sizeBin)
# Instantiate lists here to sum temporally
ALL = []
AA = []
BB = []
AB = []

for iii in range(start, end):
    
    # Mesh array
    binParts = [ [ [] for b in range(nBins) ] for a in range(nBins) ]

    # Easier accessors
    pos = positions[iii]
    pos = np.delete(pos, 2, 1)
    typ = types[iii]
    
    # Instantiate pair-wise lists (no time-average)
#    ALL = []
#    AA = []
#    BB = []
#    AB = []

    # Put particles in their respective bins
    for jjj in range(0, part_num):
        # Get mesh indices
        tmp_posX = pos[jjj][0] + h_box
        tmp_posY = pos[jjj][1] + h_box
        x_ind = int(tmp_posX / sizeBin)
        y_ind = int(tmp_posY / sizeBin)
        # Append particle id to appropriate bin
        binParts[x_ind][y_ind].append(jjj)

    # Compute distance, each pair will be counted twice
    for jjj in range(0, part_num):
        # Get mesh indices
        tmp_posX = pos[jjj][0] + h_box
        tmp_posY = pos[jjj][1] + h_box
        x_ind = int(tmp_posX / sizeBin)
        y_ind = int(tmp_posY / sizeBin)
        # Get index of surrounding bins
        l_bin = x_ind - 1   # index of left bins
        r_bin = x_ind + 1   # index of right bins
        b_bin = y_ind - 1   # index of bottom bins
        t_bin = y_ind + 1   # index of top bins
        if r_bin == nBins:
            r_bin -= nBins  # adjust if wrapped
        if t_bin == nBins:
            t_bin -= nBins  # adjust if wrapped
        h_list = [l_bin, x_ind, r_bin]  # list of horizontal bin indices
        v_list = [b_bin, y_ind, t_bin]  # list of vertical bin indices

        # Loop through all bins
        for h in range(0, len(h_list)):
            for v in range(0, len(v_list)):
                # Take care of periodic wrapping for position
                wrapX = 0.0
                wrapY = 0.0
                if h == 0 and h_list[h] == -1:
                    wrapX -= l_box
                if h == 2 and h_list[h] == 0:
                    wrapX += l_box
                if v == 0 and v_list[v] == -1:
                    wrapY -= l_box
                if v == 2 and v_list[v] == 0:
                    wrapY += l_box
                # Compute distance between particles
                for b in range(0, len(binParts[h_list[h]][v_list[v]])):
                    ref = binParts[h_list[h]][v_list[v]][b]
                    r = getDistance(pos[jjj],
                                    pos[ref][0] + wrapX,
                                    pos[ref][1] + wrapY)
                    r = round(r, 4) # round value to 4 decimal places
                    
                    # If LJ potential is on, store into a list (omit self)
                    if 0.1 < r <= r_cut:
                        ALL.append(format(r, '.4f'))                           # All particles
                        if typ[jjj] == 0 and typ[ref] == 0:     # AA distance
                            AA.append(format(r, '.4f'))
                        elif typ[jjj] == 1 and typ[ref] == 1:   # BB distance
                            BB.append(format(r, '.4f'))
                        else:                                   # AB distance
                            AB.append(format(r, '.4f'))

# Create the headers for my text output
myOut = b_file + '.txt'
f = open(myOut, 'w') # write file headings
f.write('ALL'.center(5) + '\t' +\
        'AA'.center(5) + '\t' +\
        'AB'.center(5) + '\t' +\
        'BB'.center(5) + '\n')
f.close()
# This will write my lists to a text file
import itertools as it
lis=[ALL, AA, AB, BB]
with open(myOut, 'a') as f:
    for x in it.izip_longest(*lis, fillvalue=''):
        f.write("{0}\t{1}\t{2}\t{3}\n".format(*x))
