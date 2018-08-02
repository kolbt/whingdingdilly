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
from scipy.interpolate import griddata
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
import math

inFile = sys.argv[2]
gsdPath = sys.argv[3]
myFrame = sys.argv[4]

sys.path.append(gsdPath)
import gsd
from gsd import hoomd

sigma=1.0
# The old model
epsilon = 1.0

# The new model
def computeEps(activity):
    "Given particle activity, output repulsion well depth"
    epsilon = activity * sigma / 24.0
    return epsilon

# Function that'll grab my parameters from the filenames
def getFromTxt(fname, first, last):
    """Takes a string, text before and after desired text, outs text between"""
    start = fname.index( first ) + len( first )
    end = fname.index( last, start )
    myTxt = fname[start:end]
    return float(myTxt)

# Computes distance
def getDistance(point1, point2x, point2y):
    """"Find the distance between two points"""
    distance = np.sqrt((point2x - point1[0])**2 + (point2y - point1[1])**2)
    return distance

# Computes force given distance with LJ potential
def computeLJForce(r, eps_in):
    """Given a distance, computes the force"""
    forceLJ = 4*eps_in*((12*(sigma**12)*(r**-13))-(6*(sigma**12)*(r**-7)))
    return forceLJ

# If epsilon is not in filename
peA = getFromTxt(inFile, "pa", "_pb")
peB = getFromTxt(inFile, "pb", "_xa")
xA = getFromTxt(inFile, "xa", ".gsd")

# Only if epsilon is in the filename
if type(xA) is str:
    peA = getFromTxt(inFile, "pa", "_pb")
    peB = getFromTxt(inFile, "pb", "_xa")
    xA = getFromTxt(inFile, "xa", "_ep")
    ep = getFromTxt(inFile, "ep", ".gsd")

try:
    peR = float(peA) / float(peB)       # Compute activity ratio
except:
    peR = 0

epsilonA = computeEps(peA)
epsilonB = computeEps(peB)
epsHS = epsilonA if epsilonA > epsilonB else epsilonB
epsHS = 1
partFracA = xA / 100.0    # Particle fraction

# There are two options for the infile name
fname = "pa" + str(peA) +\
"_pb" + str(peB) +\
"_xa" + str(xA) +\
"_ep" + str(ep) +\
".gsd"

fname = "pa" + str(peA) +\
"_pb" + str(peB) +\
"_xa" + str(xA) +\
".gsd"

# need to import hoomd
import hoomd
hoomd.init.read_gsd(fname, myFrame)