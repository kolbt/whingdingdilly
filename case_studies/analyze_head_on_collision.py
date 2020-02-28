'''
#                           This is an 80 character line                       #
What is the equilibrium separation of a head on collision?
'''

import sys

gsd_path = '/nas/longleaf/home/kolbt/programs/gsd/build'
sys.path.append(gsd_path)       # ensure gsd is in your python path
gsd_path = '/Users/kolbt/Desktop/compiled/gsd/build'
sys.path.append(gsd_path)       # ensure gsd is in your python path

import gsd
from gsd import hoomd
from gsd import pygsd

import freud
from freud import parallel
from freud import box
from freud import density
from freud import cluster

import numpy as np
import math
from scipy import stats

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
    
# Base filename
inFile = 'head_on_collision.gsd'
f = hoomd.open(name=inFile, mode='rb')

dump = int(f.__len__())         # get number of timesteps dumped
start = 0                       # gives first frame to read
end = dump                      # gives last frame to read
box_data = np.zeros((1), dtype=np.ndarray)        # box dimensions

# Access and read .gsd data
with hoomd.open(name=inFile, mode='rb') as t:
    snap = t[-1]
    first_tstep = snap.configuration.step
    box_data = snap.configuration.box
    # Box dimensions
    x_box = box_data[0]
    y_box = box_data[1]
    hx_box = x_box / 2.0
    hy_box = y_box / 2.0
    a_box = x_box * y_box
    
    # Get the positions
    pos = snap.particles.position[:]
    delta = []
    for i in range(1, len(pos)):
        delta.append(pos[i][0] - pos[i-1][0])
    print("Separation distances: ", delta)

# Now let's see what force this corresponds to
def ljForce(r, eps=1., sigma=1.):
    div = (sigma / r)
    dU = (24. * eps / r) * (2*(div**12) - (div**6))
    return dU
    
ljF = []
for i in range(0, len(delta)):
    ljF.append(ljForce(delta[i]))
print("The Lennard-Jones forces: ", ljF)
