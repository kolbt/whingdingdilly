'''
#                           This is an 80 character line                       #
This file is intended to track the movement of the ballistic or active particle
over time.

We want to know:
    -if reorientation allows for enhanced mobility in the HCP phase
    -the critical activity for intrusion into HCP phase
    -the critical lattice spacing
    -relationship between mobility and activity
'''

import sys

pe = float(sys.argv[1])         # activity of active particle
lat = float(sys.argv[2])        # lattice spacing of HCP phase

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

def computeA(diameter):
    """Computes area of circle"""
    radius = diameter / 2.0
    return np.pi * (radius**2)

def getDistance(point1, point2x, point2y):
    """Find the distance between two points"""
    distance = np.sqrt((point2x - point1[0])**2 + (point2y - point1[1])**2)
    return distance
    
def computeTauPerTstep(epsilon, mindt=0.000001):
    """Read in epsilon, output tauBrownian per timestep"""
    kBT = 1.0
    tstepPerTau = float(epsilon / (kBT * mindt))
    return 1. / tstepPerTau
    
# Base filename
inFile = "active_in_bulk_pe" + str(pe) + "_lattice" + str(lat) + ".gsd"
outFile = "active_in_bulk_pe" + str(pe) + "_lattice" + str(lat) + ".txt"
f = hoomd.open(name=inFile, mode='rb')

dump = int(f.__len__())         # get number of timesteps dumped
start = 0                       # gives first frame to read
end = dump                      # gives last frame to read

hcpPos = np.zeros((end), dtype=np.ndarray)
actPos = np.zeros((end), dtype=np.ndarray)          # array of active positions
types = np.zeros((end), dtype=np.ndarray)           # particle types
box_data = np.zeros((1), dtype=np.ndarray)          # box dimensions
timesteps = np.zeros((end), dtype=np.float64)       # timesteps

# Access and read .gsd data
with hoomd.open(name=inFile, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    for iii in range(start, end):
        snap = t[iii]                                   # snapshot of frame
        actPos[iii] = snap.particles.position[0]        # get active particle positions
        timesteps[iii] = snap.configuration.step        # get timestep

# Normalize timesteps
tauPerDT = computeTauPerTstep(epsilon=1.)
dtInTauB = int(1. / (tauPerDT))
looping = 0
while looping == 0:
    if dtInTauB % 3 != 0 and tInTauB % 2 != 0:
        dtInTauB += 1
    else:
        looping = 1
dtInTauR = dtInTauB / 3
dtHalfTauR = dtInTauR / 2
timesteps -= timesteps[0]
dt = timesteps[-1] - timesteps[-2]
# Brownian time
timesteps *= tauPerDT
# Rotational time units
#timesteps *= (3. * tauPerDT)
dtau = timesteps[-1] - timesteps[-2]

# Box dimensions
x_box = box_data[0]
y_box = box_data[1]
hx_box = x_box / 2.0
hy_box = y_box / 2.0
a_box = x_box * y_box

# Extent of HCP phase
#instV = np.zeros((end), dtype=np.float64)
instV = []
halfV = []
taurV = []
taubV = []
# Compute the MSD in the HCP phase
for i in range(start + 1, end):
    # Compute the distance traveled between timesteps (include periodicity)
    dx = abs(actPos[i][0] - actPos[i-1][0])
    dy = abs(actPos[i][1] - actPos[i-1][1])
    if dx > hx_box:
        dx = x_box - dx
    if dy > hy_box:
        dy = y_box - dy
    instV.append(np.sqrt((dx**2)+(dy**2))) / dtau
    
    # Is this a multiple of a larger timescale (0.5 * tauR)?
    if i % dtHalfTauR == 0:
        print("At half tauR")
        dx = abs(actPos[i][0] - actPos[i-dtHalfTauR][0])
        dy = abs(actPos[i][1] - actPos[i-dtHalfTauR][1])
        if dx > hx_box:
            dx = x_box - dx
        if dy > hy_box:
            dy = y_box - dy
        halfV.append(np.sqrt((dx**2)+(dy**2))) / (dtau * dtHalfTauR)
        # Is this a multiple of tauR
        if i % dtInTauR == 0:
            print("At tauR")
            dx = abs(actPos[i][0] - actPos[i-dtInTauR][0])
            dy = abs(actPos[i][1] - actPos[i-dtInTauR][1])
            if dx > hx_box:
                dx = x_box - dx
            if dy > hy_box:
                dy = y_box - dy
            taurV.append(np.sqrt((dx**2)+(dy**2))) / (dtau * dtInTauR)
            # Is this a multiple of tauB
            if i % dtInTauB == 0:
                print("At tauB")
                dx = abs(actPos[i][0] - actPos[i-dtInTauB][0])
                dy = abs(actPos[i][1] - actPos[i-dtInTauB][1])
                if dx > hx_box:
                    dx = x_box - dx
                if dy > hy_box:
                    dy = y_box - dy
                taubV.append(np.sqrt((dx**2)+(dy**2))) / (dtau * dtInTauB)
                
# We can plot the free and dense phase active particle velocity
plt.plot(timesteps[start+1, end], instV, c='b', ls='--', label=r'$\nu_{inst.}$')
plt.plot(timesteps[start+1, end], halfV, c='r', ls='--', label=r'$\nu_{\tau_{r}/2}$')
plt.plot(timesteps[start+1, end], taurV, c='g', ls='--', label=r'$\nu_{\tau_{r}}$')
plt.plot(timesteps[start+1, end], taubV, c='k', ls='--', label=r'$\nu_{\tau_{B}}$')
plt.xlabel(r'Time $(\tau_{B})$')
plt.ylabel(r'Velocity $\left(\frac{d\sigma}{d\tau_{B}}\right)$')
plt.legend(loc=1)
plt.xlim(0, max(timesteps))
plt.ylim(0, pe)
plt.show()

# Compute and overlay the average
avgInstV = mean(instV)
avgHalfV = mean(halfV)
avgTaurV = mean(taurV)
avgTaubV = mean(taubV)

# Plot the data
mysz = 15.
x, y, z = zip(*actPos[start:end])
# Plot the time resolved active particle
plt.scatter(x, y, c=(timesteps[start:end]/max(timesteps)), s=mysz, cmap='jet')
plt.xlim(-hx_box, hx_box)
plt.ylim(-hy_box, hy_box)
ax = plt.gca()
ax.set_aspect('equal')
cbar = plt.colorbar(ticks=[], fraction=0.02675, pad=0.05, label='Time')
plt.text(1.085, 0.015, r'$\tau_{initial}$', transform=ax.transAxes, fontsize=14)
plt.text(1.085, 0.95, r'$\tau_{final}$', transform=ax.transAxes, fontsize=14)
plt.show()

# Make text file for relevant quantities
f = open(outFile, 'w')
f.write('Activity'.center(15) + ' ' +\
        'Lattice'.center(15) + ' ' +\
        'Instant_v'.center(15) + ' ' +\
        'half_taur_v'.center(15) + ' ' +\
        'taur_v'.center(15) + ' ' +\
        'taub_v'.center(15) + '\n')
f.close()

# Populate text file
f = open(outFile, 'a')
f.write('{0:.1f}'.format(pe).center(15) + ' ' +\
        '{0:.2f}'.format(lat).center(15) + ' ' + \
        '{0:.4f}'.format(avgInstV).center(15) + ' ' + \
        '{0:.4f}'.format(avgHalfV).center(15) + ' ' + \
        '{0:.4f}'.format(avgTaurV).center(15) + ' ' + \
        '{0:.4f}'.format(avgTaubV).center(15) + '\n')
f.close()
