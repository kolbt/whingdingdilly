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

isBallistic = str(sys.argv[1])  # ballistic or active simulation
pe = float(sys.argv[2])         # activity of active particle
lat = float(sys.argv[3])        # lattice spacing of HCP phase

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
inFile = "pe" + str(pe) + "_lattice" + str(lat) + ".gsd"
outFile = "pe" + str(pe) + "_lattice" + str(lat) + ".txt"
# File name for ballistic simulation
if isBallistic == "y":
    f = hoomd.open(name=inFile, mode='rb')
# File name for active simulation
else:
    inFile = "active_" + inFile
    outFile = "active_" + outFile
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
        hcpPos[iii] = snap.particles.position[:-1]      # get hcp particle positions
        actPos[iii] = snap.particles.position[-1]       # get active particle positions
        timesteps[iii] = snap.configuration.step        # get timestep

# Normalize timesteps
tauPerDT = computeTauPerTstep(epsilon=1.)
timesteps -= timesteps[0]
dt = timesteps[1] - timesteps[0]
# Brownian time
timesteps *= tauPerDT
# Rotational time units
#timesteps *= (3. * tauPerDT)
dtau = timesteps[1] - timesteps[0]

# Box dimensions
x_box = box_data[0]
y_box = box_data[1]
hx_box = x_box / 2.0
hy_box = y_box / 2.0
a_box = x_box * y_box

# Extent of HCP phase
hcpMin = min(hcpPos[0][:,0])
hcpMax = max(hcpPos[0][:,0])
actVelHCP = np.zeros((end), dtype=np.float64)
actVelFree = np.zeros((end), dtype=np.float64)
infil = 'n'

# Compute the MSD in the HCP phase
for i in range(start + 1, end):
    # Compute the distance traveled between timesteps (include periodicity)
    dx = abs(actPos[i][0] - actPos[i-1][0])
    dy = abs(actPos[i][1] - actPos[i-1][1])
    if dx > hx_box:
        dx = x_box - dx
    if dy > hy_box:
        dy = y_box - dy
        
    # Is/was the particle in the HCP phase
    if hcpMin < actPos[i][0] < hcpMax and hcpMin < actPos[i-1][0] < hcpMax:
        actVelHCP[i] = (np.sqrt((dx**2)+(dy**2))) / dtau
        infil='y'
    else:
        actVelFree[i] = (np.sqrt((dx**2)+(dy**2))) / dtau
  
## We can plot the free and dense phase active particle velocity
#plt.plot(timesteps, actVelFree, c='r', label='Free')
#plt.plot(timesteps, actVelHCP, c='b', ls='--', label='HCP')
#plt.xlabel(r'Time $(\tau_{B})$')
#plt.ylabel(r'Velocity $\left(\frac{d\sigma}{d\tau_{B}}\right)$')
#plt.legend(loc=1)
#plt.xlim(0, max(timesteps))
#plt.ylim(0, pe)
#plt.show()

# Compute and overlay the average
if np.any(actVelHCP):
    avgVelHCP = actVelHCP[actVelHCP!=0].mean()
else:
    avgVelHCP = 0.
avgVelFree = actVelFree[actVelFree!=0].mean()
print('Mean HCP velocity: {}').format(avgVelHCP)
print('Mean free velocity: {}').format(avgVelFree)

#print(actPos)            # all timesteps and coords
#print(actPos[0])         # timestep = 0
#print(actPos[0][0])      # timestep = 0, x-coord
#print(hcpPos)           # all timesteps and coords
#print(hcpPos[-1])       # last timestep
#print(hcpPos[-1][:,0])  # last timestep x coordinates

## Plot the data
#mysz = 15.
#x, y, z = zip(*actPos)
## Plot the HCP phase
#plt.scatter(hcpPos[-1][:,0], hcpPos[-1][:,1], c='#D3D3D3', s=mysz)
## Plot the time resolved active particle
#plt.scatter(x, y, c=(timesteps/timesteps[-1]), s=mysz, cmap='jet')
#plt.xlim(-hx_box, hx_box)
#plt.ylim(-hy_box, hy_box)
#ax = plt.gca()
#ax.set_aspect('equal')
#cbar = plt.colorbar(ticks=[], fraction=0.02675, pad=0.05, label='Time')
#plt.text(1.085, 0.015, r'$\tau_{initial}$', transform=ax.transAxes, fontsize=14)
#plt.text(1.085, 0.95, r'$\tau_{final}$', transform=ax.transAxes, fontsize=14)
#plt.show()

# Make text file for relevant quantities
f = open(outFile, 'w') # write file headings
f.write('Activity'.center(10) + ' ' +\
        'Lattice'.center(10) + ' ' +\
        'Infiltrate'.center(10) + ' ' +\
        'Velocity'.center(10) + '\n')
f.close()

f = open(outFile, 'a')
f.write('{0:.1f}'.format(pe).center(10) + ' ' +\
        '{0:.2f}'.format(lat).center(10) + ' ' + \
        str(infil).center(10) + ' ' +\
        '{0:.2f}'.format(avgVelHCP).center(10) + '\n')
f.close()
