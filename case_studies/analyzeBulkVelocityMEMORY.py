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
box_data = np.zeros((1), dtype=np.ndarray)        # box dimensions

# Normalization of timesteps
tauPerDT = computeTauPerTstep(epsilon=1.)
# Compute smallest unit to enforce larger timeaverages are multiples of it
dtHalfTauR = int(1. / (6. * tauPerDT))
dtInTauR = dtHalfTauR * 2
dtInTauB = dtHalfTauR * 6

# Timescales to output velocity
instV = 0.
instC = 0
instXY = []
halfV = 0.
halfC = 0
halfXY = []
taurV = 0.
taurC = 0
taurXY = []
taubV = 0.
taubC = 0
taubXY = []

# Access and read .gsd data
with hoomd.open(name=inFile, mode='rb') as t:
    snap = t[0]
    first_tstep = snap.configuration.step
    box_data = snap.configuration.box
    # Box dimensions
    x_box = box_data[0]
    y_box = box_data[1]
    hx_box = x_box / 2.0
    hy_box = y_box / 2.0
    a_box = x_box * y_box
    # Put the first entry in each array
    instXY.append(snap.particles.position[0])
    halfXY.append(snap.particles.position[0])
    taurXY.append(snap.particles.position[0])
    taubXY.append(snap.particles.position[0])
    # Need to get dt (so I don't have to check an if statement)
    snap = t[1]
    second_tstep = snap.configuration.step
    dt = second_tstep - first_tstep
    dtau = dt * tauPerDT
    for i in range(start + 1, end):
        snap = t[i]                                 # snapshot of frame
        instXY.append(snap.particles.position[0])   # get active particle positions
        instXY = instXY[-2:]                        # you only want to keep two entries
        timestep = snap.configuration.step          # get timestep
        timestep -= first_tstep                     # adjust timestep
        # Brownian time
        timestep *= tauPerDT

        # Compute the distance traveled between timesteps (include periodicity)
        dx = abs(instXY[-1][0] - instXY[-2][0])
        dy = abs(instXY[-1][1] - instXY[-2][1])
        if dx > hx_box:
            dx = x_box - dx
        if dy > hy_box:
            dy = y_box - dy
        instV += ( (np.sqrt((dx**2)+(dy**2))) / dtau )
        instC += 1
        
        # Is this a multiple of a larger timescale (0.5 * tauR)?
        if i % dtHalfTauR == 0:
            print("At half tauR")
            halfXY.append(snap.particles.position[0])
            halfXY = halfXY[-2:]
            dx = abs(halfXY[-1][0] - halfXY[-2][0])
            dy = abs(halfXY[-1][1] - halfXY[-2][1])
            if dx > hx_box:
                dx = x_box - dx
            if dy > hy_box:
                dy = y_box - dy
            # dtau gives the answer in units of dr/dtau_B (dtHalfTauR gives # of timesteps)
            halfV += ( (np.sqrt((dx**2)+(dy**2))) / (dtau * dtHalfTauR) )
            halfC += 1

            # Is this a multiple of tauR
            if i % dtInTauR == 0:
                print("At tauR")
                taurXY.append(snap.particles.position[0])
                taurXY = taurXY[-2:]
                dx = abs(taurXY[-1][0] - taurXY[-2][0])
                dy = abs(taurXY[-1][1] - taurXY[-2][1])
                if dx > hx_box:
                    dx = x_box - dx
                if dy > hy_box:
                    dy = y_box - dy
                taurV += ( (np.sqrt((dx**2)+(dy**2))) / (dtau * dtInTauR) )
                taurC += 1

                # Is this a multiple of tauB
                if i % dtInTauB == 0:
                    print("At tauB")
                    taubXY.append(snap.particles.position[0])
                    taubXY = taubXY[-2:]
                    dx = abs(taubXY[-1][0] - taubXY[-2][0])
                    dy = abs(taubXY[-1][1] - taubXY[-2][1])
                    if dx > hx_box:
                        dx = x_box - dx
                    if dy > hy_box:
                        dy = y_box - dy
                    taubV += ( (np.sqrt((dx**2)+(dy**2))) / (dtau * dtInTauB) )
                    taubC += 1
                
## We can plot the free and dense phase active particle velocity
#plt.plot(timesteps[start+1:end], instV, c='b', ls='--', label=r'$\nu_{inst.}$')
##plt.plot(timesteps[start+1:end], halfV, c='r', ls='--', label=r'$\nu_{\tau_{r}/2}$')
##plt.plot(timesteps[start+1:end], taurV, c='g', ls='--', label=r'$\nu_{\tau_{r}}$')
##plt.plot(timesteps[start+1:end], taubV, c='k', ls='--', label=r'$\nu_{\tau_{B}}$')
#plt.xlabel(r'Time $(\tau_{B})$')
#plt.ylabel(r'Velocity $\left(\frac{d\sigma}{d\tau_{B}}\right)$')
#plt.legend(loc=1)
#plt.xlim(0, max(timesteps))
#plt.ylim(0, pe)
#plt.show()

# Compute and overlay the average
avgInstV = float(instV) / float(instC)
avgHalfV = float(halfV) / float(halfC)
avgTaurV = float(taurV) / float(taurC)
avgTaubV = float(taubV) / float(taubC)

## Plot the data
#mysz = 15.
#x, y, z = zip(*actPos[start:end])
## Plot the time resolved active particle
#plt.scatter(x, y, c=(timesteps[start:end]/max(timesteps)), s=mysz, cmap='jet')
#plt.xlim(-hx_box, hx_box)
#plt.ylim(-hy_box, hy_box)
#ax = plt.gca()
#ax.set_aspect('equal')
#cbar = plt.colorbar(ticks=[], fraction=0.02675, pad=0.05, label='Time')
#plt.text(1.085, 0.015, r'$\tau_{initial}$', transform=ax.transAxes, fontsize=14)
#plt.text(1.085, 0.95, r'$\tau_{final}$', transform=ax.transAxes, fontsize=14)
#plt.show()

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
