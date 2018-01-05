'''
I want to compute the active pressure at every timestep.
Why?
    1.) Gives scalar per-particle active pressure
How?
    1.) Load particles: orientations and positions (for a tstep)
    2.) Pick an absolute reference point (0,0)
    3.) Loop through particles, compute:
        a. sigma_x = -R_x,Absolute * F_x,swim
        b. sigma_y = ...
        c. sigma = sigma_x + sigma_y
        d. pressure_active = 1/2 sigma
    4.) Now each particle has a scalar active pressure,
        plot it with a colormap at it's given position
Issue?
    Isn't this just random? Orientation is stochastic...
'''

# Imports and loading the .gsd file
import sys

hoomd_path = str(sys.argv[4])
gsd_path = str(sys.argv[5])

# need to extract values from filename (pa, pb, xa) for naming
part_perc_a = int(sys.argv[3])
part_frac_a = float(part_perc_a) / 100.0
pe_a = int(sys.argv[1])
pe_b = int(sys.argv[2])

sys.path.append(hoomd_path)

import hoomd
from hoomd import md
from hoomd import deprecated
from mypy import load_bar
sys.path.append(gsd_path)
import gsd
from gsd import hoomd
from gsd import pygsd
import numpy as np
import math

def quat_to_theta(quat):
    x = quat[1]
    y = quat[2]
    rad = math.atan2(y,x) # gives values from [-pi,pi]
    return rad

def computeKE(v):
    vx = v[0]
    vy = v[1]
    return vx**2 + vy**2

dumps = 1
boltz = 1
temp = 1
trans_diff = 1
drag = boltz * temp / trans_diff
my_dt = 0.000001

# Reading in the GSD file
myfile = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"
f = hoomd.open(name=myfile, mode='rb')
start = 625
dumps = f.__len__()
stop = 628
dumps = stop
position_array = np.zeros((dumps), dtype=np.ndarray)        # array of position arrays
type_array = np.zeros((dumps), dtype=np.ndarray)            # particle types
orientations = np.zeros((dumps), dtype=np.ndarray)          # orientations
velocities = np.zeros((dumps), dtype=np.ndarray)            # velocities
box_data = np.zeros((dumps), dtype=np.ndarray)              # box dimensions

with hoomd.open(name=myfile, mode='rb') as t:               # open for reading
    snap = t[0]                                             # snap oth snapshot
    box_data = snap.configuration.box                       # get box dimensions
    for i in range(start, dumps):
        snap = t[i]
        type_array[i] = snap.particles.typeid
        position_array[i] = snap.particles.position         # store all particle positions
        orientations[i] = snap.particles.orientation        # store the orientation of all particles
        velocities[i] = snap.particles.velocity             # store all particle velocities

part_num = len(position_array[start])

# PER PARTICLE QUANTITIES
active_pressure = np.zeros((dumps, part_num), dtype=np.float32) # holds active pressure
kinetic_nrg = np.zeros((dumps, part_num), dtype = np.float32)   # holds the binned KE
rads = np.zeros((dumps, part_num), dtype=np.float32)            # orientations as angles
tot_force = np.zeros((dumps-1, part_num, 3), dtype=np.float32)  # total force from velocity
pressure = np.zeros((dumps-1, part_num), dtype=np.float32)      # total pressure from force
diameter = 1

for iii in range(start,dumps):
    for jjj in range(0,part_num):
        rads[iii][jjj] = quat_to_theta(orientations[iii][jjj])
        if iii < (dumps - 1):
            tot_force[iii][jjj][0] = (velocities[iii+1][jjj][0] - velocities[iii][jjj][0]) / my_dt
            tot_force[iii][jjj][1] = (velocities[iii+1][jjj][1] - velocities[iii][jjj][1]) / my_dt
            #tot_force[iii][jjj][2] = (velocities[iii+1][jjj][2] - velocities[iii][jjj][2]) / my_dt
            pressure[iii][jjj] = (tot_force[iii][jjj][0] * position_array[iii][jjj][0]) +\
                (tot_force[iii][jjj][1] * position_array[iii][jjj][1])

for iii in range(start,dumps):
    for jjj in range(part_num):
        kinetic_nrg[iii][jjj] = computeKE(velocities[iii][jjj])
        if type_array[iii][jjj] == 0:
            sigma_x = drag * pe_a * math.cos(rads[iii][jjj]) * position_array[iii][jjj][0]
            sigma_y = drag * pe_a * math.sin(rads[iii][jjj]) * position_array[iii][jjj][1]
            active_pressure[iii][jjj] = -0.5 * (sigma_x + sigma_y)
        else:
            sigma_x = drag * pe_b * math.cos(rads[iii][jjj]) * position_array[iii][jjj][0]
            sigma_y = drag * pe_b * math.sin(rads[iii][jjj]) * position_array[iii][jjj][1]
            active_pressure[iii][jjj] = -0.5 * (sigma_x + sigma_y)

# Plot the data (per particle as a snapshot)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

min = -(float(box_data[0]/2))
max = (float(box_data[0]/2))

# Make a better colormap
import matplotlib.colors as col

for mmm in range(start,dumps):
    xs = np.zeros((part_num), dtype = np.float32)
    ys = np.zeros((part_num), dtype = np.float32)
    for iii in range(0,part_num):
        xs[iii] = position_array[mmm][iii][0]
        ys[iii] = position_array[mmm][iii][1]

#    plt.scatter(xs, ys, s=0.20, c=active_pressure[mmm], edgecolors='none')
#    plt.colorbar()
#    plt.xlim(min, max)
#    plt.ylim(min, max)
#    plt.savefig('active_pressure_pa'+
#                str(pe_a) +
#                '_pb'+
#                str(pe_b)+
#                '_xa'+
#                str(part_perc_a)
#                +'_mvout_'+
#                str(mmm)+
#                '.png',
#                dpi=1000)
#    plt.close()

#    plt.scatter(xs, ys, s=0.50, c=kinetic_nrg[mmm], cmap='viridis_r', edgecolors='none')
#    plt.colorbar()
#    plt.clim(0,10)
#    plt.xlim(min, max)
#    plt.ylim(min, max)
#    plt.savefig('kinetic_nrg_pa'+
#                str(pe_a) +
#                '_pb'+
#                str(pe_b)+
#                '_xa'+
#                str(part_perc_a)
#                +'_mvout_'+
#                str(mmm)+
#                '.png',
#                dpi=1000)
#    plt.close()

    if mmm < dumps - 1:
        plt.scatter(xs, ys, s=0.50, c=pressure[mmm], cmap='viridis_r', edgecolors='none')
        plt.colorbar()
#        plt.clim(0,10)
        plt.xlim(min, max)
        plt.ylim(min, max)
        plt.savefig('tot_press_pa'+
                    str(pe_a) +
                    '_pb'+
                    str(pe_b)+
                    '_xa'+
                    str(part_perc_a)
                    +'_mvout_'+
                    str(mmm)+
                    '.png',
                    dpi=1000)
        plt.close()







