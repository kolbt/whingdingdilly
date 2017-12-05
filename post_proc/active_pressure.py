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

dumps = 1
boltz = 1
temp = 1
trans_diff = 1
drag = boltz * temp / trans_diff

# Reading in the GSD file
myfile = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"
f = hoomd.open(name=myfile, mode='rb')
dumps = f.__len__()
position_array = np.zeros((dumps), dtype=np.ndarray)        # array of position arrays
type_array = np.zeros((dumps), dtype=np.ndarray)        # particle types
orientations = np.zeros((dumps), dtype=np.ndarray)          # orientations
box_data = np.zeros((dumps), dtype=np.ndarray)              # box dimensions

with hoomd.open(name=myfile, mode='rb') as t:           # open for reading
    snap = t[100]                                       # snap 100th snapshot
    box_data = snap.configuration.box                   # get box dimensions
    type_array[0] = snap.particles.typeid
    position_array[0] = snap.particles.position         # store all particle positions
    orientations[0] = snap.particles.orientation        # store the orientation of all particles

part_num = len(position_array[0])
rads = np.zeros((dumps, part_num), dtype=np.float32)    # orientations as angles
diameter = 1

#for iii in range(0,dumps):
for jjj in range(0,part_num):
    rads[0][jjj] = quat_to_theta(orientations[0][jjj])
#print(rads[0][0])
#print(math.cos(rads[0][0]))
#print(math.sin(rads[0][0]))

active_pressure = np.zeros((part_num), dtype=np.ndarray)

for iii in range(part_num):
    if type_array[0][iii] == 0:
        sigma_x = drag * pe_a * math.cos(rads[0][iii]) * position_array[0][iii][0]
        sigma_y = drag * pe_a * math.sin(rads[0][iii]) * position_array[0][iii][1]
        active_pressure[iii] = -0.5 * (sigma_x + sigma_y)
    else:
        sigma_x = drag * pe_b * math.cos(rads[0][iii]) * position_array[0][iii][0]
        sigma_y = drag * pe_b * math.sin(rads[0][iii]) * position_array[0][iii][1]
        active_pressure[iii] = -0.5 * (sigma_x + sigma_y)

# Plot the data (per particle as a snapshot)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

min = -(float(box_data[0]/2))
max = (float(box_data[0]/2))

# Make a better colormap
import matplotlib.colors as col

#white = '#ffffff'
#black = '#000000'
#red = '#ff0000'
#blue = '#0000ff'
#anglemap = col.LinearSegmentedColormap.from_list('anglemap',
#                                                 [black, red, blue, black],
#                                                 N=265,
#                                                 gamma=1)
#plt.register_cmap(cmap=anglemap)

for mmm in range(0,dumps):
    xs = np.zeros((part_num), dtype = np.float32)
    ys = np.zeros((part_num), dtype = np.float32)
    for iii in range(0,part_num):
        xs[iii] = position_array[mmm][iii][0]
        ys[iii] = position_array[mmm][iii][1]
    
    #plt.scatter(xs, ys, s=0.5, c=rads[mmm], cmap = anglemap)
    plt.scatter(xs, ys, s=0.20, c=active_pressure, edgecolors='none')
#                cmap=plt.get_cmap('hsv'))
    plt.colorbar()
    plt.xlim(min, max)
    plt.ylim(min, max)
    plt.savefig('active_pressure_pa'+
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







