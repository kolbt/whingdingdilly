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
#from mypy import load_bar
sys.path.append(gsd_path)
import gsd
from gsd import hoomd
from gsd import pygsd
import numpy as np
import math
import gtar

def quat_to_theta(quat):
    x = quat[1]
    y = quat[2]
    rad = math.atan2(y,x) # gives values from [-pi,pi]
    return rad

def computeKE(v):
    vx = v[0]
    vy = v[1]
    return vx**2, vy**2

def computePressure(ke_x, ke_y, xx, yy):
    x_dir = (ke_x + xx)
    y_dir = (ke_y + yy)
    return x_dir + y_dir

dumps = 1
boltz = 1
temp = 1
trans_diff = 1
drag = boltz * temp / trans_diff
my_dt = 0.000001

# Reading in the GSD file
mygsd = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"
mysqlite = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".sqlite"
f = hoomd.open(name=mygsd, mode='rb')
start = 0
dumps = f.__len__()
stop = dumps
dumps = stop
position_array = np.zeros((dumps), dtype=np.ndarray)        # array of position arrays
type_array = np.zeros((dumps), dtype=np.ndarray)            # particle types
orientations = np.zeros((dumps), dtype=np.ndarray)          # orientations
velocities = np.zeros((dumps), dtype=np.ndarray)            # velocities
box_data = np.zeros((dumps), dtype=np.ndarray)              # box dimensions

with hoomd.open(name=mygsd, mode='rb') as t:                # open for reading
    snap = t[0]                                             # snap oth snapshot
    box_data = snap.configuration.box                       # get box dimensions
    for i in range(start, dumps):
        snap = t[i]
        type_array[i] = snap.particles.typeid
        position_array[i] = snap.particles.position         # store all particle positions
        orientations[i] = snap.particles.orientation        # store the orientation of all particles
        velocities[i] = snap.particles.velocity             # store all particle velocities

part_num = len(position_array[start])
pressure = np.zeros((dumps, part_num), dtype=np.float32)

count = 0
#load_bar.printLoadBar(0, dumps, prefix = "Progress:", suffix = "Complete")
with gtar.GTAR(mysqlite, 'r') as traj:
    for (frame, (virial, velocity)) in traj.recordsNamed(['virial', 'velocity']):
        for part in range(len(velocity)):
            xx = float(virial[part][0])
            yy = float(virial[part][3])
            v = velocity[part]
            ke_x, ke_y = computeKE(v)
            pressure[count][part] = computePressure(ke_x, ke_y, xx, yy)
#        load_bar.printLoadBar(count+1, dumps, prefix = "Progress:", suffix = "Complete")
        count += 1

pressure /= ((float(box_data[0]))**2)

# Plot the data (per particle as a snapshot)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

min = -(float(box_data[0]/2))
max = (float(box_data[0]/2))

# Make a better colormap
import matplotlib.colors as col

#load_bar.printLoadBar(0, dumps, prefix = "Progress:", suffix = "Complete")
for mmm in range(start,dumps):
    xs = np.zeros((part_num), dtype = np.float32)
    ys = np.zeros((part_num), dtype = np.float32)
    
    for iii in range(0,part_num):
        xs[iii] = position_array[mmm][iii][0]
        ys[iii] = position_array[mmm][iii][1]

    plt.scatter(xs, ys, s=0.5, c=pressure[mmm], cmap='viridis_r', edgecolors='none')
#plt.colorbar()
    plt.clim(0,0.025)
    plt.xlim(min, max)
    plt.ylim(min, max)
    plt.axes().get_xaxis().set_visible(False)
    plt.axes().get_yaxis().set_visible(False)
    plt.axes().set_aspect('equal')
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
    #load_bar.printLoadBar(mmm+1, dumps, prefix = "Progress:", suffix = "Complete")
