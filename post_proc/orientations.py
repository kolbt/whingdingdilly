import sys

hoomd_path = str(sys.argv[4])
gsd_path = str(sys.argv[5])
file_path = str(sys.argv[6])

# need to extract values from filename (pa, pb, xa) for naming
part_perc_a = int(sys.argv[3])
part_frac_a = float(part_perc_a) / 100.0
pe_a = int(sys.argv[1])
pe_b = int(sys.argv[2])

sys.path.append(hoomd_path)

import hoomd
from hoomd import md
from hoomd import deprecated

#initialize system randomly, can specify GPU execution here

################################################################################
############################# Begin Data Analysis ##############################
################################################################################

sys.path.append(gsd_path)
import gsd
from gsd import hoomd
from gsd import pygsd
import numpy as np
import math

def quat_to_theta(quat):
    
    x = quat[1]
    y = quat[2]
    rad = math.atan2(y,x)/np.pi # gives values from [-1,1]

    return rad

myfile = str(file_path) + "/" + "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"
print(myfile)

f = hoomd.open(name=myfile, mode='rb')
dumps = f.__len__()
size_min = 10                                           # minimum size of cluster

position_array = np.zeros((dumps), dtype=np.ndarray)    # array of position arrays
type_array = np.zeros((dumps), dtype=np.ndarray)        # particle types
box_data = np.zeros((1), dtype=np.ndarray)              # box dimensions
timesteps = np.zeros((dumps), dtype=np.float64)         # timesteps
orientations = np.zeros((dumps), dtype=np.ndarray)      # orientations
velocities = np.zeros((dumps), dtype=np.ndarray)

start = 288
stop = 348

with hoomd.open(name=myfile, mode='rb') as t:           # open for reading
    snap = t[0]                                         # snap 0th snapshot
    box_data = snap.configuration.box                   # get box dimensions
    for i in range(start, stop):
        snap = t[i]                                     # take snap of each dump
        type_array[i] = snap.particles.typeid
        position_array[i] = snap.particles.position     # store all particle positions
        timesteps[i] = snap.configuration.step          # store tstep for plotting purposes
        orientations[i] = snap.particles.orientation    # store the orientation of all particles
        velocities[i] = snap.particles.velocity         #store the velocity of all particles

part_num = len(type_array[start])
rads = np.zeros((dumps, part_num), dtype=np.float32)    # orientations as angles

for iii in range(start, stop):
    for jjj in range(0,part_num):
        rads[iii][jjj] = quat_to_theta(orientations[iii][jjj])

part_a = part_num * part_frac_a         # get the total number of A particles
part_a = int(part_a)
part_b = part_num - part_a              # get the total number of B particles
part_b = int(part_b)

timesteps -= timesteps[0]

### I want to write a png with a circle at each position, the circle should be colored by it's orientation ###

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

min = -(float(box_data[0]/2))
max = (float(box_data[0]/2))

# Make a better colormap
import matplotlib.colors as col

white = '#ffffff'
black = '#000000'
red = '#ff0000'
blue = '#0000ff'
anglemap = col.LinearSegmentedColormap.from_list('anglemap', [black, red, blue, black], N=265, gamma=1)
plt.register_cmap(cmap=anglemap)

for mmm in range(start, stop):
    xs = np.zeros((part_num), dtype = np.float32)
    ys = np.zeros((part_num), dtype = np.float32)
    for iii in range(0,part_num):
        xs[iii] = position_array[mmm][iii][0]
        ys[iii] = position_array[mmm][iii][1]

    #plt.scatter(xs, ys, s=0.5, c=rads[mmm], cmap = anglemap)
    plt.scatter(xs, ys, s=0.25, c=rads[mmm], cmap=plt.get_cmap('hsv'), edgecolors='none')
#    plt.colorbar()
    plt.xlim(min, max)
    plt.ylim(min, max)
    plt.axes().get_xaxis().set_visible(False)
    plt.axes().get_yaxis().set_visible(False)
    plt.axes().set_aspect('equal')
    plt.savefig('orientation_pa'+ str(pe_a) +'_pb'+ str(pe_b) +'_xa'+ str(part_perc_a) +'_mvout_'+ str(mmm) +'.png', dpi=1000)
    plt.close()

#ffmpeg -framerate 8 -i orientation_%d.png -vcodec libx264 -pix_fmt yuv420p test_orientation.mov
#ffmpeg -framerate 8 -i orientation_pa80_pb300_xa50_mvout_%d.png -vcodec libx264 -pix_fmt yuv420p orientation_hsv.mov

