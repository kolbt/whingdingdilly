import sys

pe_a = int(sys.argv[1])                     # activity A
pe_b = int(sys.argv[2])                     # activity B
part_perc_a = int(sys.argv[3])              # percentage A particles
part_frac_a = float(part_perc_a) / 100.0    # fraction A particles
hoomd_path = str(sys.argv[4])               # local path to hoomd-blue
gsd_path = str(sys.argv[5])                 # local path to gsd

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
import math
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def quat_to_theta(quat):
    """Take a quaternion, output a value between [-1, 1]"""
    x = quat[1]
    y = quat[2]
    rad = math.atan2(y,x)/np.pi # gives values from [-1,1]
    return rad

def slowSort(array):
    """Sort an array the slow (but certain) way"""
    cpy = np.copy(array)
    ind = np.arange(0, len(array))
    for i in xrange(len(cpy)):
        for j in xrange(len(cpy)):
            if cpy[i] > cpy[j] and i < j:
                # Swap the copy array values
                tmp = cpy[i]
                cpy[i] = cpy[j]
                cpy[j] = tmp
                # Swap the corresponding indices
                tmp = ind[i]
                ind[i] = ind[j]
                ind[j] = tmp
    return ind

def indSort(arr1, arr2):
    """Take sorted index array, use to sort array"""
    # arr1 is array to sort
    # arr2 is index array
    cpy = np.copy(arr1)
    for i in xrange(len(arr1)):
        arr1[i] = cpy[arr2[i]]

def chkSort(array):
    """Make sure sort actually did its job"""
    for i in xrange(len(array)-2):
        if array[i] > array[i+1]:
            print("{} is not greater than {} for indices=({},{})").format(array[i+1], array[i], i, i+1)
            return False
    return True

try:
    # This is for the long timescale data
    long_file = "pa" + str(pe_a) + \
                "_pb" + str(pe_b) + \
                "_xa" + str(part_perc_a) + \
                ".gsd"
    # This is for the fine timescale data
    short_file = "log_pa" + str(pe_a) + \
                 "_pb" + str(pe_b) + \
                 "_xa" + str(part_perc_a) + \
                 ".gsd"
    # File to write all data to
    out_file = "orientation_pa" + str(pe_a) + \
               "_pb" + str(pe_b) + \
               "_xa" + str(part_perc_a) + \
               "_fm"
    f = hoomd.open(name=long_file, mode='rb') # open gsd file with hoomd
    g = hoomd.open(name=short_file, mode='rb') # open gsd file with hoomd
except:
    try:
        eps = str(sys.argv[6])
    except:
        eps = 1

    long_file = "pa" + str(pe_a) + \
                "_pb" + str(pe_b) + \
                "_xa" + str(part_perc_a) + \
                "_ep" + str(eps) + \
                ".gsd"
    short_file = "log_pa" + str(pe_a) + \
                 "_pb" + str(pe_b) + \
                 "_xa" + str(part_perc_a) + \
                 "_ep" + str(eps) + \
                 ".gsd"
    # File to write all data to
    out_file = "orientation_pa" + str(pe_a) + \
                "_pb" + str(pe_b) + \
                "_xa" + str(part_perc_a) + \
                "_ep" + str(eps) +\
                "_fm"
    f = hoomd.open(name=long_file, mode='rb')  # open gsd file with hoomd
    g = hoomd.open(name=short_file, mode='rb')  # open gsd file with hoomd

dump_long = int(f.__len__())                # get number of timesteps dumped
dump_short = int(g.__len__())               # get number of timesteps dumped

start = 0                       # gives first frame to read
end = dump_long + dump_short    # gives last frame to read

positions = np.zeros((end), dtype=np.ndarray)       # array of positions
types = np.zeros((end), dtype=np.ndarray)           # particle types
box_data = np.zeros((1), dtype=np.ndarray)          # box dimensions
timesteps = np.zeros((end), dtype=np.float64)       # timesteps
orientations = np.zeros((end), dtype=np.ndarray)      # orientations

# These make testing the script faster
stopOne = start + 2
stopTwo = dump_short + 2

# Get relevant data from short.gsd file
with hoomd.open(name=short_file, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    #    for iii in range(start, dump_short):
    for iii in range(start, stopOne):
        snap = t[iii]                                   # snapshot of frame
        types[iii] = snap.particles.typeid              # get types
        positions[iii] = snap.particles.position        # get positions
        timesteps[iii] = snap.configuration.step        # get timestep
        orientations[iii] = snap.particles.orientation  # get orientation
# Get relevant data from long.gsd file
with hoomd.open(name=long_file, mode='rb') as t:
    snap = t[0]
    for iii in range(dump_short, end):
        #    for iii in range(dump_short, stopTwo):
        snap = t[iii - dump_short]                      # snapshot of frame
        types[iii] = snap.particles.typeid              # get types
        positions[iii] = snap.particles.position        # get positions
        timesteps[iii] = snap.configuration.step        # get timestep
        orientations[iii] = snap.particles.orientation  # get orientation

partNum = len(types[start])
part_A = int(partNum * part_frac_a)
part_B = partNum - part_A

rads = np.zeros((end, partNum), dtype=np.float32)    # orientations as angles
for iii in range(dump_short, end):
    for jjj in range(0,partNum):
        rads[iii][jjj] = quat_to_theta(orientations[iii][jjj])

min = -(float(box_data[0]/2))
max = (float(box_data[0]/2))

for mmm in range(dump_short, end):
    xs = np.zeros((partNum), dtype = np.float32)
    ys = np.zeros((partNum), dtype = np.float32)
    for iii in range(0,partNum):
        xs[iii] = positions[mmm][iii][0]
        ys[iii] = positions[mmm][iii][1]

    myplt = plt.scatter(xs, ys, s=0.5, c=rads[mmm], cmap=plt.get_cmap('hsv'), edgecolors='none')
    plt.xlim(min, max)
    plt.ylim(min, max)
    ax = plt.gca()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_aspect('equal')
    cb = plt.colorbar(myplt, orientation='vertical')
    cb.set_ticks([])
    plt.savefig(out_file + str(mmm) +'.png', dpi=1000)
    plt.close()
