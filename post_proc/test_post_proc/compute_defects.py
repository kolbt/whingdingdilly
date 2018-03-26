'''
#                           This is an 80 character line                       #
This file will allow us to examine local divergence (both sources and sinks) 
more clearly:

    This script utilizes scipy's interpolate.griddata module to take our 
    discrete points and put them on a mesh.
    
    Then it employs numpy's gradient function to compute the divergence of the 
    vector field
'''

# Imports and loading the .gsd file
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
from scipy.interpolate import griddata

import matplotlib.pyplot as plt

import math

def quat_to_theta(quat):
    "Takes quaternion, returns angle"
    x = quat[1]
    y = quat[2]
    rad = math.atan2(y,x)/np.pi # gives values from [-1,1]
    return rad

def computeDivergence(field):
    "return the divergence of a n-D field"
    return np.sum(np.gradient(field),axis=0)

# File to read from
in_file = "pa"+str(pe_a)+\
"_pb"+str(pe_b)+\
"_xa"+str(part_perc_a)+\
".gsd"

f = hoomd.open(name=in_file, mode='rb') # open gsd file with hoomd
dumps = f.__len__()                     # get number of timesteps dumped

start = 0       # gives first frame to read
end = dumps     # gives last frame to read
start = 155
end = 160

positions = np.zeros((end), dtype=np.ndarray)       # array of positions
types = np.zeros((end), dtype=np.ndarray)           # particle types
box_data = np.zeros((1), dtype=np.ndarray)          # box dimensions
timesteps = np.zeros((end), dtype=np.float64)       # timesteps
orient = np.zeros((end), dtype=np.ndarray)          # orientations

# Get relevant data from .gsd file
with hoomd.open(name=in_file, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    for iii in range(start, end):
        snap = t[iii]                               # snapshot of frame
        types[iii] = snap.particles.typeid          # get types
        positions[iii] = snap.particles.position    # get positions
        orient[iii] = snap.particles.orientation    # get orientation
        timesteps[iii] = snap.configuration.step    # get timestep

timesteps -= timesteps[0]       # get rid of brownian run time

# Get number of each type of particle
part_num = len(types[start])
part_A = int(part_num * part_frac_a)
part_B = part_num - part_A

# Feed data into freud analysis software
l_box = box_data[0]
h_box = l_box / 2.0
a_box = l_box * l_box
f_box = box.Box(Lx = l_box, Ly = l_box, is2D = True)    # make freud box

# Make the meshgrid we'll feed in
interp = 100
xx, yy = np.mgrid[-h_box:h_box:interp*1j, -h_box:h_box:interp*1j]

for iii in range(start, end):
    
    # Convert orientation from quaternion to angle
    rads = np.zeros((part_num), dtype=np.float32)
    for jjj in range(0,part_num):
        rads[jjj] = quat_to_theta(orient[iii][jjj])
    
    # Get data from arrays
    pos = positions[iii]
    pos = np.delete(pos, 2, 1)
    typ = types[iii]
    tst = timesteps[iii]

    grid_z0 = griddata(pos, rads[:], (xx, yy), method='nearest')
    grid_z1 = griddata(pos, rads[:], (xx, yy), method='linear')
    grid_z2 = griddata(pos, rads[:], (xx, yy), method='cubic')
    
#    plt.subplot(221)
#    plt.imshow(pos.T, extent=(-h_box,h_box,-h_box,h_box), origin='lower')
#    plt.plot(pos[:,0], pos[:,1], 'k.', ms=0.1)
#    plt.title('Original')
#
#    plt.subplot(222)

    # Nearest interpolation looks like the best fit
    plt.imshow(grid_z0.T, extent=(-h_box,h_box,-h_box,h_box), origin='lower',
               cmap=plt.get_cmap('hsv'))
    plt.title('Nearest Interpolation of Discrete Data')
    
#    plt.subplot(223)
#    plt.imshow(grid_z1.T, extent=(-h_box,h_box,-h_box,h_box), origin='lower',
#               cmap=plt.get_cmap('hsv'))
#    plt.title('Linear')
#
#    plt.subplot(224)
#    plt.imshow(grid_z2.T, extent=(-h_box,h_box,-h_box,h_box), origin='lower',
#               cmap=plt.get_cmap('hsv'))
#    plt.title('Cubic')
#    plt.gcf().set_size_inches(6, 6)
#    plt.show()

    plt.savefig('interp' + str(interp) +
                '_pa'+ str(pe_a) +
                '_pb'+ str(pe_b) +
                '_xa'+ str(part_perc_a) +
                '_mvout_'+ str(iii) +
                '.png', dpi=1000)
    plt.close()

    '''Okay, cool. Now that you have the extrapolated orientation you need to
        compute the divergence on this mesh
        
        There are a few parameters you can tinker with to get this to work.
        1.) interp: how fine the mesh is
        2.) ck_rad: how many adjacent sites to compute divergence over
    '''

    diverge = (computeDivergence(grid_z0))
    plt.imshow(diverge.T, extent=(-h_box,h_box,-h_box,h_box), origin='lower',
               clim=(-3.0, 3.0))
    plt.title('Divergence of Interpolated Orientation field')
    plt.colorbar()
    plt.show()




