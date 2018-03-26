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

def quat_to_vector(quat):
    "Takes quaternion, returns orientation vector"
    x = quat[1]
    y = quat[2]
    my_orient = (x, y)
    return my_orient

def theta_to_vec(angle):
    "Takes angle [-pi, pi], and returns (x, y) vector"
    x = np.cos(angle)
    y = np.sin(angle)
    vec = (x, y)
    return vec

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

# TEST SECTION: MAKE SAMPLE VECTOR FIELD AND COMPUTE DIVERGENCE #

# Make the vector field
#T = np.arctan2(yy, xx)
#R = 10 + np.sqrt((yy) ** 2 + (xx) ** 2)
#U, V = R * np.cos(T), R * np.sin(T)
#plt.quiver(xx, yy, U, V, R, alpha=.5)
#plt.quiver(xx, yy, U, V, edgecolor='k', facecolor='None', linewidth=.5)
#plt.xlim(-h_box, h_box)
#plt.xticks(())
#plt.ylim(-h_box, h_box)
#plt.yticks(())
#plt.show()
#
#diverge = (computeDivergence(grid_z0))
#plt.imshow(diverge.T, extent=(-h_box,h_box,-h_box,h_box), origin='lower',
#           clim=(-3.0, 3.0))
#plt.title('Divergence of Interpolated Orientation field')
#plt.colorbar()
#plt.show()

################### END TEST SECTION ############################

for iii in range(start, end):
    
    # Convert orientation from quaternion to angle
    rads = np.zeros((part_num), dtype=np.float32)
    for jjj in range(0,part_num):
        rads[jjj] = quat_to_theta(orient[iii][jjj])
    
    vecs = np.zeros((part_num), dtype=np.ndarray)
    for jjj in range(0,part_num):
        vecs[jjj] = quat_to_vector(orient[iii][jjj])
    
    # Get data from arrays
    pos = positions[iii]
    pos = np.delete(pos, 2, 1)
    typ = types[iii]
    tst = timesteps[iii]

    grid_z0 = griddata(pos, rads[:], (xx, yy), method='nearest')
    grid_z1 = griddata(pos, rads[:], (xx, yy), method='linear')
    grid_z2 = griddata(pos, rads[:], (xx, yy), method='cubic')

#    grid_z0 = griddata(pos, vecs[:], (xx, yy), method='nearest')
#    x_comp = []
#    y_comp = []
#    for i in range(0, interp):
#        for j in range(0, interp):
#            print(grid_z0[i][j][0])
#            x_comp.append(grid_z0[i][j][0])
#            y_comp.append(grid_z0[i][j][1])
#    plt.quiver(xx, yy, x_comp, y_comp)
#    plt.title('Nearest Interpolation of Discrete Data')
#    plt.show()
#    grid_z1 = griddata(pos, vecs[:], (xx, yy), method='linear')
#    grid_z2 = griddata(pos, vecs[:], (xx, yy), method='cubic')

#    plt.subplot(221)
#    plt.imshow(pos.T, extent=(-h_box,h_box,-h_box,h_box), origin='lower')
#    plt.plot(pos[:,0], pos[:,1], 'k.', ms=0.1)
#    plt.title('Original')
#
#    plt.subplot(222)

    # Nearest interpolation looks like the best fit
#    plt.imshow(grid_z0.T, extent=(-h_box,h_box,-h_box,h_box), origin='lower',
#               cmap=plt.get_cmap('hsv'))


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

#    plt.savefig('interp' + str(interp) +
#                '_pa'+ str(pe_a) +
#                '_pb'+ str(pe_b) +
#                '_xa'+ str(part_perc_a) +
#                '_mvout_'+ str(iii) +
#                '.png', dpi=1000)
#    plt.close()

    '''Okay, cool. Now that you have the extrapolated orientation you need to
        compute the divergence on this mesh
        
        There are a few parameters you can tinker with to get this to work.
        1.) interp: how fine the mesh is
        2.) ck_rad: how many adjacent sites to compute divergence over
    '''
    # Get gradient of angle scalar field
    grad_radx, grad_rady = np.gradient(grid_z0)
    
    # Take the magnitude
    grad_mag = grad_radx + grad_rady


    # Convert back to vectors
    grad_vecx = np.zeros((interp, interp), dtype=np.float32)
    grad_vecy = np.zeros((interp, interp), dtype=np.float32)
    defect = np.zeros((interp, interp), dtype=np.int8)
    for i in range(0, interp):
        for j in range(0, interp):
#            if grad_mag[i][j] > 1.0:
#                grad_mag[i][j] -= 2.0
#            if grad_mag[i][j] < -1.0:
#                grad_mag[i][j] += 2.0
            grad_vecx[i][j] = theta_to_vec(grad_mag[i][j]*np.pi)[0]
            grad_vecy[i][j] = theta_to_vec(grad_mag[i][j]*np.pi)[1]
            if -0.001 < grad_mag[i][j] < 0.001:
                defect[i][j] = 1
            if -0.9*2.0 > grad_mag[i][j] or grad_mag[i][j] > 0.9*2.0:
                defect[i][j] = -1

    diverge = np.pi * grad_mag
#    plt.quiver(xx, yy, grad_vecx, grad_vecy)
#    plt.show()

#    plt.imshow(grad_rad.T, extent=(-h_box,h_box,-h_box,h_box), origin='lower')
#    plt.show()

#    dx, dy = np.gradient(grid_z0)
#    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True)
#    for ax in axes.flat:
#        im = ax.imshow(dx)
#        im = ax.imshow(dy)
#    plt.subplot(121)    # row columns number
#    im = ax1.imshow(dx)
#    plt.subplot(122)    # row columns number
#    im = ax2.imshow(dy)
#    fig.colorbar(im, ax=axes.ravel().tolist())
#    plt.show()

#    diverge = (computeDivergence(grid_z0))
    plt.subplot(221)
    plt.quiver(xx, yy, grad_vecx, grad_vecy)
    plt.subplot(222)
    plt.imshow(diverge.T, extent=(-h_box,h_box,-h_box,h_box), origin='lower')
               #clim=(-3.0, 3.0))
#    plt.title('Divergence of Interpolated Orientation field')
    plt.colorbar()
    plt.subplot(223)
    plt.imshow(defect.T, extent=(-h_box,h_box,-h_box,h_box), origin='lower')
    plt.subplot(224)
    plt.plot(pos[:,0], pos[:,1], 'k.', ms=0.1)
    plt.show()




