'''
#                           This is an 80 character line                       #
This is going to be SUPER easy
    -look at particle type
    -get orientation
    -vectorize active force
    -sum force in each bin (this accounts for orientation)
    -take magnitude of summed quantity
    
    We'll plot this along with some other contributing heatmaps
    1. Original sim
    1. Binned orientation
    3. Binned active force magnitudes
    4. Magnitude of binned vector forces
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

def quat_to_vector(quat, type):
    "Takes quaternion, returns orientation vector"
    if type == 0:
        mag = pe_a
    else:
        mag = pe_b
    x = quat[1] * mag
    y = quat[2] * mag
    act_vec = (x, y)
    return act_vec

def getMagnitude(vecF):
    x = vecF[0]
    y = vecF[1]
    magF = np.sqrt((x**2)+(y**2))
    return magF

# Make particle colormap
cdict1 = {'0':  '#FF00FF',
          '1':  '#39FF14'
          }
plt.register_cmap(name='Particles', data=cdict1)


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

# Make the mesh
nBins = 100
sizeBin = l_box / nBins



for iii in range(start, end):

    # Array to hold the magnitude of the binned force
    binnedF = np.zeros((nBins, nBins), dtype=np.float32)
    
    # Easier accessors
    pos = positions[iii]
    pos = np.delete(pos, 2, 1)
    typ = types[iii]
    dir = orient[iii]
    act_vec = np.zeros((part_num, 2), dtype=np.float32)
    mesh = np.zeros((nBins, nBins, 2), dtype=np.float32)

    # Computes the x and y components of each particles' active force
    for jjj in range(0, part_num):
        act_vec[jjj] = quat_to_vector(dir[jjj], typ[jjj])

    # Take vector sum in each bin (will have x and y components)
    for jjj in range(0, part_num):
        tmp_posX = pos[jjj][0] + h_box
        tmp_posY = pos[jjj][1] + h_box
        x_ind = int(tmp_posX / sizeBin)
        y_ind = int(tmp_posY / sizeBin)
        mesh[x_ind][y_ind][0] += act_vec[jjj][0]
        mesh[x_ind][y_ind][1] += act_vec[jjj][1]

    # Take magnitude of each bin (removes directional component)
    for jjj in range(0, nBins):
        for mmm in range(0, nBins):
            binnedF[jjj][mmm] += getMagnitude(mesh[jjj][mmm])

    fig = plt.figure()
    cmap = plt.get_cmap('Particles')
    # Plot the original simulation
    ax = fig.add_subplot(221)
    ax.scatter(pos[:,0], pos[:,1], c=typ, cmap=cmap, s=0.05, edgecolors='none')
    plt.xticks(())
    plt.yticks(())
    plt.xlim(-h_box, h_box)
    plt.ylim(-h_box, h_box)
    ax.set_aspect('equal')

    # Plot binned orientation
    ax = fig.add_subplot(222)
    ax.scatter(pos[:,0], pos[:,1], c=typ, s=0.05, edgecolors='none')
    plt.xticks(())
    plt.yticks(())
    plt.xlim(-h_box, h_box)
    plt.ylim(-h_box, h_box)
    ax.set_aspect('equal')

    # Plot binned summed force magnitudes
    ax = fig.add_subplot(223)
    ax.scatter(pos[:,0], pos[:,1], c=typ, s=0.05, edgecolors='none')
    plt.xticks(())
    plt.yticks(())
    plt.xlim(-h_box, h_box)
    plt.ylim(-h_box, h_box)
    ax.set_aspect('equal')

    # Plot binned data using imshow
    ax = fig.add_subplot(224)
    plt.imshow(binnedF.T,
               extent=(0,nBins,0,nBins),
               origin='lower',
               clim=(0,10000))
    plt.xticks(())
    plt.yticks(())
    cb = plt.colorbar()
    cb.set_ticks([])
    ax.set_aspect('equal')

    # Figure name
#    plt.show()
    plt.savefig('nBins' + str(nBins) +
                '_pa'+ str(pe_a) +
                '_pb'+ str(pe_b) +
                '_xa'+ str(part_perc_a) +
                '_step_'+ str(iii) +
                '.png', dpi=1000)
    plt.close()

#ffmpeg -framerate 10 -i nBins100_pa${pa}_pb${pb}_xa${xa}_step_%d.png\
# -vcodec libx264 -s 1000x1000 -pix_fmt yuv420p -threads 1\
# pressure_pa${pa}_pb${pb}_xa${xa}.mp4

#ffmpeg -framerate 10 -i nBins100_pa150_pb500_xa50_step_%d.png\
# -vcodec libx264 -s 1000x1000 -pix_fmt yuv420p -threads 1\
# nBins100_pa150_pb500_xa50.mp4
