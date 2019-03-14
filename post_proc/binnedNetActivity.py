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
import math
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors

def netActivity(numA, numB, PeA, PeB):
    "Takes binned number and type of particles, outputs net activity"
    if numA == 0 and numB == 0:
        return 0.0
    total = numA + numB
    fracA = float(numA) / float(total)
    fracB = 1.0 - fracA
    return (PeA * fracA) + (PeB * fracB)

def quatToVector(quat, type):
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
    "Take force vector, output magnitude"
    x = vecF[0]
    y = vecF[1]
    magF = np.sqrt((x**2)+(y**2))
    return magF

def quatToAngle(quat):
    "Take quaternion, output angle between [0, 2pi]"
    x = quat[1]
    y = quat[2]
    theta = math.atan2(y, x)
    theta += np.pi
    return theta

def vecToAngle(vec):
    "Take vector, output angle between [-pi, pi]"
    x = vec[0]
    y = vec[1]
    theta = math.atan2(y, x)
    return theta

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

# Make particle colormap
#colorsList = [(255,0,255),(50,205,50)]
colorsList = ['#FF00FF','#39FF14']
my_cmap = colors.ListedColormap(colorsList)

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
orient = np.zeros((end), dtype=np.ndarray)      # orientations

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
        orient[iii] = snap.particles.orientation        # get orientation
# Get relevant data from long.gsd file
with hoomd.open(name=long_file, mode='rb') as t:
    snap = t[0]
    for iii in range(dump_short, end):
        #    for iii in range(dump_short, stopTwo):
        snap = t[iii - dump_short]                      # snapshot of frame
        types[iii] = snap.particles.typeid              # get types
        positions[iii] = snap.particles.position        # get positions
        timesteps[iii] = snap.configuration.step        # get timestep
        orient[iii] = snap.particles.orientation        # get orientation

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
binWidth = 5.0
nBins = int(l_box / binWidth) + 1

for iii in range(dump_short, end):

    # Empty bin for net activity
    netPe = np.zeros((nBins, nBins), dtype=np.float32)
    binTyp = np.zeros((nBins, nBins, 2), dtype=np.float32)
    # Array to hold the magnitude of the binned force
    binnedF = np.zeros((nBins, nBins), dtype=np.float32)
    # And binned orientation
    binnedO = np.zeros((nBins, nBins), dtype=np.float32)
    # And binned force magnitude
    binnedM = np.zeros((nBins, nBins), dtype=np.float32)
    # And occupancy of each index
    binnedD = np.zeros((nBins, nBins), dtype=np.float32)
    
    # Easier accessors
    pos = positions[iii]
    pos = np.delete(pos, 2, 1)
    typ = types[iii]
    dir = orient[iii]
    act_vec = np.zeros((part_num, 2), dtype=np.float32)
    mesh = np.zeros((nBins, nBins, 2), dtype=np.float32)

    # Computes the x and y components of each particles' active force
    for jjj in range(0, part_num):
        act_vec[jjj] = quatToVector(dir[jjj], typ[jjj])

    # Take vector sum in each bin (will have x and y components)
    for jjj in range(0, part_num):
        # Get mesh indices
        tmp_posX = pos[jjj][0] + h_box
        tmp_posY = pos[jjj][1] + h_box
        x_ind = int(tmp_posX / binWidth)
        y_ind = int(tmp_posY / binWidth)
        # Sum vector active force
        mesh[x_ind][y_ind][0] += act_vec[jjj][0]
        mesh[x_ind][y_ind][1] += act_vec[jjj][1]
        # Sum magnitude of each force
        if typ[jjj] == 0:
            binnedM[x_ind][y_ind] += pe_a
            binTyp[x_ind][y_ind][0] += 1            # add to A count in bin
        else:
            binnedM[x_ind][y_ind] += pe_b
            binTyp[x_ind][y_ind][1] += 1            # add to B count in bin
        # Get occupancy of each index
        binnedD[x_ind][y_ind] += 1

    # Take magnitude of each bin (removes directional component)
    for jjj in range(0, nBins):
        for mmm in range(0, nBins):
            binnedF[jjj][mmm] += getMagnitude(mesh[jjj][mmm])
            binnedO[jjj][mmm] = vecToAngle(mesh[jjj][mmm])
            netPe[jjj][mmm] = netActivity(binTyp[jjj][mmm][0], binTyp[jjj][mmm][1], pe_a, pe_b)

    fig = plt.figure()
    # Plot the original simulation
    ax = fig.add_subplot(221)
    ax.scatter(pos[:,0], pos[:,1], c=typ, cmap=my_cmap, s=0.05, edgecolors='none')
    plt.xticks(())
    plt.yticks(())
    plt.xlim(-h_box, h_box)
    plt.ylim(-h_box, h_box)
    ax.set_aspect('equal')
    plt.title('Original')

    # Plot binned orientation
    ax = fig.add_subplot(222)
    plt.imshow(netPe.T,
               extent=(0,nBins,0,nBins),
               origin='lower'
               )
    plt.xticks(())
    plt.yticks(())
    cb = plt.colorbar()
#    cb.set_ticks([])
    ax.set_aspect('equal')
    plt.title('Net Activity')

    # Plot binned summed force magnitudes
    ax = fig.add_subplot(223)
    plt.imshow(binnedM.T,
               extent=(0,nBins,0,nBins),
               origin='lower')
    plt.xticks(())
    plt.yticks(())
    cb = plt.colorbar()
#    cb.set_ticks([])
    ax.set_aspect('equal')
    plt.title('Summed ||Force||')

    # Plot binned data using imshow
    ax = fig.add_subplot(224)
    plt.imshow(binnedF.T,
               extent=(0,nBins,0,nBins),
               origin='lower',
               clim=(0,10000))
    plt.xticks(())
    plt.yticks(())
    cb = plt.colorbar()
#    cb.set_ticks([])
    ax.set_aspect('equal')
    plt.title('Summed Vector Force')

    # Figure name
#    plt.savefig('nBins' + str(nBins) +
#                '_pa'+ str(pe_a) +
#                '_pb'+ str(pe_b) +
#                '_xa'+ str(part_perc_a) +
#                '_step_'+ str(iii) +
#                '.png', dpi=500)
    plt.savefig('forces'
                '_pa'+ str(pe_a) +
                '_pb'+ str(pe_b) +
                '_xa'+ str(part_perc_a) +
                '_eps' + str(eps) +
                '_fm'+ str(iii) +
                '.png', dpi=500)
    plt.close()

    # Plot binned orientation
    plt.imshow(netPe.T,
               extent=(0,nBins,0,nBins),
               origin='lower'
               )
    plt.xticks(())
    plt.yticks(())
    cb = plt.colorbar()
    #    cb.set_ticks([])
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.title('Net Activity')
    plt.savefig('netPe' +
                '_pa' + str(pe_a) +
                '_pb' + str(pe_b) +
                '_xa' + str(part_perc_a) +
                '_eps' + str(eps) +
                '_fm'+ str(iii) +
                '.png', dpi=500)
    plt.close()

#ffmpeg -framerate 10 -i nBins100_pa${pa}_pb${pb}_xa${xa}_step_%d.png\
# -vcodec libx264 -s 1000x1000 -pix_fmt yuv420p -threads 1\
# pressure_pa${pa}_pb${pb}_xa${xa}.mp4

#ffmpeg -framerate 10 -i nBins100_pa150_pb500_xa50_step_%d.png\
# -vcodec libx264 -s 1000x1000 -pix_fmt yuv420p -threads 1\
# nBins100_pa150_pb500_xa50.mp4

