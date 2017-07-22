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

#initialize system randomly, can specify GPU execution here

part_num = 15000

part_a = part_num * part_frac_a         # get the total number of A particles
part_a = int(part_a)
part_b = part_num - part_a              # get the total number of B particles
part_b = int(part_b)

################################################################################
############################# Begin Data Analysis ##############################
################################################################################

sys.path.append(gsd_path)
import gsd
from gsd import hoomd
from gsd import pygsd
import numpy as np
import scipy.spatial as spatial
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes=True)


myfile = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"

f = hoomd.open(name=myfile, mode='rb')
dumps = f.__len__()

position_array = np.zeros((dumps), dtype=np.ndarray)    # array of position arrays
type_array = np.zeros((dumps), dtype=np.ndarray)        # particle types
box_data = np.zeros((1), dtype=np.ndarray)              # box dimensions
timesteps = np.zeros((dumps), dtype=np.float64)         # timesteps

with hoomd.open(name=myfile, mode='rb') as t:           # open for reading
    snap = t[0]                                         # snap 0th snapshot
    box_data = snap.configuration.box                   # get box dimensions
    for i in range(0,dumps):
        snap = t[i]                                     # take snap of each dump
        type_array[i] = snap.particles.typeid
        position_array[i] = snap.particles.position     # store all particle positions
        timesteps[i] = snap.configuration.step          # store tstep for plotting purposes

timesteps -= timesteps[0]
msd_time = timesteps[1:]

from freud import parallel, box, density, cluster
parallel.setNumThreads(1)                               # don't run multiple threads

l_box = box_data[0]                                     # get box dimensions (square here)
f_box = box.Box(Lx=l_box,
                Ly=l_box,
                is2D=True)                              # initialize freud box

# initialize A/B_pos arrays
pos_all = np.zeros((part_num, 2), dtype=np.float64)
A_pos = np.zeros((part_a, 2), dtype=np.float64)
B_pos = np.zeros((part_b, 2), dtype=np.float64)

# analyze all particles
for j in range(dumps-1, dumps):
    
    l_pos = position_array[j]
    a_count = 0
    b_count = 0
    
    for i in range(0, part_num):
        pos_all[i][0] = l_pos[i][0]
        pos_all[i][1] = l_pos[i][1]
        if type_array[j][i] == 0:
            A_pos[a_count][0]=l_pos[i][0]
            A_pos[a_count][1]=l_pos[i][1]
            a_count += 1
        else:
            B_pos[b_count][0]=l_pos[i][0]
            B_pos[b_count][1]=l_pos[i][1]
            b_count += 1

    tree = spatial.KDTree(pos_all)                      # tree of all points
    a_tree = spatial.KDTree(A_pos)                      # tree of A-type particles
    b_tree = spatial.KDTree(B_pos)                      # tree of B-type particles
    radius = 1.0
    
    a_neighbors = tree.query_ball_tree(a_tree, radius)
    b_neighbors = tree.query_ball_tree(b_tree, radius)

    num_a_neigh = np.array(map(len, a_neighbors))
    num_b_neigh = np.array(map(len, b_neighbors))

    ################################################################################
    #################### Plot the individual and total data ########################
    ################################################################################

    plt_name  = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a)
    plt_name1 = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + "A"
    plt_name2 = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + "B"

    if part_perc_a != 0 and part_perc_a != 100:
        # plot some junk
        fig, ax = plt.subplots()
        fig.set_facecolor('black')
        plt.subplots_adjust(top = 0.995, bottom = 0.005, right = 0.995, left = 0.005)
        x = pos_all[:, 0]
        y = pos_all[:, 1]
        z_a = num_a_neigh
        plt.scatter(x, y, c=z_a, s=1.5, cmap='hot')
        ax.get_xlim()
        ax.get_ylim()
        #ax.set_xlim([-70,70])
        #ax.set_ylim([-70,70])
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        plt.colorbar()
        plt.savefig('heatmap_' + plt_name1 + '.png', facecolor=fig.get_facecolor(), transparent=True, dpi=1000, box_inches = 'tight', edgecolor='none')
        plt.close()

        fig, ax = plt.subplots()
        fig.set_facecolor('black')
        plt.subplots_adjust(top = 0.995, bottom = 0.005, right = 0.995, left = 0.005)
        x = pos_all[:, 0]
        y = pos_all[:, 1]
        z_b = num_b_neigh
        plt.scatter(x, y, c=z_b, s=1.5, cmap='hot')
        ax.get_xlim()
        ax.get_ylim()
        #ax.set_xlim([-70,70])
        #ax.set_ylim([-70,70])
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        plt.colorbar()
        plt.savefig('heatmap_' + plt_name2 + '.png', facecolor=fig.get_facecolor(), transparent=True, dpi=1000, box_inches = 'tight', edgecolor='none')
        plt.close()

    else:                                                           # if monodisperse plot total values
        # plot some other junk
        print("Shit")
