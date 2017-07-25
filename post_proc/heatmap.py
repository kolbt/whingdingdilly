import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#import sphviewer as sph
#def myplot(x, y, nb=32, xsize=500, ysize=500):
#    xmin = np.min(x)
#    xmax = np.max(x)
#    ymin = np.min(y)
#    ymax = np.max(y)
#    
#    x0 = (xmin+xmax)/2.
#    y0 = (ymin+ymax)/2.
#    
#    pos = np.zeros([3, len(x)])
#    pos[0,:] = x
#    pos[1,:] = y
#    w = np.ones(len(x))
#    
#    P = sph.Particles(pos, w, nb=nb)
#    S = sph.Scene(P)
#    S.update_camera(r='infinity', x=x0, y=y0, z=0,
#                    xsize=xsize, ysize=ysize)
#    R = sph.Render(S)
#    R.set_logscale()
#    img = R.get_image()
#    extent = R.get_extent()
#    for i, j in zip(xrange(4), [x0,x0,y0,y0]):
#        extent[i] += j
#    return img, extent

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

import scipy.spatial as spatial

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

# Figuring out how to mesh grid my data
n_divisions = 200
find_grid = float(l_box / 2)
split_x = split_y = np.linspace(-find_grid, find_grid, n_divisions)
#mesh = np.zeros((len(split_x), len(split_y), 2), dtype = np.float64)
flat_mesh = np.zeros((len(split_x)*len(split_y), 2), dtype = np.float64)
mesh_count = 0
for iii in range(0, len(split_x)):
    for kkk in range(0, len(split_y)):
        #mesh[iii, kkk, 0] = split_x[iii]
        #mesh[iii, kkk, 1] = split_y[kkk]
        flat_mesh[mesh_count, 0] = split_x[iii]
        flat_mesh[mesh_count, 1] = split_y[kkk]
        mesh_count += 1



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

    num_a_neigh = np.array(map(len, a_neighbors), dtype=np.float64)
    #num_a_neigh = float(num_a_neigh)
    num_a_neigh /= 7.0
    num_b_neigh = np.array(map(len, b_neighbors), dtype=np.float64)
    #num_b_neigh = float(num_b_neigh)
    num_b_neigh /= 7.0

    # let's use a mesh as a reference instead
    mesh_tree = spatial.KDTree(flat_mesh)
    a_neigh = mesh_tree.query_ball_tree(a_tree, radius)
    b_neigh = mesh_tree.query_ball_tree(b_tree, radius)
    mesh_a_num = np.array(map(len, a_neigh), dtype=np.float64)
    mesh_b_num = np.array(map(len, b_neigh), dtype=np.float64)

    ################################################################################
    #################### Plot the individual and total data ########################
    ################################################################################

    plt_name  = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a)
    plt_name1 = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + "A"
    plt_name2 = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + "B"

#    from astropy.convolution import convolve
#    from astropy.convolution.kernels import Gaussian2DKernel

    if part_perc_a != 0 and part_perc_a != 100:
        # plot some junk
        
        fig, ax = plt.subplots()
        fig.set_facecolor('black')
        plt.subplots_adjust(top = 0.99, bottom = 0.01, right = 0.995, left = 0.005)
        x = flat_mesh[:, 0]
        y = flat_mesh[:, 1]
        z_a = mesh_a_num
#        N = int(len(z_a)**.5)
#        z = z_a.reshape(N, N)
#        plt.imshow(z, cmap='plasma')
        plt.scatter(x, y, c=z_a, s=1.0, marker=',', cmap='plasma')
        ax.get_xlim()
        ax.get_ylim()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        plt.colorbar()
        plt.savefig('mesh_heatmap_' + plt_name1 + '.png', facecolor=fig.get_facecolor(), transparent=True, dpi=1000, box_inches = 'tight', edgecolor='none')
        plt.close()
        
        fig, ax = plt.subplots()
        fig.set_facecolor('black')
        plt.subplots_adjust(top = 0.995, bottom = 0.005, right = 0.995, left = 0.005)
        x = pos_all[:, 0]
        y = pos_all[:, 1]
        z_a = num_a_neigh
        plt.scatter(x, y, c=z_a, s=1.5, cmap='plasma')
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
        plt.scatter(x, y, c=z_b, s=1.5, cmap='plasma')
        ax.get_xlim()
        ax.get_ylim()
        #ax.set_xlim([-70,70])
        #ax.set_ylim([-70,70])
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        plt.colorbar()
        plt.savefig('heatmap_' + plt_name2 + '.png', facecolor=fig.get_facecolor(), transparent=True, dpi=1000, box_inches = 'tight', edgecolor='none')
        plt.close()
    
#        sns.jointplot(x=A_pos[:,0], y=A_pos[:,1], kind='hex')
#        plt.savefig('sbn_heatmap.png', transparent=True, dpi=1000, box_inches = 'tight', edgecolor='none')
#        plt.close()

#        lm = sns.kdeplot(A_pos[:,0], A_pos[:,1], cmap='plasma', n_levels=10, shade=True)
#        lm.set(xlim=(-find_grid, find_grid), ylim=(-find_grid, find_grid))
#        plt.savefig('sbn_heatmap.png', transparent=True, dpi=1000, box_inches = 'tight', edgecolor='none')
#        plt.close()

#        heatmap_16, extent_16 = myplot(A_pos[:,0], A_pos[:,1], nb=32)
#        plt.imshow(heatmap_16, extent=extent_16, origin='lower', cmap='plasma', aspect='auto')
#        plt.colorbar()
#        plt.savefig('heatmap_test.png', dpi = 1000)

#        heatmap, xedges, yedges = np.histogram2d(A_pos[:,0], A_pos[:,1], bins=100)
#        #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#
#        #plt.clf()
#        plt.imshow(convolve(heatmap, Gaussian2DKernel(stddev=2.0)), cmap='plasma', interpolation='none')
#        #plt.pcolormesh(xedges, yedges, heatmap.T, cmap='hot')
#        plt.colorbar()
#        plt.savefig('heatmap_test.png', dpi = 1000)
#        plt.close()

    else:                                                           # if monodisperse plot total values
        # plot some other junk
        print("Shit")
