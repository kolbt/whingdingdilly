import sys
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

from freud import parallel, box, density, cluster
parallel.setNumThreads(1)                               # don't run multiple threads

box_data = np.zeros((1), dtype=np.ndarray)              # box dimensions
with hoomd.open(name=myfile, mode='rb') as t:           # open for reading
    snap = t[0]                                         # snap 0th snapshot
    box_data = snap.configuration.box                   # get box dimensions

l_box = box_data[0]                                     # get box dimensions (square here)
left = -(l_box/2)
right = (l_box/2)

f_box = box.Box(Lx=l_box,
                Ly=l_box,
                is2D=True)                              # initialize freud box

my_clusters = cluster.Cluster(box=f_box,
                              rcut=1.0)                 # initialize class
cluster_props = cluster.ClusterProperties(box=f_box)
ids = np.zeros((dumps), dtype=np.ndarray)
size_clusters = np.zeros((dumps), dtype=np.ndarray)

size_min = int(sys.argv[6])

position_array = np.zeros((1), dtype=np.ndarray)        # array of position arrays
type_array = np.zeros((1), dtype=np.ndarray)            # particle types
timesteps = np.zeros((dumps), dtype=np.float64)         # timesteps

with hoomd.open(name=myfile, mode='rb') as t:           # open for reading
    # analyze all particles
    for j in range(0, dumps):
        snap = t[j]
        type_array = snap.particles.typeid
        position_array = snap.particles.position        # store all particle positions
        timesteps[j] = snap.configuration.step          # store tstep for plotting purposes

        part_num = len(type_array)

        # Loading bars are fun :D
        time.sleep(1)
        percent = float(j)/float(dumps-1) * 100.0
        if percent == 100.0:
            sys.stdout.write("\r%5.1f%%\n" % percent)
            sys.stdout.flush()
        else:
            sys.stdout.write("\r%5.1f%%" % percent)
            sys.stdout.flush()
    
        l_pos = position_array
        my_clusters.computeClusters(l_pos)
        
        ids = my_clusters.getClusterIdx()                   # get cluster ids
        cluster_props.computeProperties(l_pos, ids)         # compute cluster properties
        size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each
        how_many = my_clusters.getNumClusters()             # how many clusters are there?
        
        sort_id = np.sort(ids)                              # array of IDs sorted small to large
        q_clust = np.zeros((how_many), dtype=np.ndarray)    # my binary 'is it clustered?' array
        
        index = 0                                           # index of the sorted array to look at
        for a in range(0,len(q_clust)):
            add_clust = 0
            while 1:
                add_clust += 1
                if index == part_num:                       # break if index is too large
                    break
                if sort_id[index] != a:                     # break if ID changes
                    break
                if add_clust == 1:                          # all particles appear once
                    q_clust[a] = 0
                if add_clust > size_min:                    # only multiple ids appear twice
                    q_clust[a] = 1
                index += 1                                  # increment index

        # This will get me the length of each array
        all_count = 0
        for c in range(0, part_num):
            if q_clust[ids[c]] == 1:
                all_count += 1

        str_j = str(j)
        pad_j = str_j.zfill(3)

        if all_count != 0:
            loop_count = 0
            all_pos = np.zeros((all_count, 2), dtype=np.float64)
            for c in range(0, part_num):
                if q_clust[ids[c]] == 1:
                    all_pos[loop_count][0] = l_pos[c][0]
                    all_pos[loop_count][1] = l_pos[c][1]
                    loop_count += 1
        
            # This is just for checking the cluster algorithm with a visual
            fig, ax = plt.subplots()
            fig.set_facecolor('black')
            plt.subplots_adjust(top = 0.99, bottom = 0.01, right = 0.995, left = 0.005)
            x_tot = l_pos[:, 0]
            y_tot = l_pos[:, 1]
            plt.scatter(x_tot, y_tot, s=1.5, c='w')
            x = all_pos[:, 0]
            y = all_pos[:, 1]
            plt.scatter(x, y, s=1.5, c='c')
            ax.set_xlim([left, right])
            ax.set_ylim([left, right])
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_aspect(aspect='equal')
            plt.savefig('sm_' + str(size_min) + '_opt_'+ pad_j +'.png',
                        facecolor=fig.get_facecolor(),
                        transparent=True,
                        dpi=72,
                        box_inches = 'tight',
                        edgecolor='none')
            plt.close()
    
        else:
            fig, ax = plt.subplots()
            fig.set_facecolor('black')
            plt.subplots_adjust(top = 0.99, bottom = 0.01, right = 0.995, left = 0.005)
            x_tot = l_pos[:, 0]
            y_tot = l_pos[:, 1]
            plt.scatter(x_tot, y_tot, s=1.5, c='w')
            ax.set_xlim([left, right])
            ax.set_ylim([left, right])
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_aspect(aspect='equal')
            plt.savefig('sm_' + str(size_min) + '_opt_'+ pad_j +'.png',
                        facecolor=fig.get_facecolor(),
                        transparent=True,
                        dpi=72,
                        box_inches = 'tight',
                        edgecolor='none')
            plt.close()


#timesteps -= timesteps[0]
#msd_time = timesteps[1:]
