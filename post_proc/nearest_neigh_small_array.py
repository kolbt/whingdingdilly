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

# initialize A/B_pos arrays
#pos_all = np.zeros((part_num, 2), dtype=np.float64)
#A_pos = np.zeros((part_a, 2), dtype=np.float64)
#B_pos = np.zeros((part_b, 2), dtype=np.float64)

size_min = 1000

# arrays to hold the avg neighbor data
avg_aa = np.zeros(dumps, dtype=np.float64)
avg_ab = np.zeros(dumps, dtype=np.float64)
avg_ba = np.zeros(dumps, dtype=np.float64)
avg_bb = np.zeros(dumps, dtype=np.float64)

position_array = np.zeros((1), dtype=np.ndarray)        # array of position arrays
type_array = np.zeros((1), dtype=np.ndarray)            # particle types
timesteps = np.zeros((dumps), dtype=np.float64)         # timesteps

with hoomd.open(name=myfile, mode='rb') as t:           # open for reading
    # analyze all particles
    for j in range(0, dumps):
#        start_tot = time.clock()
        snap = t[j]
        type_array = snap.particles.typeid
        position_array = snap.particles.position        # store all particle positions
        timesteps[j] = snap.configuration.step          # store tstep for plotting purposes

        part_num = len(type_array)
#        part_a = part_num * part_frac_a         # get the total number of A particles
#        part_a = int(part_a)
#        part_b = part_num - part_a              # get the total number of B particles
#        part_b = int(part_b)

        # Loading bars are fun :D
        time.sleep(1)
        percent = float(j) / float(dumps) * 100.0
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
        Aliq_count = 0
        Bliq_count = 0
        all_count = 0
        for c in range(0, part_num):
            if q_clust[ids[c]] == 1:
                all_count += 1
                if type_array[c] == 0:
                    Aliq_count += 1
                else:
                    Bliq_count += 1

    #    if all_count != 0:
    #        loop_count = 0
    #        all_pos = np.zeros((all_count, 2), dtype=np.float64)
    #        for c in range(0, part_num):
    #            if q_clust[ids[c]] == 1:
    #                all_pos[loop_count][0] = l_pos[c][0]
    #                all_pos[loop_count][1] = l_pos[c][1]
    #                loop_count += 1

        if Aliq_count != 0:
            loop_count = 0
            Aliq_pos = np.zeros((Aliq_count, 2), dtype=np.float64)
            for c in range(0, part_num):
                if q_clust[ids[c]] == 1:
                    if type_array[c] == 0:
                        Aliq_pos[loop_count, 0] = l_pos[c][0]
                        Aliq_pos[loop_count, 1] = l_pos[c][1]
                        loop_count += 1

        if Bliq_count != 0:
            loop_count = 0
            Bliq_pos = np.zeros((Bliq_count, 2), dtype=np.float64)
            for c in range(0, part_num):
                if q_clust[ids[c]] == 1:
                    if type_array[c] == 1:
                        Bliq_pos[loop_count, 0] = l_pos[c][0]
                        Bliq_pos[loop_count, 1] = l_pos[c][1]
                        loop_count += 1
    
        a_tree = spatial.KDTree(Aliq_pos)                      # tree of A-type particles
        b_tree = spatial.KDTree(Bliq_pos)                      # tree of B-type particles
        radius = 1.0
        
        # Let's look at the dense phase:
        # -how many A neighbors the avg A particle has and,
        '''
            97% of the total compuational time in this loop is spent here.
            The query between two trees is something I don't currently know
            how to outsource, plus, it only performs slowly on very large
            data sets (N~10^6).  So, for the time being, this clunky slow
            code will stay.
        '''
#        start = time.clock() ####
        ###########################################
        aa = a_tree.query_ball_tree(a_tree, radius)
        ###########################################
#        part_1 = time.clock() - start
        num_aa = np.array(map(len, aa), dtype=np.float64)
        num_aa -= 1.0                                       # can't reference itself
        if len(num_aa) != 0:
            avg_aa[j] = (np.sum(num_aa)/len(num_aa))
    
        # -how many B neighbors the avg A particle has and,
        ab = a_tree.query_ball_tree(b_tree, radius)
        num_ab = np.array(map(len, ab), dtype=np.float64)
        if len(num_ab) != 0:
            avg_ab[j] = (np.sum(num_ab)/len(num_ab))
    
        # -how many A neighbors the avg B particle has and,
        ba = b_tree.query_ball_tree(a_tree, radius)
        num_ba = np.array(map(len, ba), dtype=np.float64)
        if len(num_ba) != 0:
            avg_ba[j] = (np.sum(num_ba)/len(num_ba))

        # -how many B neighbors the avg B particle has.
        bb = b_tree.query_ball_tree(b_tree, radius)
        num_bb = np.array(map(len, bb), dtype=np.float64)
        num_bb -= 1.0                                       # can't reference itself
        if len(num_bb) != 0:
            avg_bb[j] = (np.sum(num_bb)/len(num_bb))

        #tot_time = time.clock() - start_tot
        #print("part_1 time is: " + str((part_1/tot_time)*100) + "%")

timesteps -= timesteps[0]
msd_time = timesteps[1:]

################################################################################
#################### Plot the individual and total data ########################
################################################################################

plt_name  = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a)
plt_name1 = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + "A"
plt_name2 = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + "B"

if part_perc_a != 0 and part_perc_a != 100:
    # plot some junk
    
    plt.plot(avg_aa, color="g")
    plt.plot(avg_ab, color="r")
    plt.plot(avg_ba, color="b")
    plt.plot(avg_bb, color="k")
    plt.savefig('avg_neighs_'+ plt_name + '.png', dpi=1000)
    plt.close()
    
#    # This is just for checking the cluster algorithm with a visual
#    fig, ax = plt.subplots()
#    fig.set_facecolor('black')
#    plt.subplots_adjust(top = 0.99, bottom = 0.01, right = 0.995, left = 0.005)
#    x = all_pos[:, 0]
#    y = all_pos[:, 1]
#    plt.scatter(x, y, s=1.5, c='b')
#    ax.set_xlim([left, right])
#    ax.set_ylim([left, right])
#    ax.get_xaxis().set_visible(False)
#    ax.get_yaxis().set_visible(False)
#    plt.savefig('all_' + plt_name1 + '.png', facecolor=fig.get_facecolor(), transparent=True, dpi=1000, box_inches = 'tight', edgecolor='none')
#    plt.close()

    fig, ax = plt.subplots()
    fig.set_facecolor('black')
    plt.subplots_adjust(top = 0.99, bottom = 0.01, right = 0.995, left = 0.005)
    x = Aliq_pos[:, 0]
    y = Aliq_pos[:, 1]
    z_a = num_aa
    plt.scatter(x, y, c=z_a, s=1.5, cmap='plasma')
    ax.set_xlim([left, right])
    ax.set_ylim([left, right])
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.colorbar()
    plt.savefig('aa_' + plt_name1 + '.png', facecolor=fig.get_facecolor(), transparent=True, dpi=1000, box_inches = 'tight', edgecolor='none')
    plt.close()

    fig, ax = plt.subplots()
    fig.set_facecolor('black')
    plt.subplots_adjust(top = 0.99, bottom = 0.01, right = 0.995, left = 0.005)
    x = Aliq_pos[:, 0]
    y = Aliq_pos[:, 1]
    z_a = num_ab
    plt.scatter(x, y, c=z_a, s=1.5, cmap='plasma')
    ax.set_xlim([left, right])
    ax.set_ylim([left, right])
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.colorbar()
    plt.savefig('ab_' + plt_name1 + '.png', facecolor=fig.get_facecolor(), transparent=True, dpi=1000, box_inches = 'tight', edgecolor='none')
    plt.close()

    fig, ax = plt.subplots()
    fig.set_facecolor('black')
    plt.subplots_adjust(top = 0.99, bottom = 0.01, right = 0.995, left = 0.005)
    x = Bliq_pos[:, 0]
    y = Bliq_pos[:, 1]
    z_a = num_ba
    plt.scatter(x, y, c=z_a, s=1.5, cmap='plasma')
    ax.set_xlim([left, right])
    ax.set_ylim([left, right])
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.colorbar()
    plt.savefig('ba_' + plt_name1 + '.png', facecolor=fig.get_facecolor(), transparent=True, dpi=1000, box_inches = 'tight', edgecolor='none')
    plt.close()

    fig, ax = plt.subplots()
    fig.set_facecolor('black')
    plt.subplots_adjust(top = 0.99, bottom = 0.01, right = 0.995, left = 0.005)
    x = Bliq_pos[:, 0]
    y = Bliq_pos[:, 1]
    z_a = num_bb
    plt.scatter(x, y, c=z_a, s=1.5, cmap='plasma')
#        ax.get_xlim()
#        ax.get_ylim()
    ax.set_xlim([left, right])
    ax.set_ylim([left, right])
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.colorbar()
    plt.savefig('bb_' + plt_name1 + '.png', facecolor=fig.get_facecolor(), transparent=True, dpi=1000, box_inches = 'tight', edgecolor='none')
    plt.close()

else:                                                           # if monodisperse plot total values
    # plot some other junk
    print("Shit")
