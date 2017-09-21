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
my_dt = 0.000001

################################################################################
############################# Begin Data Analysis ##############################
################################################################################

sys.path.append(gsd_path)
import gsd
from gsd import hoomd
from gsd import pygsd
import numpy as np

myfile = "MSD_pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"

f = hoomd.open(name=myfile, mode='rb')
dumps = f.__len__()
#size_min = 1000                                         # minimum size of cluster
size_min = 2

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
msd_time *= my_dt

part_num = len(type_array[0])

part_a = part_num * part_frac_a         # get the total number of A particles
part_a = int(part_a)
part_b = part_num - part_a              # get the total number of B particles
part_b = int(part_b)

pos_A = np.zeros((dumps), dtype=np.ndarray)             # type A positions
pos_B = np.zeros((dumps), dtype=np.ndarray)             # type B positions
tmpA = np.zeros((part_a, 3), dtype=np.float32)          # temporary storage arrays
tmpB = np.zeros((part_b, 3), dtype=np.float32)

from freud import parallel, box, density, cluster
parallel.setNumThreads(1)                               # don't run multiple threads

my_density = density.LocalDensity(r_cut=2.5,
                                  volume=0.79,
                                  diameter=1.0)         # initiate class, use area of circle

l_box = box_data[0]                                     # get box dimensions (square here)
f_box = box.Box(Lx=l_box,
                Ly=l_box,
                is2D=True)                               # initialize freud box

my_clusters = cluster.Cluster(box=f_box,
                              rcut=1.0)                 # initialize class
cluster_props = cluster.ClusterProperties(box=f_box)

ids = np.zeros((dumps), dtype=np.ndarray)
size_clusters = np.zeros((dumps), dtype=np.ndarray)
tot_size = np.zeros((dumps), dtype=np.ndarray)          # number of particles in clusters
tot_num = np.zeros((dumps), dtype=np.ndarray)           # total number of clusters
MCS = np.zeros((dumps), dtype=np.ndarray)               # Mean cluster size
GF = np.zeros((dumps), dtype=np.ndarray)                # Gas fraction
A_ids = np.zeros((part_a), dtype=np.ndarray)            # type A ids
B_ids = np.zeros((part_b), dtype=np.ndarray)            # type B ids
largest = np.zeros((dumps), dtype=np.ndarray)           # read out largest cluster at each tstep

# analyze all particles
for j in range(0, dumps):
    
    l_pos = position_array[j]
    my_clusters.computeClusters(l_pos)
    ids = my_clusters.getClusterIdx()                   # get cluster ids
    cluster_props.computeProperties(l_pos, ids)
    size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each

    #####################################################################
    ### Find avg cluster size, gas fraction, and largest cluster size ###
    #####################################################################
    l_clust = 0                                             # int size of largest cluster
    for k in range(0, len(size_clusters[j])):
        # the size minimum is a very important value to consider
        if size_clusters[j][k] > size_min and size_clusters[j][k] < part_num:
            tot_size[j] += size_clusters[j][k]
            tot_num[j] += 1
            if size_clusters[j][k] > l_clust:           # if larger cluster is found
                l_clust = size_clusters[j][k]           # set l_clust to that size
                    
    largest[j] = l_clust                                # save largest cluster size for tstep

    f_largest = "largest_pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".txt"
    if j == 0:
        a_w = 'w'
    else:
        a_w = 'a'
    f = open(f_largest, a_w)
    f.write(str(l_clust) + '\n')
    f.close()

    mcs_text = "MCS" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".txt"
    gf_text = "GF" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".txt"
    if tot_num[j] > 0:
        MCS[j] = float(tot_size[j]/tot_num[j])/float(part_num)
        GF[j] = float(part_num - tot_size[j]) / float(part_num)
    else:
        MCS[j] = 0.0001
        GF[j] = 1

    f = open(mcs_text, a_w)
    f.write(str(timesteps[j]) + ' ' + str(MCS[j]) + '\n')
    f.close()

    f = open(gf_text, a_w)
    f.write(str(timesteps[j]) + ' ' + str(GF[j]) + '\n')
    f.close()

################################################################################
###### perform the same analysis on species A and species B individually #######
################################################################################

if part_perc_a != 0 and part_perc_a != 100:

    tot_size_A = np.zeros((dumps), dtype=np.ndarray)          # number of particles in clusters
    tot_num_A = np.zeros((dumps), dtype=np.ndarray)           # total number of clusters
    MCS_A = np.zeros((dumps), dtype=np.ndarray)               # Mean cluster size
    GF_A = np.zeros((dumps), dtype=np.ndarray)                # Gas fraction

    tot_size_B = np.zeros((dumps), dtype=np.ndarray)          # number of particles in clusters
    tot_num_B = np.zeros((dumps), dtype=np.ndarray)           # total number of clusters
    MCS_B = np.zeros((dumps), dtype=np.ndarray)               # Mean cluster size
    GF_B = np.zeros((dumps), dtype=np.ndarray)                # Gas fraction

    for j in range(0, dumps):
        countA = 0
        countB = 0
        for g in range(0, part_num):
            if type_array[j][g] == 0:
                tmpA[countA][0] = position_array[j][g][0]
                tmpA[countA][1] = position_array[j][g][1]
                tmpA[countA][2] = position_array[j][g][2]
                countA += 1
            else:
                tmpB[countB][0] = position_array[j][g][0]
                tmpB[countB][1] = position_array[j][g][1]
                tmpB[countB][2] = position_array[j][g][2]
                countB += 1

        pos_A[j] = tmpA
        pos_B[j] = tmpB
        
        l_pos = pos_A[j]
        my_clusters.computeClusters(l_pos)
        ids = my_clusters.getClusterIdx()                   # get cluster ids
        cluster_props.computeProperties(l_pos, ids)
        size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each
        
        ####################################
        ### GF, MCS for A-A correlations ###
        ####################################
        
        for k in range(0, len(size_clusters[j])):
            # the size minimum is a very important value to consider
            if size_clusters[j][k] > size_min and size_clusters[j][k] < part_num:
                tot_size_A[j] += size_clusters[j][k]
                tot_num_A[j] += 1

        if tot_num_A[j] > 0:
            MCS_A[j] = float(tot_size_A[j]/tot_num_A[j])/float(part_a)
            GF_A[j] = float(part_a - tot_size_A[j]) / float(part_a)
        
        else:
            MCS_A[j] = 0.0001
            GF_A[j] = 1

        l_pos = pos_B[j]
        my_clusters.computeClusters(l_pos)
        ids = my_clusters.getClusterIdx()                   # get cluster ids
        cluster_props.computeProperties(l_pos, ids)
        size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each
        
        ####################################
        ### GF, MCS for A-A correlations ###
        ####################################

        for k in range(0, len(size_clusters[j])):
            # the size minimum is a very important value to consider
            if size_clusters[j][k] > size_min and size_clusters[j][k] < part_num:
                tot_size_B[j] += size_clusters[j][k]
                tot_num_B[j] += 1

        if tot_num_B[j] > 0:
            MCS_B[j] = float(tot_size_B[j]/tot_num_B[j])/float(part_b)
            GF_B[j] = float(part_b - tot_size_B[j]) / float(part_b)
        
        else:
            MCS_B[j] = 0.0001
            GF_B[j] = 1

################################################################################
#################### Plot the individual and total data ########################
################################################################################

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes=True)

plt_name  = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a)
plt_name1 = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + "A"
plt_name2 = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + "B"

if part_perc_a != 0 and part_perc_a != 100:

    plt.plot(msd_time, MCS[1:], color="g")
    plt.plot(msd_time, MCS_A[1:], color="r")
    plt.plot(msd_time, MCS_B[1:], color="b")
    #plt.ylim((0,1))
    plt.ylim(ymin=0.0001)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Time (tau)')
    plt.ylabel('MCS')
    plt.savefig('MCS_'+ plt_name + '.png', dpi=1000)
    plt.close()

    plt.plot(msd_time, GF[1:], color="g")
    plt.plot(msd_time, GF_A[1:], color="r")
    plt.plot(msd_time, GF_B[1:], color="b")
    plt.ylim((0,1))
    plt.xlabel('Time (tau)')
    plt.ylabel('GF')
    plt.savefig('GF_'+plt_name+'.png', dpi=1000)
    plt.close()

    plt.plot(msd_time, largest[1:], color="g")
    plt.xlabel('Time (tau)')
    plt.ylabel('Largest Cluster')
    plt.savefig('Largest_clust_'+plt_name+'.png', dpi=1000)
    plt.close()

else:                                                           # if monodisperse plot total values
    plt.plot(msd_time, MCS[1:], color="g")
    plt.ylim(ymin=0.0001)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Time (tau)')
    plt.ylabel('MCS')
    plt.savefig('MCS_'+ plt_name + '.png', dpi=1000)
    plt.close()

    plt.plot(msd_time, GF[1:], color="g")
    plt.ylim((0,1))
    plt.xlabel('Time (tau)')
    plt.ylabel('GF')
    plt.savefig('GF_'+plt_name+'.png', dpi=1000)
    plt.close()

    plt.plot(msd_time, largest[1:], color="g")
    plt.xlabel('Time (tau)')
    plt.ylabel('Largest Cluster')
    plt.savefig('Largest_clust_'+plt_name+'.png', dpi=1000)
    plt.close()
