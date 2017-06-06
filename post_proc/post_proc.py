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

#########################################################################
########################## Begin Data Analysis ##########################
#########################################################################

sys.path.append(gsd_path)
import gsd
from gsd import hoomd
from gsd import pygsd
import numpy as np

myfile = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"

f = hoomd.open(name=myfile, mode='rb')
dumps = f.__len__()

position_array = np.zeros((dumps), dtype=np.ndarray)    # array of position arrays
type_array = np.zeros((dumps), dtype=np.ndarray)        # particle types
box_data = np.zeros((1), dtype=np.ndarray)              # box dimensions

with hoomd.open(name=myfile, mode='rb') as t:           # open for reading
    snap = t[0]                                         # snap 0th snapshot
    box_data = snap.configuration.box                   # get box dimensions
    for i in range(0,dumps):
        snap = t[i]                                     # take snap of each dump
        type_array[i] = snap.particles.typeid
        position_array[i] = snap.particles.position     # store all particle positions

pos_A = np.zeros((dumps), dtype=np.ndarray)             # type A positions
pos_B = np.zeros((dumps), dtype=np.ndarray)             # type B positions
tmpA = np.zeros((part_a, 3), dtype=np.float32)          # temporary storage arrays
tmpB = np.zeros((part_b, 3), dtype=np.float32)

from freud import parallel, box, density, cluster, kspace
parallel.setNumThreads(1)                               # don't run multiple threads

my_density = density.LocalDensity(r_cut=2.5,
                                  volume=0.79,
                                  diameter=1.0)         # initiate class, use area of circle

l_box = box_data[0]                                     # get box dimensions (square here)
f_box = box.Box(Lx=l_box,
                Ly=l_box,
                is2D=True)                               # initialize freud box

sf_box = box.Box(Lx=l_box,
                 Ly=l_box,
                 Lz=0.000005)

my_clusters = cluster.Cluster(box=f_box,
                              rcut=1.0)                 # initialize class
cluster_props = cluster.ClusterProperties(box=f_box)
sfactor = kspace.SFactor3DPoints(box=sf_box,
                                 g=1)

number_clusters = np.zeros((dumps), dtype=np.ndarray)   # arrays to store things
ids = np.zeros((dumps), dtype=np.ndarray)
size_clusters = np.zeros((dumps), dtype=np.ndarray)
tot_size = np.zeros((dumps), dtype=np.ndarray)          # number of particles in clusters
tot_num = np.zeros((dumps), dtype=np.ndarray)           # total number of clusters
MCS = np.zeros((dumps), dtype=np.ndarray)               # Mean cluster size
GF = np.zeros((dumps), dtype=np.ndarray)                # Gas fraction
A_ids = np.zeros((part_a), dtype=np.ndarray)            # type A ids
B_ids = np.zeros((part_b), dtype=np.ndarray)            # type B ids
percent_A = np.zeros((dumps), dtype=np.ndarray)         # composition A at each timestep
largest = np.zeros((dumps), dtype=np.ndarray)           # read out largest cluster at each tstep
MSD = np.zeros((dumps - 1, part_num), dtype=np.ndarray) # array of individual particle MSDs
MSD_A = np.zeros((dumps - 1, part_a), dtype=np.ndarray) # array for a particles
MSD_B = np.zeros((dumps - 1, part_b), dtype=np.ndarray) # array for a particles
tmp_struct_fact = np.zeros((100), dtype=np.ndarray)     # static sf at each timestep
struct_fact = np.zeros((100), dtype=np.ndarray)         # get magnitude keep copy
sf_count = 0

# I've altered this to take the last 100 dumps, change this as you see fit
# Note: there's something screwy going on with the sf, took the det of the matrix

# analyze all particles
for j in range(0, dumps):
    
    l_pos = position_array[j]
    my_clusters.computeClusters(l_pos)
    number_clusters[j] = my_clusters.getNumClusters()   # find number of clusters
    ids = my_clusters.getClusterIdx()                   # get cluster ids
    cluster_props.computeProperties(l_pos, ids)
    size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each
    
    if j >= dumps - 100:
        sfactor.compute(l_pos)
        tmp_struct_fact[sf_count] = sfactor.getS()
        struct_fact[sf_count] = np.linalg.det(tmp_struct_fact[sf_count][0])
        sf_count += 1

#    how_many = my_clusters.getNumClusters()
#    
#    A_id_count = 0
#    B_id_count = 0
#    for h in range(0, part_num):
#        if type_array[j][h] == 0:
#            A_ids[A_id_count] = ids[h]                  # store the cluster ids for A type
#            A_id_count += 1                             # IMPROVE: sort while placing?
#        else:
#            B_ids[B_id_count] = ids[h]                  # store the cluster ids for B type
#            B_id_count += 1                             # could put ids in order ...
#
#    clust_dat = np.zeros((how_many), dtype = np.ndarray)
#    clust_dat_A = np.zeros((how_many), dtype = np.ndarray)
#    clust_dat_B = np.zeros((how_many), dtype = np.ndarray)
#    numerator_A = 0
#    denominator_tot = 0
#    
#    for m in range(0, how_many):
#        clust_dat_A[m] = (A_ids == m).sum()             # sum all A type particles in a cluster
#        clust_dat_B[m] = (B_ids == m).sum()
#        clust_dat[m] = clust_dat_A[m] + clust_dat_B[m]  # find total number of particles in cluster
#        if clust_dat[m] > 15:
#            numerator_A += clust_dat_A[m]
#            denominator_tot += clust_dat[m]
#    # get the total percent of A particles in all clusters
#    if denominator_tot != 0:
#        percent_A[j] =  float(numerator_A) / float(denominator_tot)
#
#    l_clust = 0                                             # int size of largest cluster
#    for k in range(0, len(size_clusters[j])):
#        # the size minimum is a very important value to consider
#        if size_clusters[j][k] > 15 and size_clusters[j][k] < part_num:
#            tot_size[j] += size_clusters[j][k]
#            tot_num[j] += 1
#            if size_clusters[j][k] > l_clust:           # if larger cluster is found
#                l_clust = size_clusters[j][k]           # set l_clust to that size
#                    
#    largest[j] = l_clust                                # save largest cluster size for tstep
#
#    if tot_num[j] > 0:
#        MCS[j] = float(tot_size[j]/tot_num[j])/float(part_num)
#        GF[j] = float(part_num - tot_size[j]) / float(part_num)
#    
#    else:
#        MCS[j] = 0
#        GF[j] = 1
#
#    # let's start by getting the MSD for all particles (don't care about type)
#    if j != dumps - 1:
#        msda_count = 0
#        msdb_count = 0
#        for w in range(0,part_num):
#            MSD[j][w] = np.sqrt(((position_array[j+1][w][0] - position_array[j][w][0])**2) +
#                                ((position_array[j+1][w][1] - position_array[j][w][1])**2) +
#                                ((position_array[j+1][w][2] - position_array[j][w][2])**2))
#                                
#            if type_array[j][w] == 0:
#                MSD_A[j][msda_count] = np.sqrt(((position_array[j+1][w][0] - position_array[j][w][0])**2) +
#                                               ((position_array[j+1][w][1] - position_array[j][w][1])**2) +
#                                               ((position_array[j+1][w][2] - position_array[j][w][2])**2))
#                msda_count += 1
#            else:
#                MSD_B[j][msdb_count] = np.sqrt(((position_array[j+1][w][0] - position_array[j][w][0])**2) +
#                                               ((position_array[j+1][w][1] - position_array[j][w][1])**2) +
#                                               ((position_array[j+1][w][2] - position_array[j][w][2])**2))
#                msdb_count += 1
#
#
#
#def getDensityPlease(n):                                # call this function as needed
#    l_pos = position_array[n]                           # get ith position array
#    my_density.compute(f_box,
#                       l_pos,
#                       l_pos)
#    return my_density.getDensity()
#
#avg_sys_density = np.zeros((1), dtype=np.ndarray)
#
#take_last = dumps - 50
#last = dumps - 1
#msd_last = dumps - 2
#for j in range(take_last, dumps):
#    avg_sys_density[0] += getDensityPlease(j)
#
#avg_sys_density[0] /= (dumps - take_last)
#
##########################################################################
#### perform the same analysis on species A and species B individually ###
##########################################################################
#
#if part_perc_a != 0 and part_perc_a != 100:
#
#    tot_size_A = np.zeros((dumps), dtype=np.ndarray)          # number of particles in clusters
#    tot_num_A = np.zeros((dumps), dtype=np.ndarray)           # total number of clusters
#    MCS_A = np.zeros((dumps), dtype=np.ndarray)               # Mean cluster size
#    GF_A = np.zeros((dumps), dtype=np.ndarray)                # Gas fraction
#
#    tot_size_B = np.zeros((dumps), dtype=np.ndarray)          # number of particles in clusters
#    tot_num_B = np.zeros((dumps), dtype=np.ndarray)           # total number of clusters
#    MCS_B = np.zeros((dumps), dtype=np.ndarray)               # Mean cluster size
#    GF_B = np.zeros((dumps), dtype=np.ndarray)                # Gas fraction
#
#    for j in range(0, dumps):
#        countA = 0
#        countB = 0
#        for g in range(0, part_num):
#            if type_array[j][g] == 0:
#                tmpA[countA][0] = position_array[j][g][0]
#                tmpA[countA][1] = position_array[j][g][1]
#                tmpA[countA][2] = position_array[j][g][2]
#                countA += 1
#            else:
#                tmpB[countB][0] = position_array[j][g][0]
#                tmpB[countB][1] = position_array[j][g][1]
#                tmpB[countB][2] = position_array[j][g][2]
#                countB += 1
#
#        pos_A[j] = tmpA
#        pos_B[j] = tmpB
#        
#        l_pos = pos_A[j]
#        my_clusters.computeClusters(l_pos)
#        number_clusters[j] = my_clusters.getNumClusters()   # find number of clusters
#        ids = my_clusters.getClusterIdx()                   # get cluster ids
#        cluster_props.computeProperties(l_pos, ids)
#        size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each
#        
#        for k in range(0, len(size_clusters[j])):
#            # the size minimum is a very important value to consider
#            if size_clusters[j][k] > 15 and size_clusters[j][k] < part_num:
#                tot_size_A[j] += size_clusters[j][k]
#                tot_num_A[j] += 1
#
#        if tot_num_A[j] > 0:
#            MCS_A[j] = float(tot_size_A[j]/tot_num_A[j])/float(part_a)
#            GF_A[j] = float(part_a - tot_size_A[j]) / float(part_a)
#        
#        else:
#            MCS_A[j] = 0
#            GF_A[j] = 1
#
#        l_pos = pos_B[j]
#        my_clusters.computeClusters(l_pos)
#        number_clusters[j] = my_clusters.getNumClusters()   # find number of clusters
#        ids = my_clusters.getClusterIdx()                   # get cluster ids
#        cluster_props.computeProperties(l_pos, ids)
#        size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each
#        
#        for k in range(0, len(size_clusters[j])):
#            # the size minimum is a very important value to consider
#            if size_clusters[j][k] > 15 and size_clusters[j][k] < part_num:
#                tot_size_B[j] += size_clusters[j][k]
#                tot_num_B[j] += 1
#
#        if tot_num_B[j] > 0:
#            MCS_B[j] = float(tot_size_B[j]/tot_num_B[j])/float(part_b)
#            GF_B[j] = float(part_b - tot_size_B[j]) / float(part_b)
#        
#        else:
#            MCS_B[j] = 0
#            GF_B[j] = 1
#
#
#
#    def getDensityA(n):                                     # call this function as needed
#        countA = 0
#        for g in range(0, part_num):
#            if type_array[n][g] == 0:
#                tmpA[countA][0] = position_array[n][g][0]
#                tmpA[countA][1] = position_array[n][g][1]
#                tmpA[countA][2] = position_array[n][g][2]
#                countA += 1
#        pos_A[n] = tmpA
#        l_pos = pos_A[n]                                    # get ith position array
#        my_density.compute(f_box,
#                           l_pos,
#                           l_pos)
#        return my_density.getDensity()
#
#    avg_dense_A = np.zeros((1), dtype=np.ndarray)
#
#    for j in range(take_last, dumps):
#        avg_dense_A[0] += getDensityA(j)
#
#    avg_dense_A[0] /= (dumps - take_last)
#
#    def getDensityB(n):                                     # call this function as needed
#        countB = 0
#        for g in range(0, part_num):
#            if type_array[n][g] == 1:
#                tmpB[countB][0] = position_array[n][g][0]
#                tmpB[countB][1] = position_array[n][g][1]
#                tmpB[countB][2] = position_array[n][g][2]
#                countB += 1
#        pos_B[n] = tmpB
#        l_pos = pos_B[n]                                    # get ith position array
#        my_density.compute(f_box,
#                           l_pos,
#                           l_pos)
#        return my_density.getDensity()
#
#    avg_dense_B = np.zeros((1), dtype=np.ndarray)
#
#    for j in range(take_last, dumps):
#        avg_dense_B[0] += getDensityB(j)
#
#    avg_dense_B[0] /= (dumps - take_last)

##############################################
##### Plot the individual and total data #####
##############################################

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes=True)

plt_name  = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a)
plt_name1 = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + "A"
plt_name2 = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + "B"

if part_perc_a != 0 and part_perc_a != 100:
#    sns.kdeplot(avg_sys_density[0], shade = True, color="g")
#    sns.kdeplot(avg_dense_A[0], shade = True, color="r")
#    sns.kdeplot(avg_dense_B[0], shade = True, color="b")
#    plt.savefig('avg_density_' + plt_name + '.png', dpi=1000)
#    plt.close()
#
#    sns.kdeplot(getDensityPlease(last), shade = True, color="g")
#    sns.kdeplot(getDensityA(last), shade = True, color="r")
#    sns.kdeplot(getDensityB(last), shade = True, color="b")
#    plt.savefig('final_density_' + plt_name + '.png', dpi=1000)
#    plt.close()
#
#    plt.plot(MCS, color="g")
#    plt.plot(MCS_A, color="r")
#    plt.plot(MCS_B, color="b")
#    #plt.ylim((0,1))
#    plt.savefig('MCS_'+ plt_name + '.png', dpi=1000)
#    plt.close()
#
#    plt.plot(GF, color="g")
#    plt.plot(GF_A, color="r")
#    plt.plot(GF_B, color="b")
#    plt.ylim((0,1))
#    plt.savefig('GF_'+plt_name+'.png', dpi=1000)
#    plt.close()
#
#    plt.plot(percent_A, color="r")
#    #plt.ylim((0,1))
#    plt.savefig('A_comp_'+plt_name+'.png', dpi=1000)
#    plt.close()
#
#    plt.plot(largest, color="g")
#    plt.savefig('Largest_clust_'+plt_name+'.png', dpi=1000)
#    plt.close()
#
#    sns.kdeplot(MSD[msd_last], shade = True, color="g")
#    sns.kdeplot(MSD_A[msd_last], shade = True, color="r")
#    sns.kdeplot(MSD_B[msd_last], shade = True, color="b")
#    plt.savefig('MSD_'+plt_name+'.png', dpi=1000)
#    plt.close()
    plt.plot(struct_fact, color="g")
    plt.savefig('Structure_factor_'+plt_name+'.png', dpi=1000)
    plt.close()

else:                                                           # if monodisperse plot total values
#    sns.kdeplot(avg_sys_density[0], shade = True, color="g")
#    plt.savefig('avg_density_' + plt_name + '.png', dpi=1000)
#    plt.close()
#
#    sns.kdeplot(getDensityPlease(last), shade = True, color="g")
#    plt.savefig('final_density_' + plt_name + '.png', dpi=1000)
#    plt.close()
#
#    plt.plot(MCS, color="g")
#    plt.savefig('MCS_'+ plt_name + '.png', dpi=1000)
#    plt.close()
#
#    plt.plot(GF, color="g")
#    plt.ylim((0,1))
#    plt.savefig('GF_'+plt_name+'.png', dpi=1000)
#    plt.close()
#
#    plt.plot(largest, color="g")
#    plt.savefig('Largest_clust_'+plt_name+'.png', dpi=1000)
#    plt.close()
#
#    sns.kdeplot(MSD[msd_last], shade = True, color="g")
#    plt.savefig('MSD_'+plt_name+'.png', dpi=1000)
#    plt.close()

    plt.plot(struct_fact, color="g")
    plt.savefig('Structure_factor_'+plt_name+'.png', dpi=1000)
    plt.close()
