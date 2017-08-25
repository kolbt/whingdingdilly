import sys
import time
sys.path.append("/Users/kolbt/Desktop")
from wdd_pkg import my_clust
from wdd_pkg import qclust
from wdd_pkg import msd
from wdd_pkg import msd_out
from wdd_pkg import pos_by_type
from wdd_pkg import nearest_neighbor
from mypy import hey_there
from mypy import load_bar

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

myfile = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"
msdfile = "MSD_pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"

f = hoomd.open(name=myfile, mode='rb')
dumps = f.__len__()
size_min = 1000                                         # minimum size of cluster

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


#ids = np.zeros((dumps), dtype=np.ndarray)               # arrays to store things
size_clusters = np.zeros((dumps), dtype=np.ndarray)
tot_size = np.zeros((dumps), dtype=np.ndarray)          # number of particles in clusters
tot_num = np.zeros((dumps), dtype=np.ndarray)           # total number of clusters
MCS = np.zeros((dumps), dtype=np.ndarray)               # Mean cluster size
GF = np.zeros((dumps), dtype=np.ndarray)                # Gas fraction
percent_A = np.zeros((dumps), dtype=np.ndarray)         # composition A at each timestep
largest = np.zeros((dumps), dtype=np.ndarray)           # read out largest cluster at each tstep

avg_aa = np.zeros((dumps), dtype=np.float32)
avg_ab = np.zeros((dumps), dtype=np.float32)
avg_ba = np.zeros((dumps), dtype=np.float32)
avg_bb = np.zeros((dumps), dtype=np.float32)

LIQ_A = np.zeros((dumps - 1), dtype=np.ndarray)         # arrays for MSD
LIQ_B = np.zeros((dumps - 1), dtype=np.ndarray)
GAS_A = np.zeros((dumps - 1), dtype=np.ndarray)
GAS_B = np.zeros((dumps - 1), dtype=np.ndarray)
MSD_T = np.zeros((dumps - 1), dtype=np.float64)
MSD_TL = np.zeros((dumps - 1), dtype=np.ndarray)
MSD_TG = np.zeros((dumps - 1), dtype=np.ndarray)

disp_x = np.zeros((part_num), dtype=np.ndarray)         # displacement vectors
disp_y = np.zeros((part_num), dtype=np.ndarray)
disp_z = np.zeros((part_num), dtype=np.ndarray)

################################################################################################
################################# Actual analysis starts here ##################################
################################################################################################

# analyze all particles
load_bar.printLoadBar(0, dumps, prefix = "Progress:", suffix = "Complete")
for j in range(0, dumps):
    
    load_bar.printLoadBar(j+1, dumps, prefix = "Progress:", suffix = "Complete")
    
    my_clusters.computeClusters(position_array[j])
    ids = my_clusters.getClusterIdx()                   # get cluster ids
    cluster_props.computeProperties(position_array[j], ids)
    size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each
    
    how_many = my_clusters.getNumClusters()
    f_largest = "largest_pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".txt"

    l_clust,\
    MCS[j],\
    GF[j] = my_clust.clustDat(size_clusters[j],
                              part_num,
                              j,
                              f_largest,
                              size_min,
                              True)
                    
    largest[j] = l_clust                                # save largest cluster size for tstep

    #########################################################
    ### Find MSD for A, B individually, also total system ###
    #########################################################

    q_clust = qclust.binaryCluster(ids, how_many)

    if j > 0:
        
        msd_val = msd.msd(l_box,
                          position_array[j],
                          position_array[j-1],
                          disp_x,
                          disp_y,
                          disp_z)
                          
        LIQ_A[j-1],\
        LIQ_B[j-1],\
        GAS_A[j-1],\
        GAS_B[j-1],\
        MSD_TG[j-1],\
        MSD_TL[j-1],\
        percent_A[j] = msd_out.writeArrays(msd_val,
                                           q_clust,
                                           ids,
                                           type_array[j])
        MSD_T[j-1] = np.sum(msd_val)

    if j == dumps - 1:
        avg_aa[j],\
        avg_ab[j],\
        avg_ba[j],\
        avg_bb[j] = nearest_neighbor.nearNeigh(q_clust,
                                               ids,
                                               position_array[j],
                                               type_array[j])

############################
### Density caluclations ###
############################

# OUT OF LOOP

avg_sys_density = 0
last_dense = 0
take_last = dumps - 10
last = dumps - 1
msd_last = dumps - 2
for j in range(take_last, dumps):
    my_density.compute(f_box, position_array[j], position_array[j])
    avg_sys_density += my_density.getDensity()

avg_sys_density /= (dumps - take_last)

my_density.compute(f_box, position_array[last], position_array[last])
last_dense += my_density.getDensity()

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

    load_bar.printLoadBar(0, dumps, prefix = "Progress:", suffix = "Complete")
    for j in range(0, dumps):
        load_bar.printLoadBar(j+1, dumps, prefix = "Progress:", suffix = "Complete")
        
        pos_A[j],\
        pos_B[j] = pos_by_type.posByType(position_array[j],
                                         type_array[j],
                                         part_a)
        
        ####################################
        ### GF, MCS for A-A correlations ###
        ####################################
        
        my_clusters.computeClusters(pos_A[j])
        ids = my_clusters.getClusterIdx()                   # get cluster ids
        cluster_props.computeProperties(pos_A[j], ids)
        size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each
        
        l_clust,\
        MCS_A[j],\
        GF_A[j] = my_clust.clustDat(size_clusters[j],
                                    part_num,
                                    j,
                                    f_largest,
                                    size_min,
                                    False)

        ####################################
        ### GF, MCS for B-B correlations ###
        ####################################

        my_clusters.computeClusters(pos_B[j])
        ids = my_clusters.getClusterIdx()                   # get cluster ids
        cluster_props.computeProperties(pos_B[j], ids)
        size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each
        
        l_clust,\
        MCS_B[j],\
        GF_B[j] = my_clust.clustDat(size_clusters[j],
                                    part_num,
                                    j,
                                    f_largest,
                                    size_min,
                                    False)

    avg_dense_A = 0
    last_denseA = 0
    for j in range(take_last, dumps):
        my_density.compute(f_box, position_array[j], pos_A[j])
        avg_dense_A += my_density.getDensity()

    avg_dense_A /= (dumps - take_last)
    my_density.compute(f_box, position_array[last], pos_A[last])
    last_denseA += my_density.getDensity()

    avg_dense_B = 0
    last_denseB = 0
    for j in range(take_last, dumps):
        my_density.compute(f_box, position_array[j], pos_B[j])
        avg_dense_B += my_density.getDensity()

    avg_dense_B /= (dumps - take_last)
    my_density.compute(f_box, position_array[last], pos_B[last])
    last_denseB += my_density.getDensity()

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
    sns.kdeplot(avg_sys_density, shade = True, color="g")
    sns.kdeplot(avg_dense_A, shade = True, color="r")
    sns.kdeplot(avg_dense_B, shade = True, color="b")
    plt.savefig('avg_density_' + plt_name + '.png', dpi=1000)
    plt.close()

    sns.kdeplot(last_dense, shade = True, color="g", label='all')
    sns.kdeplot(last_denseA, shade = True, color="r", label='A')
    sns.kdeplot(last_denseB, shade = True, color="b", label='B')
    plt.legend(loc='upper left')
    plt.savefig('final_density_' + plt_name + '.png', dpi=1000)
    plt.close()

    plt.plot(MCS, color="g")
    plt.plot(MCS_A, color="r")
    plt.plot(MCS_B, color="b")
    #plt.ylim((0,1))
    plt.savefig('MCS_'+ plt_name + '.png', dpi=1000)
    plt.close()

    plt.plot(GF, color="g")
    plt.plot(GF_A, color="r")
    plt.plot(GF_B, color="b")
    plt.ylim((0,1))
    plt.savefig('GF_'+plt_name+'.png', dpi=1000)
    plt.close()

    plt.plot(percent_A, color="r")
    #plt.ylim((0,1))
    plt.savefig('A_comp_'+plt_name+'.png', dpi=1000)
    plt.close()

    plt.plot(largest, color="g")
    plt.savefig('Largest_clust_'+plt_name+'.png', dpi=1000)
    plt.close()

    plt.plot(msd_time, GAS_A,  color="r", marker='o', markersize=1, linestyle='None', label='Gas_A')
    plt.plot(msd_time, GAS_B,  color="b", marker='o', markersize=1, linestyle='None', label='Gas_B')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Timesteps')
    plt.ylabel('MSD')
    plt.legend(loc='upper left')
    plt.savefig('MSD_GAS_AB_' + plt_name + '.png', dpi=1000)
    plt.close()
    
    #if np.any(LIQ_A) == 1 and np.any(LIQ_B) == 1:
    plt.plot(msd_time, LIQ_A,  color="r", marker='o', markersize=1, linestyle='None', label='Liq_A')
    plt.plot(msd_time, LIQ_B,  color="b", marker='o', markersize=1, linestyle='None', label='Liq_B')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Timesteps')
    plt.ylabel('MSD')
    plt.legend(loc='upper left')
    plt.savefig('MSD_LIQ_AB_' + plt_name + '.png', dpi=1000)
    plt.close()

    plt.plot(msd_time, MSD_T,  color="g", marker='o', markersize=1, linestyle='None', label='MSD')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Timesteps')
    plt.ylabel('MSD')
    plt.legend(loc='upper left')
    plt.savefig('MSD_total_' + plt_name + '.png', dpi=1000)
    plt.close()

    #if np.any(MSD_TL) == 1:
    plt.plot(msd_time, MSD_TL,  color="b", marker='o', markersize=1, linestyle='None', label='Liq')
    plt.plot(msd_time, MSD_TG,  color="r", marker='o', markersize=1, linestyle='None', label='Gas')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Timesteps')
    plt.ylabel('MSD')
    plt.legend(loc='upper left')
    plt.savefig('MSD_LG_' + plt_name + '.png', dpi=1000)
    plt.close()

    plt.plot(avg_aa, color="g", label='avg_aa')
    plt.plot(avg_ab, color="r", label='avg_ab')
    plt.plot(avg_ba, color="b", label='avg_ba')
    plt.plot(avg_bb, color="k", label='avg_bb')
    plt.legend(loc='upper left')
    plt.savefig('avg_neighs_'+ plt_name + '.png', dpi=1000)
    plt.close()

else:                                                           # if monodisperse plot total values
    sns.kdeplot(avg_sys_density[0], shade = True, color="g")
    plt.savefig('avg_density_' + plt_name + '.png', dpi=1000)
    plt.close()

    sns.kdeplot(getDensityPlease(last), shade = True, color="g")
    plt.savefig('final_density_' + plt_name + '.png', dpi=1000)
    plt.close()

    plt.plot(MCS, color="g")
    plt.savefig('MCS_'+ plt_name + '.png', dpi=1000)
    plt.close()

    plt.plot(GF, color="g")
    plt.ylim((0,1))
    plt.savefig('GF_'+plt_name+'.png', dpi=1000)
    plt.close()

    plt.plot(largest, color="g")
    plt.savefig('Largest_clust_'+plt_name+'.png', dpi=1000)
    plt.close()

    plt.plot(msd_time, MSD_T,  color="g", marker='o', markersize=1, linestyle='None', label='MSD')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Timesteps')
    plt.ylabel('MSD')
    plt.legend(loc='upper left')
    plt.savefig('MSD_total_' + plt_name + '.png', dpi=1000)
    plt.close()

    if np.any(MSD_TL) == 1:
        plt.plot(msd_time, MSD_TL,  color="b", marker='o', markersize=1, linestyle='None', label='Liq')
    plt.plot(msd_time, MSD_TG,  color="r", marker='o', markersize=1, linestyle='None', label='Gas')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Timesteps')
    #plt.xlabel(r'Time ($\tau$)')
    plt.ylabel('MSD')
    plt.legend(loc='upper left')
    plt.savefig('MSD_LG_' + plt_name + '.png', dpi=1000)
    plt.close()
