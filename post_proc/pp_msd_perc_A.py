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

#part_num = 15000
part_num = 500000

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

msdfile = "MSD_pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"

f = hoomd.open(name=msdfile, mode='rb')
dumps = f.__len__()
size_min = 1000                                         # minimum size of cluster

position_array = np.zeros((dumps), dtype=np.ndarray)    # array of position arrays
type_array = np.zeros((dumps), dtype=np.ndarray)        # particle types
box_data = np.zeros((1), dtype=np.ndarray)              # box dimensions
timesteps = np.zeros((dumps), dtype=np.float64)         # timesteps

with hoomd.open(name=msdfile, mode='rb') as t:          # open for reading
    snap = t[0]                                         # snap 0th snapshot
    box_data = snap.configuration.box                   # get box dimensions
    for i in range(0,dumps):
        snap = t[i]                                     # take snap of each dump
        type_array[i] = snap.particles.typeid
        position_array[i] = snap.particles.position     # store all particle positions
        timesteps[i] = snap.configuration.step          # store tstep for plotting purposes
part_num = len(type_array[0])
print(part_num)
timesteps -= timesteps[0]
msd_time = timesteps[1:]

from freud import parallel, box, density, cluster
parallel.setNumThreads(1)                               # don't run multiple threads

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
percent_A = np.zeros((dumps), dtype=np.ndarray)         # composition A at each timestep

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

# analyze all particles
for j in range(0, dumps):
    
    l_pos = position_array[j]
    my_clusters.computeClusters(l_pos)
    
    ids = my_clusters.getClusterIdx()                   # get cluster ids
    cluster_props.computeProperties(l_pos, ids)
    size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each
    
    how_many = my_clusters.getNumClusters()
    
    #########################################################
    ### Find MSD for A, B individually, also total system ###
    #########################################################

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

    lq_a_count = 0
    lq_b_count = 0
    gs_a_count = 0
    gs_b_count = 0
    if j > 0:
        numerator_A = 0
        denominator_tot = 0
        for b in range(0,part_num):
            
            # check instantaneous disp. over last timestep
            dx = position_array[j][b][0] - position_array[j-1][b][0]
            dy = position_array[j][b][1] - position_array[j-1][b][1]
            dz = position_array[j][b][2] - position_array[j-1][b][2]
            
            # if it is over some threshold, then it went past a boundary
            if dx < -50:
                dx += l_box
            if dx > 50:
                dx -= l_box
            disp_x[b] += dx
            
            if dy < -50:
                dy += l_box
            if dy > 50:
                dy -= l_box
            disp_y[b] += dy
            
            if dz < -50:
                dz += l_box
            if dz > 50:
                dz -= l_box
            disp_z[b] += dz
            
            msd_val = np.sqrt(((disp_x[b])**2) + ((disp_y[b])**2) + ((disp_z[b])**2))
            MSD_T[j-1] += msd_val
            if q_clust[ids[b]] == 1:                        # check if in liquid
                MSD_TL[j-1] += msd_val                      # add to tot. lq. msd
                if type_array[j][b] == 0:                   # type A case
                    LIQ_A[j-1] += msd_val
                    lq_a_count += 1
                else:
                    LIQ_B[j-1] += msd_val
                    lq_b_count += 1
            else:                                           # else, particle is gas
                MSD_TG[j-1] += msd_val                      # add to tot. gs. msd
                if type_array[j][b] == 0:                   # type A case
                    GAS_A[j-1] += msd_val
                    gs_a_count += 1
                else:
                    GAS_B[j-1] += msd_val
                    gs_b_count += 1

        # if-gating these so we don't break our program
        if lq_a_count != 0: LIQ_A[j-1] /= lq_a_count
        if lq_b_count != 0: LIQ_B[j-1] /= lq_b_count
        if gs_a_count != 0: GAS_A[j-1] /= gs_a_count
        if gs_b_count != 0: GAS_B[j-1] /= gs_b_count
        MSD_T[j-1] /= part_num
        if lq_a_count + lq_b_count != 0: MSD_TL[j-1] /= lq_a_count + lq_b_count
        if gs_a_count + gs_b_count != 0: MSD_TG[j-1] /= gs_a_count + gs_b_count

        numerator_A = lq_a_count
        denominator_tot = lq_a_count + lq_b_count
        
        if denominator_tot != 0:
            percent_A[j] =  float(numerator_A) / float(denominator_tot)

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

    plt.plot(percent_A, color="r")
    #plt.ylim((0,1))
    plt.savefig('A_comp_'+plt_name+'.png', dpi=1000)
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
    
    plt.plot(msd_time, MSD_TL,  color="b", marker='o', markersize=1, linestyle='None', label='Liq')
    plt.plot(msd_time, MSD_TG,  color="r", marker='o', markersize=1, linestyle='None', label='Gas')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Timesteps')
    plt.ylabel('MSD')
    plt.legend(loc='upper left')
    plt.savefig('MSD_LG_' + plt_name + '.png', dpi=1000)
    plt.close()

else:                                                           # if monodisperse plot total values

    plt.plot(msd_time, MSD_T,  color="g", marker='o', markersize=1, linestyle='None', label='MSD')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Timesteps')
    plt.ylabel('MSD')
    plt.legend(loc='upper left')
    plt.savefig('MSD_total_' + plt_name + '.png', dpi=1000)
    plt.close()

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
