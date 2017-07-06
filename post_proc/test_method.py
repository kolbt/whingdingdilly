import sys

#hoomd_path = str(sys.argv[4])
#gsd_path = str(sys.argv[5])

# need to extract values from filename (pa, pb, xa) for naming
#part_perc_a = int(sys.argv[3])
#part_frac_a = float(part_perc_a) / 100.0
#pe_a = int(sys.argv[1])
#pe_b = int(sys.argv[2])

# manual input
hoomd_path = "/Users/kolbt/Desktop/compiled/hoomd-blue/build"
gsd_path = "/Users/kolbt/Desktop/compiled/gsd/build"

part_perc_a = 50
part_frac_a = float(part_perc_a) / 100.0
pe_a = 80
pe_b = 80

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
#myfile = "test.gsd"

f = hoomd.open(name=myfile, mode='rb')
dumps = f.__len__()
print(dumps)

real_time = np.arange(0,float(dumps))
#real_time *= 0.00001 * 20000                           # step_size * frequency
real_time *= 0.2

log_time = np.zeros((dumps-1), dtype=np.float64)
#value_to_dump = 0
#jumper = 1
#count = 1
#for jjj in range(1,63):
#    if (count-2) % 9 == 0 and count != 2:
#        jumper *= 10
#    value_to_dump += jumper
#    log_time[jjj] = value_to_dump
#    count += 1

#jumper = 10
#value_to_dump = 100
#count = 0
#for iii in range(1,dumps):
#    if iii < 100:
#        log_time[iii-1] = iii
#    else:
#        log_time[iii-1] = value_to_dump
#        value_to_dump += jumper
#        count += 1
#    if count % 90 == 0 and count != 0:
#        jumper *= 10
#        count = 0
#    print(log_time[iii-1])

msd_dumps = np.zeros((101), dtype=np.float64)
jumper = 5
value_to_dump = 15
count = 10
for iii in range(0,101):
    if iii <= 10:
        msd_dumps[iii] = iii
    elif count == 95:
        msd_dumps[iii] = value_to_dump
        jumper *= 10
        value_to_dump += jumper
        count = 10
    else:
        msd_dumps[iii] = value_to_dump
        value_to_dump += jumper
        count += 5

msd_dumps += 9110000



log_time *= 0.00001

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
                              rcut=0.95)                 # initialize class
cluster_props = cluster.ClusterProperties(box=f_box)

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
    
    print("On timestep " + str(j))
    
    l_pos = position_array[j]
    my_clusters.computeClusters(l_pos)
    number_clusters[j] = my_clusters.getNumClusters()   # find number of clusters
    ids = my_clusters.getClusterIdx()                   # get cluster ids
    cluster_props.computeProperties(l_pos, ids)
    size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each
    
    how_many = my_clusters.getNumClusters()
    #print(how_many)

    sort_id = np.sort(ids)                              # IDs sorted small to large, len = part_num
    q_clust = np.zeros((how_many), dtype=np.ndarray)    # my binary 'is it clustered?' array
    index = 0                                           # index of the sorted array to look at
    size_min = 15                                       # ignore clusters that aren't at least 15
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
            
            msd_val = (((disp_x[b])**2) + ((disp_y[b])**2) + ((disp_z[b])**2))
            MSD_T[j-1] += msd_val
            if q_clust[ids[b]] == 1:                        # check if in liquid
                MSD_TL[j-1] += msd_val                        # add to tot. lq. msd
                if type_array[j][b] == 0:                   # type A case
                    LIQ_A[j-1] += msd_val
                    lq_a_count += 1
                else:
                    LIQ_B[j-1] += msd_val
                    lq_b_count += 1
            else:                                           # else, particle is gas
                MSD_TG[j-1] += msd_val                        # add to tot. gs. msd
                if type_array[j][b] == 0:                   # type A case
                    GAS_A[j-1] += msd_val
                    gs_a_count += 1
                else:
                    GAS_B[j-1] += msd_val
                    gs_b_count += 1

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
        print(percent_A[j])

#print(MSD_T[j])

# take log of values?
#for w in range(0,len(LIQ_A)):
#    LIQ_A[w] = np.log10(LIQ_A[w])
#    LIQ_B[w] = np.log10(LIQ_B[w])
#    GAS_A[w] = np.log10(GAS_A[w])
#    GAS_B[w] = np.log10(GAS_B[w])
#    MSD_T[w] = np.log10(MSD_T[w])
#    MSD_TL[w] = np.log10(MSD_TL[w])
#    MSD_TG[w] = np.log10(MSD_TG[w])
#
#for y in range(0,len(real_time)):
#    real_time[y] = np.log10(real_time[y])


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

#########################################################################
### perform the same analysis on species A and species B individually ###
#########################################################################

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
from scipy import interpolate
#from scipy.signal import savgol_filter
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

#    zzz = savgol_filter(GAS_A, 41, 12)

#    plt.plot(real_time, MSD_T, color="g")
#    plt.plot(real_time, MSD_TL, color="r")
#    plt.plot(real_time, MSD_TG, color="b")
#    plt.savefig('MSD_TS' + plt_name + '.png', dpi=1000)
#    plt.close()
#    
#    plt.plot(real_time, LIQ_A, color="r")
#    plt.savefig('MSD_LA' + plt_name + '.png', dpi=1000)
#    plt.close()
#    
#    plt.plot(real_time, LIQ_B, color="b")
#    plt.savefig('MSD_LB' + plt_name + '.png', dpi=1000)
#    plt.close()
#    
#    plt.plot(real_time, GAS_A, color="r")
#    #    plt.plot(real_time, zzz, basex=10, color="g")
#    plt.savefig('MSD_GA' + plt_name + '.png', dpi=1000)
#    plt.close()
#    
#    plt.plot(real_time, GAS_B, color="b")
#    plt.savefig('MSD_GB' + plt_name + '.png', dpi=1000)
#    plt.close()

    plt.plot(msd_time, MSD_T,  color="g", marker='o', markersize=1, linestyle='None', label='MSD')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Time ($\tau$)')
    plt.ylabel('MSD')
    plt.legend(loc='upper left')
    plt.savefig('FITTED' + plt_name + '.png', dpi=1000)
    plt.close()
    
#    x = log_time
#    y = MSD_T
#    
#    x2 = np.zeros((100), dtype=np.float64)
#    y2 = np.zeros((100), dtype=np.float64)
#    
#    x3 = np.zeros((dumps - 101), dtype=np.float64)
#    y3 = np.zeros((dumps - 101), dtype=np.float64)
#    
#    for i in range(0,100):
#        x2[i] = x[i]
#        y2[i] = y[i]
#    for i in range(100,dumps-1):
#        x3[i-100] = x[i]
#        y3[i-100] = y[i]
#
#    first = np.polyfit(np.log10(x2), np.log10(y2), 1)
#    second = np.polyfit(np.log10(x3), np.log10(y3), 1)
##    poly1 = np.poly1d(first)
##    poly2 = np.poly1d(second)
##    lin1 = poly1(np.log10(x2))
##    lin2 = poly2(np.log10(x3))
#    print(first)
#    print(second)
#
#    def lin_plot(m, b, fx):
#        fy = (fx**m) * (10**b)
#        return fy
#
#    lin1 = lin_plot(first[0], first[1], x2)
#    lin2 = lin_plot(second[0], second[1], x3)
#
#    plt.plot(log_time, MSD_T,  color="g", marker='o', markersize=1, linestyle='None', label='MSD')
#    #plt.plot(x2, lin1, color="r", linestyle="solid", label='Slope = ' + str(first[0]))
#    #plt.plot(x3, lin2, color="b", linestyle="solid", label='Slope = ' + str(second[0]))
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.xlabel(r'Time ($\tau$)')
#    plt.ylabel('MSD')
#    plt.legend(loc='upper left')
#    plt.savefig('FITTED2' + plt_name + '.png', dpi=1000)
#    plt.close()

    
    
    
    
#    slopes = np.diff(np.log10(MSD_T)) / np.diff(np.log10(log_time))
#    slopes_diff = np.zeros((dumps-2), dtype=np.float64)
#    
#    for i in range(1,len(slopes)):
#        slopes_diff[i-1] = slopes[i] - slopes[i-1]
#    
#    plt.plot(slopes)
#    plt.ylim((0,2))
#    plt.savefig('SLOPES' + plt_name + '.png', dpi=1000)
#    plt.close()
#
#    plt.plot(slopes_diff)
#    plt.ylim((-2,2))
#    plt.savefig('SLOPES_DIFF' + plt_name + '.png', dpi=1000)
#    plt.close()


    

#
#    line1 = np.zeros((99), dtype=np.float64)
#    line2 = np.zeros((dumps-100), dtype=np.float64)
#
#    for i in range(1,100):
#        line1[i-1] = lin_plot(first[0], first[1], x[i])
#
#    for i in range(100,dumps):
#        line2[i-100] = lin_plot(second[0], second[1], x[i])

#    tck = interpolate.splrep(x2, y2, k=1, s=0)
#    tck2 = interpolate.splrep(x3, y3, k=1, s=0)
#    
#    y4 = interpolate.splev(x2, tck)
#    y5 = interpolate.splev(x3, tck2)

#    plt.loglog(log_time, MSD_T, basex=10, color="g")
#    plt.loglog(x2, line1, basex=10, color="r", linestyle="solid")
#    plt.loglog(x3, line2, basex=10, color="b", linestyle="solid")
#    plt.loglog(x2, line1, basex=10, color="r")
#    plt.loglog(x3, line2, basex=10, color="b")
    #plt.loglog(x, yinterp, basex=10, color="b", linestyle="solid")
#    plt.loglog(x2, y4, basex=10, color="b", linestyle="solid")
#    plt.loglog(x3, y5, basex=10, color="r", linestyle="solid")
#    plt.rcParams['backend'] = 'TkAgg'  # or, 'qt4agg'
#    plt.rcParams['agg.path.chunksize'] = 100000
    #plt.loglog(real_time, MSD_TL, basex=10, color="r")
    #plt.loglog(real_time, MSD_TG, basex=10, color="b")
#    plt.savefig('MSD_TS' + plt_name + '.png', dpi=1000)
#    plt.close()

#    plt.loglog(log_time, LIQ_A, basex=10, color="r")
#    plt.savefig('MSD_LA' + plt_name + '.png', dpi=1000)
#    plt.close()
#    
#    plt.loglog(log_time, LIQ_B, basex=10, color="b")
#    plt.savefig('MSD_LB' + plt_name + '.png', dpi=1000)
#    plt.close()
#    
#    plt.loglog(log_time, GAS_A, basex=10, color="r")
##    plt.loglog(real_time, zzz, basex=10, color="g")
#    plt.savefig('MSD_GA' + plt_name + '.png', dpi=1000)
#    plt.close()
#    
#    plt.loglog(log_time, GAS_B, basex=10, color="b")
#    plt.savefig('MSD_GB' + plt_name + '.png', dpi=1000)
#    plt.close()


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

    print("HI THERE")
