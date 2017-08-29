#Arguments in order:
#    1.) hoomd_path (string): the path to hoomd on your machine
#    2.) tsteps (int): the number of timesteps your simulation is running
#    3.) dump_freq (int): the frequency of gsd dumps
#    4.) part_frac_a (float): particle fraction of type a particles
#    5.) pe_a (float): activity of type a particles
#    6.) pe_b (float): activity of tyep b particles

import sys
import os

#hoomd_path = str(sys.argv[1])
hoomd_path = "${hoomd_path}"
#tsteps = int(sys.argv[2])
tsteps = ${tsteps}
#dump_freq = int(sys.argv[3])
dump_freq = ${dump_freq}
#part_frac_a = float(sys.argv[4])
part_perc_a = ${part_frac_a}
part_frac_a = float(part_perc_a) / 100.0
#pe_a = float(sys.argv[5])
pe_a = ${pe_a}
#pe_b = float(sys.argv[6])
pe_b = ${pe_b}
#gsd_path = str(sys.argv[7])
gsd_path = "${gsd_path}"

# calculate number of tsteps which are dumped
dumps = tsteps/dump_freq

sys.path.append(hoomd_path)

import hoomd
from hoomd import md
from hoomd import deprecated
import numpy as np

#get tsteps for msd calculations
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
#msd_dumps += 9110000
msd_dumps += tsteps - 1000000 + 110000

#initialize system randomly, can specify GPU execution here

part_num = ${part_num}

hoomd.context.initialize()
system = hoomd.deprecated.init.create_random(N = part_num,
                                             phi_p = 0.6,
                                             name = 'A',
                                             min_dist = 0.5,
                                             seed = 230958,
                                             dimensions = 2)

system.particles.types.add('B')
snapshot = system.take_snapshot()

part_a = part_num * part_frac_a         # get the total number of A particles
part_a = int(part_a)
part_b = part_num - part_a              # get the total number of B particles
part_b = int(part_b)
mid = int(part_a)                       # starting point for assigning B particles

if part_perc_a == 0:                    # take care of all b case
    mid = 0
    for i in range(mid,part_num):
        system.particles[i].type = 'B'
elif part_perc_a != 100:                # mix of each
    for i in range(mid,part_num):
        system.particles[i].type = 'B'

all = hoomd.group.all()
gA = hoomd.group.type(type = 'A', update=True)
gB = hoomd.group.type(type = 'B', update=True)
N = len(all)
Na = len(gA)
Nb = len(gB)

#define potential between pairs
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=2**(1/6), nlist=nl)
lj.set_params(mode='shift')
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('A', 'B', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0)

#integrator type
hoomd.md.integrate.mode_minimize_fire(group=all, dt=0.00001, ftol=1e-2, Etol=1e-7)
hoomd.run(10000)

#run simulation with current settings here
#hoomd.md.integrate.mode_standard(dt=0.000005)
hoomd.md.integrate.mode_standard(dt=0.0000005)
hoomd.md.integrate.brownian(group=all, kT=0.5, seed=123)
hoomd.run(100000)

#set the activity of each type

angle = np.random.rand(part_num) * 2 * np.pi

if part_perc_a != 0 and part_perc_a != 100:
    activity_a = []
    for i in range(0,mid):
        x = (np.cos(angle[i])) * pe_a
        y = (np.sin(angle[i])) * pe_a
        z = 0
        tuple = (x, y, z)
        activity_a.append(tuple)
    activity_b = []
    for i in range(mid,part_num):
        x = (np.cos(angle[i])) * pe_b
        y = (np.sin(angle[i])) * pe_b
        z = 0
        tuple = (x, y, z)
        activity_b.append(tuple)
    hoomd.md.force.active(group=gA,
                          seed=123,
                          f_lst=activity_a,
                          rotation_diff=3.0,
                          orientation_link=False)
    hoomd.md.force.active(group=gB,
                          seed=375,
                          f_lst=activity_b,
                          rotation_diff=3.0,
                          orientation_link=False)
else:
    if part_perc_a == 0:
        activity_b = []
        for i in range(0,part_num):
            x = (np.cos(angle[i])) * pe_b
            y = (np.sin(angle[i])) * pe_b
            z = 0
            tuple = (x, y, z)
            activity_b.append(tuple)
        hoomd.md.force.active(group=gB,
                              seed=375,
                              f_lst=activity_b,
                              rotation_diff=3.0,
                              orientation_link=False)
    else:
        activity_a = []
        for i in range(0,part_num):
            x = (np.cos(angle[i])) * pe_a
            y = (np.sin(angle[i])) * pe_a
            z = 0
            tuple = (x, y, z)
            activity_a.append(tuple)
        hoomd.md.force.active(group=gA,
                              seed=123,
                              f_lst=activity_a,
                              rotation_diff=3.0,
                              orientation_link=False)

#write dumps
name = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"
msd_name = "MSD_pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"

### Dump for MSD ###
def dump_spec(timestep):
    if timestep in msd_dumps:
        hoomd.dump.gsd(filename=msd_name, period=None, group=all, overwrite=False, static=[])
        os.close(2)

hoomd.analyze.callback(callback = dump_spec, period = 1)
####################

hoomd.dump.gsd(name, period=dump_freq, group=all, overwrite=True, static=[])

#run
hoomd.run(tsteps)

#########################################################################
########################## Begin Data Analysis ##########################
#########################################################################

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
    number_clusters[j] = my_clusters.getNumClusters()   # find number of clusters
    ids = my_clusters.getClusterIdx()                   # get cluster ids
    cluster_props.computeProperties(l_pos, ids)
    size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each
    
    how_many = my_clusters.getNumClusters()
    
    #############################################################
    ### This finds the cluster ids for type A and B particles ###
    #############################################################
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

    #######################################################################
    ### If clusters are greater than a threshold size, find composition ###
    #######################################################################
    
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
    
    if tot_num[j] > 0:
        MCS[j] = float(tot_size[j]/tot_num[j])/float(part_num)
        GF[j] = float(part_num - tot_size[j]) / float(part_num)
    else:
        MCS[j] = 0
        GF[j] = 1

    #########################################################
    ### Find MSD for A, B individually, also total system ###
    #########################################################

    # ?you can enhance difference between gas and liq by setting min clust requirement

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

############################
### Density caluclations ###
############################

def getDensityPlease(n):                                # call this function as needed
    l_pos = position_array[n]                           # get ith position array
    my_density.compute(f_box,
                       l_pos,
                       l_pos)
    return my_density.getDensity()

avg_sys_density = np.zeros((1), dtype=np.ndarray)

take_last = dumps - 10
last = dumps - 1
msd_last = dumps - 2
for j in range(take_last, dumps):
    avg_sys_density[0] += getDensityPlease(j)

avg_sys_density[0] /= (dumps - take_last)

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
        number_clusters[j] = my_clusters.getNumClusters()   # find number of clusters
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
            MCS_A[j] = 0
            GF_A[j] = 1

        l_pos = pos_B[j]
        my_clusters.computeClusters(l_pos)
        number_clusters[j] = my_clusters.getNumClusters()   # find number of clusters
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
            MCS_B[j] = 0
            GF_B[j] = 1



    def getDensityA(n):                                     # call this function as needed
        countA = 0
        for g in range(0, part_num):
            if type_array[n][g] == 0:
                tmpA[countA][0] = position_array[n][g][0]
                tmpA[countA][1] = position_array[n][g][1]
                tmpA[countA][2] = position_array[n][g][2]
                countA += 1
        pos_A[n] = tmpA
        l_pos = pos_A[n]                                    # get ith position array
        my_density.compute(f_box,
                           l_pos,
                           l_pos)
        return my_density.getDensity()

    avg_dense_A = np.zeros((1), dtype=np.ndarray)
    
    for j in range(take_last, dumps):
        avg_dense_A[0] += getDensityA(j)

    avg_dense_A[0] /= (dumps - take_last)

    def getDensityB(n):                                     # call this function as needed
        countB = 0
        for g in range(0, part_num):
            if type_array[n][g] == 1:
                tmpB[countB][0] = position_array[n][g][0]
                tmpB[countB][1] = position_array[n][g][1]
                tmpB[countB][2] = position_array[n][g][2]
                countB += 1
        pos_B[n] = tmpB
        l_pos = pos_B[n]                                    # get ith position array
        my_density.compute(f_box,
                           l_pos,
                           l_pos)
        return my_density.getDensity()

    avg_dense_B = np.zeros((1), dtype=np.ndarray)
    
    for j in range(take_last, dumps):
        avg_dense_B[0] += getDensityB(j)
    
    avg_dense_B[0] /= (dumps - take_last)

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
    sns.kdeplot(avg_sys_density[0], shade = True, color="g")
    sns.kdeplot(avg_dense_A[0], shade = True, color="r")
    sns.kdeplot(avg_dense_B[0], shade = True, color="b")
    plt.savefig('avg_density_' + plt_name + '.png', dpi=1000)
    plt.close()
    
    sns.kdeplot(getDensityPlease(last), shade = True, color="g")
    sns.kdeplot(getDensityA(last), shade = True, color="r")
    sns.kdeplot(getDensityB(last), shade = True, color="b")
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
