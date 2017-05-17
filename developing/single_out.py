import sys

hoomd_path = "/Users/kolbt/Desktop/compiled/hoomd-blue/build"
tsteps = 500
dump_freq = 1
part_perc_a = 20
part_frac_a = float(part_perc_a) / 100.0
pe_a = 150
pe_b = 150
gsd_path = "/Users/kolbt/Desktop/compiled/gsd/build"

dumps = tsteps/dump_freq
part_num = 15000
part_a = int(float(part_num) * part_frac_a)
part_b = part_num - part_a

sys.path.append(gsd_path)
import gsd
from gsd import hoomd
from gsd import pygsd

import numpy as np

position_array = np.zeros((dumps), dtype=np.ndarray)    # array of position arrays
type_array = np.zeros((dumps), dtype=np.ndarray)
box_data = np.zeros((1), dtype=np.ndarray)              # box dimensions
myfile = "pa150_pb150_xa20.gsd"

with hoomd.open(name=myfile, mode='rb') as t:           # open for reading
    snap = t[0]                                         # snap 0th snapshot
    box_data = snap.configuration.box                   # get box dimensions
    for i in range(0,dumps):
        snap = t[i]                                     # take snap of each dump
        type_array[i] = snap.particles.typeid
        position_array[i] = snap.particles.position     # store all particle positions

# try to do a temporary storage that I can then transfer over
pos_A = np.zeros((dumps), dtype=np.ndarray)
pos_B = np.zeros((dumps), dtype=np.ndarray)
tmpA = np.zeros((part_a, 3), dtype=np.float32)
tmpB = np.zeros((part_b, 3), dtype=np.float32)

print(tmpA.shape)
print(tmpB.shape)

for f in range(0, dumps):
    countA = 0
    countB = 0
    for g in range(0, part_num):
        if type_array[f][g] == 0:
            tmpA[countA][0] = position_array[f][g][0]
            tmpA[countA][1] = position_array[f][g][1]
            tmpA[countA][2] = position_array[f][g][2]
            countA += 1
        else:
            tmpB[countB][0] = position_array[f][g][0]
            tmpB[countB][1] = position_array[f][g][1]
            tmpB[countB][2] = position_array[f][g][2]
            countB += 1

    print(tmpA)
    pos_A[f] = tmpA
    pos_B[f] = tmpB


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

# analyze all particles
for j in range(0, dumps):
    
    l_pos = position_array[j]
    my_clusters.computeClusters(l_pos)
    number_clusters[j] = my_clusters.getNumClusters()   # find number of clusters
    ids = my_clusters.getClusterIdx()                   # get cluster ids
    cluster_props.computeProperties(l_pos, ids)
    size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each
    
    for k in range(0, len(size_clusters[j])):
        # the size minimum is a very important value to consider
        if size_clusters[j][k] > 25 and size_clusters[j][k] < part_num:
            tot_size[j] += size_clusters[j][k]
            tot_num[j] += 1

    if tot_num[j] > 0:
        MCS[j] = float(tot_size[j]/tot_num[j])/float(part_num)
        GF[j] = float(part_num - tot_size[j]) / float(part_num)
    
    else:
        MCS[j] = 0
        GF[j] = 1

def getDensityPlease(n):                                # call this function as needed
    l_pos = position_array[n]                           # get ith position array
    my_density.compute(f_box,
                       l_pos,
                       l_pos)
    return my_density.getDensity()

avg_sys_density = np.zeros((1), dtype=np.ndarray)

take_last = dumps - 500
last = dumps - 1
for j in range(take_last, dumps):
    avg_sys_density[0] += getDensityPlease(j)

avg_sys_density[0] /= (dumps - take_last)

# perform the same analysis on species A

tot_size_A = np.zeros((dumps), dtype=np.ndarray)          # number of particles in clusters
tot_num_A = np.zeros((dumps), dtype=np.ndarray)           # total number of clusters
MCS_A = np.zeros((dumps), dtype=np.ndarray)               # Mean cluster size
GF_A = np.zeros((dumps), dtype=np.ndarray)                # Gas fraction

for j in range(0, dumps):
    
    l_pos = pos_A[j]
    my_clusters.computeClusters(l_pos)
    number_clusters[j] = my_clusters.getNumClusters()   # find number of clusters
    ids = my_clusters.getClusterIdx()                   # get cluster ids
    cluster_props.computeProperties(l_pos, ids)
    size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each
    
    for k in range(0, len(size_clusters[j])):
        # the size minimum is a very important value to consider
        if size_clusters[j][k] > 25 and size_clusters[j][k] < part_num:
            tot_size_A[j] += size_clusters[j][k]
            tot_num_A[j] += 1

    if tot_num_A[j] > 0:
        MCS_A[j] = float(tot_size_A[j]/tot_num_A[j])/float(part_a)
        GF_A[j] = float(part_a - tot_size_A[j]) / float(part_a)
    
    else:
        MCS_A[j] = 0
        GF_A[j] = 1

def getDensityA(n):                                     # call this function as needed
    l_pos = pos_A[n]                                    # get ith position array
    my_density.compute(f_box,
                       l_pos,
                       l_pos)
    return my_density.getDensity()

avg_dense_A = np.zeros((1), dtype=np.ndarray)

for j in range(take_last, dumps):
    avg_dense_A[0] += getDensityA(j)

avg_dense_A[0] /= (dumps - take_last)

# perform the same analysis on species B

tot_size_B = np.zeros((dumps), dtype=np.ndarray)          # number of particles in clusters
tot_num_B = np.zeros((dumps), dtype=np.ndarray)           # total number of clusters
MCS_B = np.zeros((dumps), dtype=np.ndarray)               # Mean cluster size
GF_B = np.zeros((dumps), dtype=np.ndarray)                # Gas fraction

for j in range(0, dumps):
    
    l_pos = pos_B[j]
    my_clusters.computeClusters(l_pos)
    number_clusters[j] = my_clusters.getNumClusters()   # find number of clusters
    ids = my_clusters.getClusterIdx()                   # get cluster ids
    cluster_props.computeProperties(l_pos, ids)
    size_clusters[j] = cluster_props.getClusterSizes()  # get number of particles in each
    
    for k in range(0, len(size_clusters[j])):
        # the size minimum is a very important value to consider
        if size_clusters[j][k] > 25 and size_clusters[j][k] < part_num:
            tot_size_B[j] += size_clusters[j][k]
            tot_num_B[j] += 1

    if tot_num_B[j] > 0:
        MCS_B[j] = float(tot_size_B[j]/tot_num_B[j])/float(part_b)
        GF_B[j] = float(part_b - tot_size_B[j]) / float(part_b)
    
    else:
        MCS_B[j] = 0
        GF_B[j] = 1

def getDensityB(n):                                     # call this function as needed
    l_pos = pos_B[n]                                    # get ith position array
    my_density.compute(f_box,
                       l_pos,
                       l_pos)
    return my_density.getDensity()

avg_dense_B = np.zeros((1), dtype=np.ndarray)

for j in range(take_last, dumps):
    avg_dense_B[0] += getDensityB(j)

avg_dense_B[0] /= (dumps - take_last)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes=True)

plt_name  = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a)
plt_name1 = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + "A"
plt_name2 = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + "B"

sns.kdeplot(avg_sys_density[0], shade = True)
sns.kdeplot(avg_dense_A[0], shade = True)
sns.kdeplot(avg_dense_B[0], shade = True)
plt.savefig('avg_density_' + plt_name + '.png', dpi=1000)
plt.close()

sns.kdeplot(getDensityPlease(last), shade = True)
sns.kdeplot(getDensityA(last), shade = True)
sns.kdeplot(getDensityB(last), shade = True)
plt.savefig('final_density_' + plt_name + '.png', dpi=1000)
plt.close()

plt.plot(MCS)
plt.plot(MCS_A)
plt.plot(MCS_B)
plt.savefig('MCS_'+ plt_name + '.png', dpi=1000)
plt.close()

plt.plot(GF)
plt.plot(GF_A)
plt.plot(GF_B)
plt.savefig('GF_'+plt_name+'.png', dpi=1000)
plt.close()
