'''
Things this file should do
   1. sort the list of cluster ids
   2. initialize a 'q_clust' array
   3. fill array by testing if cluster has multiple particles
'''
import numpy as np

# Read in the clust ID array
def binaryCluster( ids, clust_num, size_min=1000 ):
    
    part_num = len(ids)
    sort = np.sort(ids)
    biclust = np.zeros((clust_num), dtype=np.bool_)
    index = 0
    for a in range(0, clust_num):
        add_clust = 0
        while 1:
            add_clust += 1
            if index == part_num:       # break if index is too large
                break
            if sort[index] != a:        # break if ID changes
                break
            if add_clust == 1:          # all particles appear once
                biclust[a] = 0
            if add_clust > size_min:    # only multiple ids appear twice
                biclust[a] = 1
            index += 1                  # increment index

    return biclust
