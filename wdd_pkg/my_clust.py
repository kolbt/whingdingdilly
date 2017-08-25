'''
This module will look at the clusters in our system, and get some
information out from that.
    1. Read in
        size_clusters[k]
        size_min (optional)
    2. Compute and write
        largest cluster size
        mean cluster size
        gas fraction
'''

import numpy as np

def clustDat(size, part_num, step, f_largest, size_min=1000, write_flag='False'):
    
    tot_clust = 0   # number of particles in dense phase
    num_clust = 0   # number of unique clusters
    l_clust = 0     # size of largest cluster
    
    for c in range(0, len(size)):
        if size[c] > size_min:      # if cluster is large enough, accept
            tot_clust += size[c]    # contributes to dense phase
            num_clust += 1          # contributes to number of clusters
            if size[c] > l_clust:   # is cluster new largest cluster
                l_clust = size[c]   # write to largest cluster var

    # compute mean cluster size and gas fraction
    if num_clust > 0:
        mcs = float(tot_clust/num_clust) / float(part_num)
        gf  = float(part_num-tot_clust) / float(part_num)
    else:
        mcs = 0
        gf  = 1

    # write largest cluster size to text file
    if write_flag == True:
        if step == 0:
            a_w = 'w'
        else:
            a_w = 'a'
        f = open(f_largest, a_w)
        f.write(str(l_clust) + '\n')
        f.close()

    return l_clust, mcs, gf
