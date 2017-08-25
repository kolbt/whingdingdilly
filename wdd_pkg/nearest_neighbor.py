'''
The purpose of this file is to find the number of
nearest neighbors for each type of particle
correlation: a with a neighbors, ab, ba, and aa.
    1. Read in
        q_clust[j]
        ids_tot[j]
        type_array[j]
        a_liq_count
        b_liq_count
    2. Compute nearest neighbors with kdtrees
    3. Returns
        avg number of nearest neighbors for each
        correlation
        optional: dense phase A and B positions
        can heatmap for indiv. tstep
'''

import numpy as np
import scipy.spatial as spatial

def nearNeigh(q_clust, ids, pos, type):

    part_num = len(ids)
    
    a_num = 0
    b_num = 0
    for f in range(0, part_num):
        if q_clust[ids[f]] == 1:
            if type[f] == 0:
                a_num += 1
            else:
                b_num += 1
    
    # both a and b particles are in dense phase
    if a_num != 0 and b_num != 0:
        a_count = 0
        b_count = 0
        A_dpos = np.zeros((a_num, 2), dtype=np.float32)
        B_dpos = np.zeros((b_num, 2), dtype=np.float32)
        # if in the dense phase, get position by type
        for e in range(0, part_num):
            if q_clust[ids[e]] == 1:
                if type[e] == 0:
                    A_dpos[a_count, 0] = pos[e][0]
                    A_dpos[a_count, 1] = pos[e][1]
                    a_count += 1
                else:
                    B_dpos[b_count, 0] = pos[e][0]
                    B_dpos[b_count, 1] = pos[e][1]
                    b_count += 1
        
        a_tree = spatial.KDTree(A_dpos)
        b_tree = spatial.KDTree(B_dpos)
        radius = 1.0
        
        # How many A neighbors does the avg A particle have
        aa = a_tree.query_ball_tree(a_tree, radius)
        num_aa = np.array(map(len, aa), dtype=np.float32)
        num_aa -= 1.0           # can't reference itself
        if len(num_aa) != 0:
            avg_aa = (np.sum(num_aa)/len(num_aa))

        # How many B neighbors does the avg A particle have
        ab = a_tree.query_ball_tree(b_tree, radius)
        num_ab = np.array(map(len, ab), dtype=np.float32)
        if len(num_ab) != 0:
            avg_ab = (np.sum(num_ab)/len(num_ab))

        # How many A neighbors does the avg B particle have
        ba = b_tree.query_ball_tree(a_tree, radius)
        num_ba = np.array(map(len, ba), dtype=np.float32)
        if len(num_ba) != 0:
            avg_ba = (np.sum(num_ba)/len(num_ba))

        # How many B neighbors does the avg B particle have
        bb = b_tree.query_ball_tree(b_tree, radius)
        num_bb = np.array(map(len, bb), dtype=np.float32)
        num_bb -= 1.0           # can't reference itself
        if len(num_bb) != 0:
            avg_bb = (np.sum(num_bb)/len(num_bb))

        return avg_aa, avg_ab, avg_ba, avg_bb


    elif a_num == 0 and b_num == 0:
        
        return 0, 0, 0, 0


    elif a_num == 0 and b_num != 0:
        b_count = 0
        B_dpos = np.zeros((b_num, 2), dtype=np.float32)
        for e in range(0, part_num):
            if q_clust[ids[e]] == 1:
                B_dpos[b_count, 0] = pos[e][0]
                B_dpos[b_count, 1] = pos[e][1]
                b_count += 1

        b_tree = spatial.KDTree(B_dpos)
        radius = 1.0
        bb = b_tree.query_ball_tree(b_tree, radius)
        num_bb = np.array(map(len, bb), dtype=np.float32)
        num_bb -= 1.0           # can't reference itself
        if len(num_bb) != 0:
            avg_bb = (np.sum(num_bb)/len(num_bb))

        return 0, 0, 0, avg_bb


    elif a_num != 0 and b_num == 0:
        a_count = 0
        A_dpos = np.zeros((a_num, 2), dtype=np.float32)
        for e in range(0, part_num):
            if q_clust[ids[e]] == 1:
                A_dpos[a_count, 0] = pos[e][0]
                A_dpos[a_count, 1] = pos[e][1]
                a_count += 1
    
        a_tree = spatial.KDTree(A_dpos)
        radius = 1.0
        aa = a_tree.query_ball_tree(a_tree, radius)
        num_aa = np.array(map(len, aa), dtype=np.float32)
        num_aa -= 1.0           # can't reference itself
        if len(num_aa) != 0:
            avg_aa = (np.sum(num_aa)/len(num_aa))

        return avg_aa, 0, 0, 0
