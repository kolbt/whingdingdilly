'''
This file will:
    1. Read in:
        msd
        q_clust
        type_array
    2. Sort msd data into appropriate arrays
    3. Write out MSD for:
        A   : liquid, gas
        B   : liquid, gas
        All : liquid, gas
'''
import numpy as np

def writeArrays(msd_val, q_clust, ids, type):
    
    part_num = len(type)
    lq_all = 0
    lq_a = 0
    lq_b = 0
    gs_all = 0
    gs_a = 0
    gs_b = 0
    lq_a_count = 0
    lq_b_count = 0
    gs_a_count = 0
    gs_b_count = 0
    num = 0
    den = 0
    perc_A = 0
    
    for b in range(0, part_num):
        
        if q_clust[ids[b]] == 1:            # check if in liquid
            lq_all += msd_val[b]            # add to tot. lq. msd
            if type[b] == 0:                # type A case
                lq_a += msd_val[b]
                lq_a_count += 1
            else:
                lq_b += msd_val[b]
                lq_b_count += 1
        else:                               # else, particle is gas
            gs_all += msd_val[b]            # add to tot. gs. msd
            if type[b] == 0:                # type A case
                gs_a += msd_val[b]
                gs_a_count += 1
            else:
                gs_b += msd_val[b]
                gs_b_count += 1

    if lq_a_count != 0: lq_a /= lq_a_count
    if lq_b_count != 0: lq_b /= lq_b_count
    if gs_a_count != 0: gs_a /= gs_a_count
    if gs_b_count != 0: gs_a /= gs_b_count
    if ( lq_a_count + lq_b_count ) != 0: lq_all /= (lq_a_count + lq_b_count)
    if ( gs_a_count + gs_b_count ) != 0: gs_all /= (gs_a_count + gs_b_count)

    num = lq_a_count
    den = ( lq_a_count + lq_b_count )

    if den != 0:
        perc_A = float(num) / float(den)

    return lq_a, lq_b, gs_a, gs_b, gs_all, lq_all, perc_A
