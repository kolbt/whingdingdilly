'''
This file is intended to sort positions by type.  So,
we'll end up with pos_A and pos_B arrays of position
data.
    1. Read in
        position_array[step]
        type_array[step]
    2. Sort
    3. Write out
        pos_A
        pos_B
'''

import numpy as np

def posByType(pos, type, A):

    part_num = len(type)
    countA = 0
    countB = 0
    pos_A = np.zeros((A, 3), dtype=np.float32)
    pos_B = np.zeros((part_num - A, 3), dtype=np.float32)

    for d in range(0, part_num):
        if type[d] == 0:
            pos_A[countA][0] = pos[d][0]
            pos_A[countA][1] = pos[d][1]
            pos_A[countA][2] = pos[d][2]
            countA += 1
        else:
            pos_B[countB][0] = pos[d][0]
            pos_B[countB][1] = pos[d][1]
            pos_B[countB][2] = pos[d][2]
            countB += 1

    return pos_A, pos_B
