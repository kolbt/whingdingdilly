''' This file will read in the Peclet number and give the predicted swim diffusivity '''

import sys
import os
import numpy as np
import math

pe = float(sys.argv[1])

def swimDiffusivity(peclet):
    '''In our implementation Pe = v_p and kT = 1.0'''
    D_s = (peclet ** 2 * 1.0) / 6.0
    return D_s
