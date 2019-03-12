'''Plot relationship of N to the ratio of persistence length : box length'''

import sys
import os
import numpy as np
import math
import matplotlib.pyplot as plt

def areaCircle(r):
    return np.pi * (r ** 2)

sigma = 1.0
def nToLBox(N):
    A_box = (N * areaCircle(sigma / 2.0)) / 0.6
    L_box = np.sqrt(A_box)
    return L_box

def lp(Pe):
    return Pe / 3.0
