''' This file will read in the Peclet number and give the predicted swim diffusivity '''

import sys
import os
import numpy as np
import math
import matplotlib.pyplot as plt

#pe = float(sys.argv[1])

def swimDiffusivity(peclet):
    '''In our implementation Pe = v_p and kT = 1.0'''
    '''D_{Swim} = v_{p}^{2} * tau_{R} / 6 '''
    D_s = (peclet ** 2) * (1.0 / 3.0) / 6.0
    return D_s

xs = np.arange(0, 500.0, 0.001)
ys = np.zeros_like(xs)
for i in xrange(len(xs)):
    ys[i] = swimDiffusivity(xs[i])

plt.plot(xs, ys, lw=2)
plt.xlabel(r'Peclet Number')
plt.ylabel(r'Swim Diffusivity $(D_{Swim})$')
plt.show()
