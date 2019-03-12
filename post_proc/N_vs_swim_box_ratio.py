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

def lpBoxRatio(lp, box):
    return float(lp) / float(box)

pes = np.arange(0.0, 500.0, 0.01)
ys1k = np.zeros_like(pes)
ys10k = np.zeros_like(pes)
ys100k = np.zeros_like(pes)
ys1m = np.zeros_like(pes)

l_box1k = nToLBox(1000)
l_box10k = nToLBox(10000)
l_box100k = nToLBox(100000)
l_box1m = nToLBox(1000000)
for i in xrange(len(pes)):
    ys1k[i] = lpBoxRatio(lp(pes[i]), l_box1k)
    ys10k[i] = lpBoxRatio(lp(pes[i]), l_box10k)
    ys100k[i] = lpBoxRatio(lp(pes[i]), l_box100k)
    ys1m[i] = lpBoxRatio(lp(pes[i]), l_box1m)

# We want to plot this relationship at constant N values
plt.plot(pes, ys1k, label=r'$N=1000$')
plt.plot(pes, ys10k, label=r'$N=10000$')
plt.plot(pes, ys100k, label=r'$N=100000$')
plt.plot(pes, ys1m, label=r'$N=1000000$')
plt.xlabel(r'Peclet Number $(Pe)$')
plt.ylabel(r'Ratio of $l_{p}:l_{box}$')
plt.xlim(0.0, 500.0)
plt.ylim(0.0, 0.5)
plt.legend()
plt.savefig('N_vs_lp_to_box_ratio.png', dpi=1000)
plt.close()
