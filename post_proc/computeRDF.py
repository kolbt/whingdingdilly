'''
#                           This is an 80 character line                       #
First: obtain the pair correlation function (using Freud)
Then: take fourier transform of it (obtain structure factor)
And then: take inverse first moment (obtain coarsening length)
Finally: plot all that shit
'''

# Imports and loading the .gsd file
import sys

pe_a = int(sys.argv[1])                     # activity A
pe_b = int(sys.argv[2])                     # activity B
part_perc_a = int(sys.argv[3])              # percentage A particles
part_frac_a = float(part_perc_a) / 100.0    # fraction A particles
hoomd_path = str(sys.argv[4])               # local path to hoomd-blue
gsd_path = str(sys.argv[5])                 # local path to gsd
eps = int(sys.argv[6])

sys.path.append(hoomd_path)     # ensure hoomd is in your python path
sys.path.append(gsd_path)       # ensure gsd is in your python path

import hoomd
from hoomd import md
from hoomd import deprecated

import gsd
from gsd import hoomd
from gsd import pygsd

import freud
from freud import parallel
from freud import box
from freud import density
from freud import cluster

import numpy as np
from scipy.fftpack import fft, ifft

import matplotlib.pyplot as plt
from matplotlib import colors

import math

# File to read from
in_file = "pa"+str(pe_a)+\
"_pb"+str(pe_b)+\
"_xa"+str(part_perc_a)+\
"_ep"+str(eps)+\
".gsd"

out = "pa"+str(pe_a)+\
"_pb"+str(pe_b)+\
"_xa"+str(part_perc_a)+\
"_ep"+str(eps)

f = hoomd.open(name=in_file, mode='rb') # open gsd file with hoomd
dumps = f.__len__()                     # get number of timesteps dumped

start = 0       # gives first frame to read
end = dumps     # gives last frame to read

positions = np.zeros((end), dtype=np.ndarray)       # array of positions
types = np.zeros((end), dtype=np.ndarray)           # particle types
box_data = np.zeros((1), dtype=np.ndarray)          # box dimensions
timesteps = np.zeros((end), dtype=np.float64)       # timesteps
orient = np.zeros((end), dtype=np.ndarray)          # orientations

# Get relevant data from .gsd file
with hoomd.open(name=in_file, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    for iii in range(start, end):
        snap = t[iii]                               # snapshot of frame
        types[iii] = snap.particles.typeid          # get types
        positions[iii] = snap.particles.position    # get positions
        orient[iii] = snap.particles.orientation    # get orientation
        timesteps[iii] = snap.configuration.step    # get timestep

timesteps -= timesteps[0]       # get rid of brownian run time

# Get number of each type of particle
part_num = len(types[start])
part_A = int(part_num * part_frac_a)
part_B = part_num - part_A

# Feed data into freud analysis software
l_box = box_data[0]
h_box = l_box / 2.0
a_box = l_box * l_box
f_box = box.Box(Lx = l_box, Ly = l_box, is2D = True)    # make freud box

nBins = 1000
widthBin = 0.005
searchRange = nBins * widthBin
radialDF = freud.density.RDF(searchRange, widthBin)

#r = np.arange(0.0, searchRange, widthBin)
#k = np.arange(0.0, )

N = nBins                       # number of samples
T = widthBin                    # spacing between samples
r = np.linspace(0.0, N*T, N)    # 0 through searchRange with Nbins
# // is floor division, adjusts to left in number line
k = np.linspace(0.0, 1.0/(2.0*T), N//2)

for iii in range(start, end):
    
    # Easier accessors
    pos = positions[iii]
    typ = types[iii]
    dir = orient[iii]

    # Compute RDF for all particles
    radialDF.compute(f_box, pos, pos)
    myRDF = radialDF.getRDF()

    plt.plot(r, myRDF, label='All')

    posA = np.zeros((part_A, 3), dtype=np.float64)
    posB = np.zeros((part_B, 3), dtype=np.float64)
    countA = 0
    countB = 0
    for j in xrange(len(pos)):
        if typ[j] == 0:
            posA[countA][0] = pos[j][0]
            posA[countA][1] = pos[j][1]
            countA += 1
        else:
            posB[countB][0] = pos[j][0]
            posB[countB][1] = pos[j][1]
            countB += 1

    # Compute RDF for AA
    # radialDF.reset()
    radialDF.compute(f_box, posA, posA)
    rdfAA = radialDF.getRDF()

    plt.plot(r, rdfAA, label='AA')

    # Compute RDF for BB
    # radialDF.reset()
    radialDF.compute(f_box, posB, posB)
    rdfBB = radialDF.getRDF()

    plt.plot(r, rdfBB, label='BB')

    # Compute RDF for AB
    # radialDF.reset()
    radialDF.compute(f_box, posA, posB)
    rdfAB = radialDF.getRDF()

    plt.plot(r, rdfAB, label='AB')

    # Compute RDF for BA
    # radialDF.reset()
    radialDF.compute(f_box, posA, posB)
    rdfBA = radialDF.getRDF()

    plt.plot(r, rdfBA, label='BA')

    plt.xlim(0.0, searchRange)
    # plt.ylim(0.0, 10.0)
    plt.xlabel(r'r $(\sigma)$')
    plt.ylabel(r'g(r)')
    plt.legend()
    plt.savefig('RDF_' + out + '_fm' + str(iii) + '.png', dpi=1000)
    plt.close()
