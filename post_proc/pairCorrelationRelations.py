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
".gsd"

f = hoomd.open(name=in_file, mode='rb') # open gsd file with hoomd
dumps = f.__len__()                     # get number of timesteps dumped

start = 0       # gives first frame to read
end = dumps     # gives last frame to read
start = 2
end = 3

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

    radialDF.compute(f_box, pos, pos)
    myRDF = radialDF.getRDF()
    sFact = fft(myRDF)

    fig = plt.figure()
    ax = fig.add_subplot(231)
    plt.plot(r, myRDF)
    plt.xlim(0.0, searchRange)
    plt.ylim(0.0)
    plt.xlabel(r'r $(\sigma)$')
    plt.ylabel(r'g(r)')
    plt.title("Pair correlation function")

    ax = fig.add_subplot(232)
    plt.plot(k, 2.0/N * np.abs(sFact[0:N//2]))
    plt.xlabel(r'$q\sigma$')
    plt.ylabel(r'S(q)')
    plt.title("SciPy linspace Method")

    ax = fig.add_subplot(233)
    plt.loglog(1/r, np.abs(sFact))
    plt.xlabel(r'$q\sigma^{-1}$')
    plt.ylabel(r'S(q)')
    plt.title("Wavevector = inverse distance (loglog)")

    ax = fig.add_subplot(234)
    plt.plot(1/r, np.abs(sFact))
    plt.xlabel(r'$q\sigma^{-1}$')
    plt.ylabel(r'S(q)')
    plt.title("Wavevector = inverse distance")

    ax = fig.add_subplot(235)
    plt.plot(r, np.abs(sFact))
    plt.xlabel(r'$q\sigma$')
    plt.ylabel(r'S(q)')
    plt.title("Wavevector = distance")

    ax = fig.add_subplot(236)
    plt.loglog(r, np.abs(sFact))
    plt.xlabel(r'$q\sigma$')
    plt.ylabel(r'S(q)')
#    plt.ylim(0, 200)
    plt.title("Wavevector = distance (loglog)")
    plt.show()

    plt.plot(1/r, np.abs(sFact))
    plt.xlabel(r'$q(\sigma^{-1})$')
    plt.ylabel(r'S(q)')
#    plt.xlim(0, 10)
#    plt.ylim(0, 200)
    plt.title("Wavevector = inverse distance")
    plt.show()


#ffmpeg -framerate 10 -i nBins100_pa150_pb500_xa50_step_%d.png\
# -vcodec libx264 -s 1000x1000 -pix_fmt yuv420p -threads 1\
# nBins100_pa150_pb500_xa50.mp4

