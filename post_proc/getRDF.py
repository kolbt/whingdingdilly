'''
#                           This is an 80 character line                       #
Output RDF of the liquid phase
'''

import sys

# Run locally
sys.path.append('/Users/kolbt/Desktop/compiled/hoomd-blue/build')
sys.path.append('/Users/kolbt/Desktop/compiled/gsd/build')
# Run on the cpu
sys.path.append('/nas/longleaf/home/kolbt/programs/cpu-hoomd/hoomd-blue/build')
# Run on the gpu

sys.path.append('/nas/longleaf/home/kolbt/programs/hoomd_2.2.1/hoomd-blue/build')
sys.path.append('/nas/longleaf/home/kolbt/programs/gsd/build')

import gsd
from gsd import hoomd
from gsd import pygsd

import freud
from freud import parallel
from freud import box
from freud import density
from freud import cluster

import numpy as np
import math
import random
from scipy import stats

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.collections
from matplotlib import colors
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.ticker as ticker

# Here are my rc parameters for matplotlib
mpl.rc('font', serif='Helvetica Neue')
mpl.rcParams.update({'font.size': 9})
#mpl.rcParams['figure.figsize'] = 3.2, 2.8
#mpl.rcParams['figure.dpi'] = 2000
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['axes.linewidth'] = 1.5
    
def distComps(point1, point2x, point2y):
    '''Given points output x, y and r'''
    dx = point2x - point1[0]
    dy = point2y - point1[1]
    r = np.sqrt((dx**2) + (dy**2))
    return dx, dy, r
    
def computeFLJ(r, dx, dy, eps):
    sig = 1.
    f = (24. * eps / r) * ( (2*((sig/r)**12)) - ((sig/r)**6) )
    fx = f * (dx) / r
    fy = f * (dy) / r
    return fx, fy

def computeTauPerTstep(epsilon, mindt=0.000001):
    '''Read in epsilon, output tauBrownian per timestep'''
#    if epsilon != 1.:
#        mindt=0.00001
    kBT = 1.0
    tstepPerTau = float(epsilon / (kBT * mindt))
    return 1. / tstepPerTau

def roundUp(n, decimals=0):
    '''Round up size of bins to account for floating point inaccuracy'''
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier
    
def getNBins(length, minSz=(2**(1./6.))):
    "Given box size, return number of bins"
    initGuess = int(length) + 1
    NBins = initGuess
    # This loop only exits on function return
    while True:
        if length / NBins > minSz:
            return NBins
        else:
            NBins -= 1
            
# Get lattice spacing for particle size
def ljForce(r, eps, sigma=1.):
    div = (sigma/r)
    dU = (24. * eps / r) * ((2*(div**12)) - (div)**6)
    return dU
    
def avgCollisionForce(pe, power=1.):
    '''Computed from the integral of possible angles'''
    magnitude = np.sqrt(28.)
    return (magnitude * (pe)) / (np.pi)
    
def conForRClust(pe, eps):
    out = []
    r = 1.112
    skip = [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001]
    for j in skip:
        while ljForce(r, eps) < avgCollisionForce(pe):
            r -= j
        r += j
    out = r
    return out

# Get infile and open
inFile = str(sys.argv[1])
if inFile[0:7] == "cluster":
    add = 'cluster_'
else:
    add = ''
outF = inFile[:-4]
    
f = hoomd.open(name=inFile, mode='rb')
# Inside and outside activity from command line
peA = float(sys.argv[2])
peB = float(sys.argv[3])
parFrac = float(sys.argv[4])
eps = float(sys.argv[5])
try:
    phi = float(sys.argv[6])
    intPhi = int(phi)
    phi /= 100.
except:
    phi = 0.6
    intPhi = 60
try:
    dtau = float(sys.argv[7])
except:
    dtau = 0.000001
    
lat = conForRClust(peA, eps)
    
# Get number of particles for textfile name
with hoomd.open(name=inFile, mode='rb') as t:
    snap = t[0]
    typ = snap.particles.typeid
    partNum = len(typ)

# Outfile to write data to
base = add + 'hmap_pressure_pa' + str(peA) +\
       '_pb' + str(peB) +\
       '_xa' + str(parFrac) +\
       '_phi' + str(intPhi) +\
       '_ep' + '{0:.3f}'.format(eps) +\
       '_N' + str(partNum)

start = 0                   # first frame to process
dumps = int(f.__len__())    # get number of timesteps dumped
end = dumps                 # final frame to process
start = end - 1

box_data = np.zeros((1), dtype=np.ndarray)  # box dimension holder
r_cut = 2**(1./6.)                          # potential cutoff
tauPerDT = computeTauPerTstep(epsilon=eps)  # brownian time per timestep

# Make a text file for this (so I can overlay everything)
outTxt = 'rdf_' + outF + '.txt'
g = open(outTxt, 'w') # write file headings
g.write('tauB'.center(30) + ' ' +\
        'r'.center(30) + ' ' +\
        'g(r)'.center(30) + '\n')
g.close()

with hoomd.open(name=inFile, mode='rb') as t:
    snap = t[0]
    first_tstep = snap.configuration.step
    box_data = snap.configuration.box
    l_box = box_data[0]
    h_box = l_box / 2.
    typ = snap.particles.typeid
    partNum = len(typ)
    
    # Width, in distance units, of bin
    wBins = 0.2
    # Number of bins given this distance
    nBins = l_box / wBins
    # Distance to compute RDF for
    rstop = 10.
    # x-axis for plots of rdf
    r = np.arange(0., rstop, wBins)
    
    # Instantiate cluster computation
    f_box = box.Box(Lx=l_box, Ly=l_box, is2D=True)
    my_clust = cluster.Cluster()
    c_props = cluster.ClusterProperties()
    
    # Loop through each timestep
    for j in range(start, end):
        snap = t[j]
        # Easier accessors
        pos = snap.particles.position               # position
        pos[:,-1] = 0.0
        xy = np.delete(pos, 2, 1)
        typ = snap.particles.typeid                 # type
        tst = snap.configuration.step               # timestep
        tst -= first_tstep                          # normalize by first timestep
        tst *= dtau                                 # convert to Brownian time
        
        # Compute clusters for this timestep
        system = freud.AABBQuery(f_box, f_box.wrap(pos))
#        # Compute neighbor list for only largest cluster
#        my_clust.compute(system, neighbors={'r_max': 1.0})
#        ids = my_clust.cluster_idx              # get id of each cluster
#        c_props.compute(system, ids)            # find cluster properties
#        clust_size = c_props.sizes              # find cluster sizes
#        lcID = np.where(clust_size == np.amax(clust_size))
        
        # Get positions of particles in the largest cluster
        rdf = freud.density.RDF(bins=nBins, r_max=rstop, r_min=0.1, normalize=True)
        rdf.compute(system=system, query_points=pos)
        
        # Get frame for output
        pad = str(j).zfill(4)
        
        # Steal the data from the INANE in-module plotting
        fig = plt.figure()
        ax = plt.gca()
        myrdf = rdf.plot(ax)
        line = ax.lines[0]
        xs = line.get_xdata()
        ys = line.get_ydata()
        plt.close()
        
        # Plot it how I want to
        fig = plt.figure()
        ax = plt.gca()
        ax.plot(xs/lat, ys, lw=1.5)
        ax.set_xlim(0.1, 10.)
        ax.set_ylim(0,)
        # Set the x and y minor ticks
        loc = ticker.MultipleLocator(base=1.)
        ax.xaxis.set_major_locator(loc)
        loc = ticker.MultipleLocator(base=0.2)
        ax.xaxis.set_minor_locator(loc)
        # Tick width and height
        ax.xaxis.set_tick_params(width=1.5, size=5.)
        ax.yaxis.set_tick_params(width=1.5, size=5.)
        ax.xaxis.set_tick_params(which='minor', width=1.25, size=3.)
        ax.yaxis.set_tick_params(which='minor', width=1.25, size=3.)
        ax.set_xlabel(r'$r/a$')
        ax.set_ylabel(r'$g(r)$')
        plt.tight_layout()
        plt.savefig('rdf_' + outF + '_' + pad + '.pdf', dpi=100)
        plt.close()
        
        # Append data to file
        g = open(outTxt, 'a')
        for k in range(0, len(xs)):
            g.write('{0:.3f}'.format(tst).center(30) + ' ')
            g.write('{0:.3f}'.format(xs[k]).center(30) + ' ')
            g.write('{0:.3f}'.format(ys[k]).center(30) + '\n')
        g.close()
