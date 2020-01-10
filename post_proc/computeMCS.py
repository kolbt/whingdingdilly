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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def computeTauPerTstep(epsilon, mindt=0.000001):
    '''Read in epsilon, output tauBrownian per timestep'''
    if epsilon != 1.:
        mindt=0.00001
    kBT = 1.0
    tstepPerTau = float(epsilon / (kBT * mindt))
    return 1. / tstepPerTau
    
# Get infile and open
inFile = str(sys.argv[1])
f = hoomd.open(name=inFile, mode='rb')
# Inside and outside activity from command line
peA = int(sys.argv[2])
peB = int(sys.argv[3])
parFrac = int(sys.argv[4])
eps = float(sys.argv[5])
try:
    phi = float(sys.argv[6])
    intPhi = int(phi)
    phi /= 100.
except:
    phi = 0.6
    intPhi = 60

# Outfile to write data to
txtFile = 'MCS_pa' + str(peA) +\
          '_pb' + str(peB) +\
          '_xa' + str(parFrac) +\
          '_phi' + str(intPhi) +\
          '.txt'
# Write file headings
g = open(txtFile, 'w')
g.write('Timestep'.center(10) + ' ' +\
        'N_clusts'.center(10) + ' ' +\
        'all_parts'.center(10) + ' ' +\
        'thr_clusts'.center(10) + ' ' +\
        'thr_parts'.center(10) + ' ' +\
        'thr_c1000'.center(10) + ' ' +\
        'thr_p1000'.center(10) + '\n')
g.close()

start = 0                   # first frame to process
dumps = int(f.__len__())    # get number of timesteps dumped
end = dumps                 # final frame to process

box_data = np.zeros((1), dtype=np.ndarray)  # box dimension holder
r_cut = 2**(1./6.)                          # potential cutoff
tauPerDT = computeTauPerTstep(epsilon=eps)  # brownian time per timestep
min_size = 20
thresh_large = 1000

with hoomd.open(name=inFile, mode='rb') as t:
    snap = t[0]
    first_tstep = snap.configuration.step
    box_data = snap.configuration.box
    l_box = box_data[0]
    typ = snap.particles.typeid
    partNum = len(typ)
    # Set up cluster computation using box
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
        tst *= tauPerDT                             # convert to Brownian time
        
        # Compute clusters for this timestep
        system = freud.AABBQuery(f_box, pos)
        my_clust.compute(system, neighbors={'r_max': 1.005})
        ids = my_clust.cluster_idx              # get id of each cluster
        c_props.compute(system, ids)            # find cluster properties
        clust_size = c_props.sizes              # find cluster sizes
        
        nClusts = len(clust_size)
        nClustThr = 0
        tot_clust = 0
        nC1000 = 0
        nP1000 = 0
        for k in range(0, nClusts):
            if clust_size[k] > min_size:
                nClustThr += 1
                tot_clust += clust_size[k]
                if clust_size[k] > thresh_large:
                    nC1000 += 1
                    nP1000 += clust_size[k]
        
        g = open(txtFile, 'a')
        g.write('{0:.1f}'.format(tst).center(10) + ' ' +\
                str(nClusts).center(10) + ' ' +\
                str(partNum).center(10) + ' ' +\
                str(nClustThr).center(10) + ' ' +\
                str(tot_clust).center(10) + ' ' +\
                str(nC1000).center(10) + ' ' +\
                str(nP1000).center(10) + '\n')
        g.close()
        
        
        
    
