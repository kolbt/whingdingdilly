'''
#                           This is an 80 character line                       #
PURPOSE:  The intent of this file is to get out a few basic values as text files

                    File   :   Column 1    Column 2    Column 3    Column 4
   gas_pea#_peb#_xa#.txt   :   time        gas_A       gas_B       gas_total
 dense_pea#_peb#_xa#.txt   :   time        dense_A     dense_B     dense_total
lclust_pea#_peb#_xa#.txt   :   time        l_clust
         
Designations Explained:
         
            gas_A : number of A-type particles in gas phase
            gas_B : number of B-type particles in gas phase
        gas_total : total number of gaseous partilces
          dense_A : number of A-type particles in dense phase
          dense_B : number of B-type particles in dense phase
      dense_total : total number of dense phase particles
          l_clust : largest individual cluster
'''

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

import numpy as np

in_file = "pa"+str(pe_a)+\
"_pb"+str(pe_b)+\
"_xa"+str(part_perc_a)+\
".gsd"

f = hoomd.open(name=in_file, mode='rb') # open gsd file with hoomd
dumps = f.__len__()                     # get number of timesteps dumped

start = 0       # gives first frame to read
end = dumps     # gives last frame to read

positions = np.zeros((end), dtype=np.ndarray)           # array of positions
types = np.zeros((end), dtype=np.ndarray)               # particle types
box_data = np.zeros((1), dtype=np.ndarray)              # box dimensions
timesteps = np.zeros((end), dtype=np.float64)           # timesteps


with hoomd.open(name=in_file, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    for iii in range(start, end):
        snap = t[iii]                               # snapshot of frame
        types[iii] = snap.particles.typeid          # get types
        positions[iii] = snap.particles.positions   # get positions
        timesteps[iii] = snap.configuration.step    # get timestep
