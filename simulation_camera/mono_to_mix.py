'''
#                           This is an 80 character line                       #
1.) Read existing monodisperse .gsd file
2.) Swap the type of first N * particle fraction particles
3.) Save single (final) single frame of .gsd
'''

import sys

inFile = str(sys.argv[1])                   # filename
pe = int(sys.argv[2])                       # activity of mono system
part_perc_a = int(sys.argv[3])              # percentage A particles
part_frac_a = float(part_perc_a) / 100.0    # fraction A particles
hoomd_path = str(sys.argv[4])               # local path to hoomd-blue
gsd_path = str(sys.argv[5])                 # local path to gsd

sys.path.append(hoomd_path)     # ensure hoomd is in your python path
sys.path.append(gsd_path)       # ensure gsd is in your python path

# This gets the length of the hoomd file
import gsd
from gsd import hoomd
from gsd import pygsd

f = gsd.hoomd.open(name=inFile, mode='rb')  # open gsd file with hoomd
dumps = int(f.__len__()) - 1                # get number of timesteps dumped

import hoomd
from hoomd import md
from hoomd import deprecated
import numpy as np

hoomd.context.initialize()
system = hoomd.init.read_gsd(filename=inFile,
                             restart=None,
                             frame=dumps,
                             time_step=None)

# Add B-type particles
system.particles.types.add('B')
snapshot = system.take_snapshot()
part_num = len(system.particles)
part_a = part_num * part_frac_a     # get the total number of A particles
part_a = int(part_a)                # make sure it is an integer
part_b = part_num - part_a          # get the total number of B particles
part_b = int(part_b)                # make sure it is an integer
mid = int(part_a)                   # starting index to assign B particles

# Assign particles a type within the snapshot
if part_perc_a == 0:                    # take care of all b case
    mid = 0
    for i in range(mid, part_num):
        system.particles[i].type = 'B'
elif part_perc_a != 100:                  # mix of each or all A
    for i in range(mid, part_num):
        system.particles[i].type = 'B'

# Assigning groups and lengths to particles
all = hoomd.group.all()
gA = hoomd.group.type(type = 'A', update=True)
gB = hoomd.group.type(type = 'B', update=True)
N = len(all)
Na = len(gA)
Nb = len(gB)

# Then use dump to write to different filename
outFile =  'final_mono_pa' + str(pe) + '_xa' + str(part_perc_a) + '.gsd'
hoomd.dump.gsd(outFile,
               period=None,
               group=all,
               overwrite=True,
               dynamic=['attribute', 'property', 'momentum'])
