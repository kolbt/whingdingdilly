#Arguments in order:
#    1.) hoomd_path (string): the path to hoomd on your machine
#    2.) tsteps (int): the number of timesteps your simulation is running
#    3.) dump_freq (int): the frequency of gsd dumps
#    4.) part_frac_a (float): particle fraction of type a particles
#    5.) pe_a (float): activity of type a particles
#    6.) pe_b (float): activity of tyep b particles

import sys


#hoomd_path = str(sys.argv[1])
hoomd_path = "/Users/kolbt/Desktop/compiled/hoomd-blue/build"
#tsteps = int(sys.argv[2])
tsteps = 10000000
#dump_freq = int(sys.argv[3])
dump_freq = 20000
#part_frac_a = float(sys.argv[4])
part_perc_a = 40
part_frac_a = float(part_perc_a) / 100.0
#pe_a = float(sys.argv[5])
pe_a = 0
#pe_b = float(sys.argv[6])
pe_b = 0

# calculate number of tsteps which are dumped
dumps = tsteps/dump_freq

###
import numpy as np

dump_list = np.zeros((70), dtype=int)
value_to_dump = 0
jumper = 1
count = 1
for jjj in range(1,70):
    if (count-2) % 9 == 0 and count != 2:
        jumper *= 10
    value_to_dump += jumper
    dump_list[jjj] = value_to_dump
    count += 1
dump_list += 110000
###

sys.path.append(hoomd_path)

import hoomd
from hoomd import md
from hoomd import deprecated

#initialize system randomly, can specify GPU execution here

part_num = 15000

hoomd.context.initialize()
system = hoomd.deprecated.init.create_random(N = part_num,
                                             phi_p = 0.6,
                                             name = 'A',
                                             min_dist = 0.5,
                                             seed = 230958,
                                             dimensions = 2)

system.particles.types.add('B')
snapshot = system.take_snapshot()

part_a = part_num * part_frac_a         # get the total number of A particles
part_a = int(part_a)
part_b = part_num - part_a              # get the total number of B particles
part_b = int(part_b)
mid = int(part_a)                       # starting point for assigning B particles

if part_perc_a == 0:                    # take care of all b case
    mid = 0
    for i in range(mid,part_num):
        system.particles[i].type = 'B'
elif part_perc_a != 100:                # mix of each
    for i in range(mid,part_num):
        system.particles[i].type = 'B'

all = hoomd.group.all()
gA = hoomd.group.type(type = 'A', update=True)
gB = hoomd.group.type(type = 'B', update=True)
N = len(all)
Na = len(gA)
Nb = len(gB)

#define potential between pairs
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=2**(1/6), nlist=nl)
lj.set_params(mode='shift')
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('A', 'B', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0)

#integrator type
hoomd.md.integrate.mode_minimize_fire(group=all, dt=0.00001, ftol=1e-2, Etol=1e-7)
hoomd.run(10000)

#run simulation with current settings here
hoomd.md.integrate.mode_standard(dt=0.00001)
hoomd.md.integrate.brownian(group=all, kT=0.5, seed=123)
hoomd.run(100000)

#set the activity of each type
import numpy as np

angle = np.random.rand(part_num) * 2 * np.pi

if part_perc_a != 0 and part_perc_a != 100:
    activity_a = []
    for i in range(0,mid):
        x = (np.cos(angle[i])) * pe_a
        y = (np.sin(angle[i])) * pe_a
        z = 0
        tuple = (x, y, z)
        activity_a.append(tuple)
    activity_b = []
    for i in range(mid,part_num):
        x = (np.cos(angle[i])) * pe_b
        y = (np.sin(angle[i])) * pe_b
        z = 0
        tuple = (x, y, z)
        activity_b.append(tuple)
    hoomd.md.force.active(group=gA,
                          seed=123,
                          f_lst=activity_a,
                          rotation_diff=3.0,
                          orientation_link=False)
    hoomd.md.force.active(group=gB,
                          seed=375,
                          f_lst=activity_b,
                          rotation_diff=3.0,
                          orientation_link=False)
else:
    if part_perc_a == 0:
        activity_b = []
        for i in range(0,part_num):
            x = (np.cos(angle[i])) * pe_b
            y = (np.sin(angle[i])) * pe_b
            z = 0
            tuple = (x, y, z)
            activity_b.append(tuple)
        hoomd.md.force.active(group=gB,
                              seed=375,
                              f_lst=activity_b,
                              rotation_diff=3.0,
                              orientation_link=False)
    else:
        activity_a = []
        for i in range(0,part_num):
            x = (np.cos(angle[i])) * pe_a
            y = (np.sin(angle[i])) * pe_a
            z = 0
            tuple = (x, y, z)
            activity_a.append(tuple)
        hoomd.md.force.active(group=gA,
                              seed=123,
                              f_lst=activity_a,
                              rotation_diff=3.0,
                              orientation_link=False)
###
def dump_spec(timestep):
    if timestep in dump_list:
        hoomd.dump.gsd(filename="test.gsd", period=None, group=all, overwrite=False, static=[])

hoomd.analyze.callback(callback = dump_spec, period = 1)
###

#write dumps
name = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"
hoomd.dump.gsd(name, period=dump_freq, group=all, overwrite=True, static=[])

#run
hoomd.run(tsteps)
