import sys
sys.path.append('/Users/kolbt/Desktop/compiled/hoomd-blue/build')
import hoomd
from hoomd import md
from hoomd import dem
from hoomd import deprecated
import numpy as np

# Simulation box mesh into grid delimit by particle diameter
# list of mesh indices random number generator to select index
# remove index from list once particle is placed

tsteps = 5000000
dump_freq = 10000
part_perc_a = 50
part_frac_a = float(part_perc_a) / float(100)
pe_a = 80
pe_b = 300
phi = 0.6
part_num = 24102
dumps = tsteps/dump_freq
diameter = 1

# find the box parameters
area_part = np.pi * ((float(diameter)/float(2))**2) * part_num
box_area = area_part / phi
side = int(np.sqrt(box_area))
side = 140
#while side % 10 != 0:                    # this is sub par... fix it
#side += 1                           # or just pick part_num so that this is okay

# initialize system randomly
hoomd.context.initialize()
part_num = 13950
part_a = part_num * part_frac_a         # get the total number of A particles
part_a = int(part_a)
part_b = part_num - part_a              # get the total number of B particles
mid = int(part_a)                       # starting point for assigning B particles
snap = hoomd.data.make_snapshot(N = part_num,
                                box = hoomd.data.boxdim(L=side,
                                                        dimensions=2),
                                particle_types = ['A', 'B'])

part = np.zeros((3))


start_y = -69.5     # box is -70:70 for x and y dimensions
sep_row = 0.90      # distance between particles along x axis
sep_col = 0.78      # distance to increment rows (maintains center to center distance)
ith = 0             # particle counter
m = 0               # incrementer for y value
row = 2             # start on an even row (this determines first x placement in row)

# Places particles in lower left quadrant (-70, -70) - (0, 0)
# while loop that increments y value
while 1:
    part[0] = start_y + m
    n = 0
    # while that increments x value (place row at constant height, y value)
    while 1:
        # ensures rows are offset from one another
        if row % 2 == 0:
            start_x = -69.50
        else:
            start_x = -69.05
        
        part[1] = start_x + n
        snap.particles.position[ith] = part
        snap.particles.typeid[ith] = 0
        ith += 1
        n += sep_row
        # placing into lower left quadrant
        if start_x + n > 0:
            break
                
    row += 1
    m += sep_col
    # ensure particles are limited to lower left quadrant
    if -69.5 + m > 0:
        break

# Places particles in upper right quadrant (0,0) - (70, 70)
m = 0
row = 2
start_y = 0.5
while 1:
    part[0] = start_y + m
    n = 0
    while 1:
        
        if row % 2 == 0:
            start_x = 0.5
        else:
            start_x = 0.95
        
        part[1] = 0.5 + n
        snap.particles.position[ith] = part
        snap.particles.typeid[ith] = 1
        ith += 1
        n += sep_row
        
        if start_x + n > 70:
            break

    row += 1
    m += sep_col
    if start_y + m > 70:
        break

print(ith)
print(ith)

# now let's get the quaternion and moment of inertia
thetas = np.random.uniform(0, 2*np.pi, (part_num,))         # generate random angles
quats = np.array([np.cos(thetas/2),
                  np.zeros_like(thetas),
                  np.zeros_like(thetas),
                  np.sin(thetas/2)]).T                      # generate quaternions from the angles
snap.particles.orientation[:] = quats

inertia = float(1)/float(16)
snap.particles.diameter[:] = 1                              # set particle diameters
snap.particles.moment_inertia[:] = (inertia, inertia, 0)    # set moment of inertia
snap.particles.types = ['A', 'B']                           # or 0, 1 in typeid vernacular

####################################
### NOW SET FORCES / INTEGRATORS ###
####################################

# initialize the system
system = hoomd.init.read_snapshot(snap)

all = hoomd.group.all()
gA = hoomd.group.type(type = 'A', update=True)
gB = hoomd.group.type(type = 'B', update=True)
N = len(all)
part_num = N
Na = len(gA)
Nb = len(gB)

print(part_num)

nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=2**(1/6), nlist=nl)
lj.set_params(mode='shift')
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('A', 'B', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0)

angle = np.random.rand(part_num) * 2 * np.pi   # random orientation of each particle

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

# minimize for no overlaps
fire=hoomd.md.integrate.mode_minimize_fire(group=all,
                                           dt=0.00001,
                                           ftol=1e-2,
                                           Etol=1e-7)
hoomd.run(1000)

# brownian integration
hoomd.md.integrate.mode_standard(dt=0.000002)
bd = hoomd.md.integrate.brownian(group=all, kT=0.5, seed=123)
bd.set_gamma('A', gamma=1.0)
bd.set_gamma_r('A', gamma_r=1.0)

#write dump
hoomd.dump.gsd("hcp_test.gsd", period=1000, group=all, overwrite=True, static=[])

#run
hoomd.run(tsteps)
