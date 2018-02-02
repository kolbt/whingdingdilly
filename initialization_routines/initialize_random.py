import sys
sys.path.append('/Users/kolbt/Desktop/compiled/hoomd-blue/build')
import hoomd
from hoomd import md
from hoomd import dem
from hoomd import deprecated
import numpy as np

tsteps = 500000
dump_freq = 1000
part_perc_a = 50
part_frac_a = float(part_perc_a) / float(100)
pe_a = 0
pe_b = 150
phi = 0.5
part_num = 15000
dumps = tsteps/dump_freq
diameter = 1

# find the box parameters
area_part = np.pi * ((float(diameter)/float(2))**2) * part_num
box_area = area_part / phi
side = int(np.sqrt(box_area))
while side % 10 != 0:                    # this is sub par... fix it
    side += 1                           # or just pick part_num so that this is okay

# find the number of A and B particles
part_a = part_num * part_frac_a         # get the total number of A particles
part_b = part_num - part_a              # get the total number of B particles
mid = int(part_a) - 1                   # starting point for assigning B particles

# write a distance formula function
def computeDistance(arg1, arg2):
    r = np.sqrt((arg2[0] - arg1[0])**2+(arg2[1] - arg1[1])**2)
    return r

# initialize system randomly
hoomd.context.initialize()

snap = hoomd.data.make_snapshot(N = part_num,
                                box = hoomd.data.boxdim(L=side,
                                                        dimensions=2),
                                particle_types = ['A', 'B'])

# make a mesh for my simulation box
# keep track of arrays by their index
spacer = side / (5 * diameter)                          # spacing between indiviudal boxes
mesh = np.zeros((spacer, spacer), dtype = np.ndarray)   # array of each grid
#pos = np.zeros((part_num, 2))
occ = 40                                                # max occupancy of a grid box
pos = np.zeros((occ, 2))
for j in range(0, spacer):
    for k in range(0, spacer):
        mesh[j][k] =  np.zeros_like(pos)

ind_cnt = mesh.shape[0] - 1                                 # final index in grid array

# generate positions for particles
# CHECK R WITH SURROUNDING LIST IF CLOSE TO BORDER
for i in range(0, part_num):
    part = np.zeros(3)                                      # tmp particle position
    while part[0] == 0 and part[1] == 0:
        
        part[0] = np.random.rand() * side                   # generate random x coord
        part[1] = np.random.rand() * side                   # generate random y coord
        part[2] = 0
        loc_x = int(part[0] / (5 * diameter))               # find x mesh index
        loc_y = int(part[1] / (5 * diameter))               # find y mesh index
        lft_bound = 5 * diameter * loc_x                    # coord of grid boundaries
        rgt_bound = lft_bound + 5 * diameter
        dwn_bound = 5 * diameter * loc_y
        top_bound = loc_y + 5 * diameter
        lft_chk = lft_bound + 1.0                           # additional check condition
        rgt_chk = rgt_bound - 1.0                           # if one away from index edge check neighbor
        dwn_chk = dwn_bound + 1.0
        top_chk = top_bound - 1.0
        lft_wrap = loc_x - 1
        rgt_wrap = loc_x + 1
        dwn_wrap = loc_y - 1
        top_wrap = loc_y + 1
        
        if loc_x == 0:                                      # wrap the system (periodic)
            lft_wrap = ind_cnt
        if loc_x == ind_cnt:
            rgt_wrap = 0
        if loc_y == 0:
            dwn_wrap = ind_cnt
        if loc_y == ind_cnt:
            top_wrap = 0

        for j in range(0, occ):
            if mesh[loc_x][loc_y][j][0] != 0:               # if non-zero entry compute distance
                r = computeDistance(mesh[loc_x][loc_y][j], part)
                if r < 1.0:                                 # if too close, do not place
                    part[0] = 0                             # reset position, re-loop
                    part[1] = 0
                    break

        if part[0] != 0 and part[1] != 0:                   # this checks adjacent grids if necessary
            for k in range(0, occ):
                if part[0] < lft_chk:                       # check other grids if too close
                    if mesh[lft_wrap][loc_y][k][0] != 0:
                        r = computeDistance(mesh[lft_wrap][loc_y][k], part)
                        if r < 1.0:
                            part[0] = 0
                            part[1] = 0
                            break
                if part[0] > rgt_chk:
                    if mesh[rgt_wrap][loc_y][k][0] != 0:
                        r = computeDistance(mesh[rgt_wrap][loc_y][k], part)
                        if r < 1.0:
                            part[0] = 0
                            part[1] = 0
                            break
                if part[1] < dwn_chk:
                    if mesh[loc_x][dwn_wrap][k][0] != 0:
                        r = computeDistance(mesh[loc_x][dwn_wrap][k], part)
                        if r < 1.0:
                            part[0] = 0
                            part[1] = 0
                            break
                if part[1] > top_chk:
                    if mesh[loc_x][top_wrap][k][0] != 0:
                        r = computeDistance(mesh[loc_x][top_wrap][k], part)
                        if r < 1.0:
                            part[0] = 0
                            part[1] = 0
                            break
                else:
                    break
                                
        if part[0] != 0 and part[1] != 0:                   # particle has no violations
            for s in range(0, occ):
                if mesh[loc_x][loc_y][s][0] == 0:
                    mesh[loc_x][loc_y][s][0] = part[0]      # add to the mesh
                    mesh[loc_x][loc_y][s][1] = part[1]
                    break
            if part[0] > (side/2):                          # compute real position
                part[0] -= side
            if part[1] > (side/2):
                part[1] -= side
            snap.particles.position[i] = part               # set particle position in hoomd
            print(snap.particles.position[i])
            print(i)


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

# set particle types
if part_perc_a == 0:                    # all b
    for i in range(0, part_num):
        snap.particles.typeid[i] = 1

elif part_perc_a != 100:                # mix of each
    for i in range(mid, part_num):
        snap.particles.typeid[i] = 1

else:
    for i in range(0, part_num):        # all a
        snap.particles.typeid[i] = 0

####################################
### NOW SET FORCES / INTEGRATORS ###
####################################

# initialize the system
system = hoomd.init.read_snapshot(snap)

all = hoomd.group.all()

nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=2**(1/6), nlist=nl)
lj.set_params(mode='shift')
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('A', 'B', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0)

hoomd.md.force.active(group=all,
                      seed=123,
                      f_lst=part_num * [(0, 0, 0)],
                      rotation_diff=3,
                      orientation_link=False)


hoomd.dump.gsd("spheres.gsd", period=1000, group=all, overwrite=True)


# minimize for no overlaps
fire=hoomd.md.integrate.mode_minimize_fire(group=all,
                                           dt=0.0000000001,
                                           ftol=1e-2,
                                           Etol=1e-7)
hoomd.run(1000000)

# brownian integration
hoomd.md.integrate.mode_standard(dt=0.00002, aniso=True)
bd = hoomd.md.integrate.brownian(group=all, kT=1.0, seed=123)
bd.set_gamma('A', gamma=1.0)
bd.set_gamma_r('A', gamma_r=1.0)

#write dump
hoomd.dump.gsd("spheres.gsd", period=1000, group=all, overwrite=True)

#run
hoomd.run(50000)
