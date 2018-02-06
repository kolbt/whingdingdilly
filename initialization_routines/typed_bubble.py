'''
SCRIPT SYNOPSIS: I want to look at bubbles in my simulations.  In order to perform
systematic studies, I need to be able to initialize these guys with some flexibility.

Here's how this will work:

           Draw smallest circle:  particle density WITHIN this circle is equal
                                  to gas density
                                  
    Draw second smallest circle:  inbetween this and the smallest circle
                                  particles should be HCP and orientationally
                                  aligned facing outward (away from (0, 0))
                                  
     Draw second largest circle:  between this and the second smallest circle
                                  particles should have random orientation
                                  but be HCP
                                  
            Draw largest circle:  between this and the second largest circle
                                  particles should be aligned facing inward
                                  (towards (0, 0)), and be HCP. Outside this
                                  circle particle density is equal to gas 
                                  density
                                  
                      Gas Phase:  after placing dense phase particles, gas
                                  phase particles will be placed in the 
                                  interior region randomly (to match the gas
                                  phase area fraction), and then randomly in
                                  the exterior gas
             
Immutable parameters:

    Area fraction:      held at 0.6
    Particle radius:    held at 0.5
    
User input:

    Positional Data
    R1:         radius of innermost circle
    R2:         radius of 2nd smallest circle
    R3:         radius of 2nd largest circle
    R4:         radius of largest circle
    box_width:  width of simulation box
    hcp_space:  spacing between particles in HCP phase
    
    Activity Data
    pe_a:       activity of species a
    pe_b:       activity of species b
    
Errors to throw:
    
    !!! R1 < R2 < R3 < R4 < BOX_WIDTH/2 !!!
    !!! part_num > part_dense !!!
'''

import sys
sys.path.append('/Users/kolbt/Desktop/compiled/hoomd-blue_11_8_17/hoomd-blue/build')
import hoomd
from hoomd import md
from hoomd import dem
from hoomd import deprecated
import numpy as np
from mypy import load_bar

# Formula to obtain hcp distances
def computePythag(hcp_space):
    a = hcp_space / 2
    c = hcp_space
    b = np.sqrt((c**2) - (a**2))
    return b

# Distance formula
def distTo(x1, y1, x2, y2):
    return np.sqrt(((x2-x1)**2)+((y2-y1)**2))

# Compute the total number of particles
def computePartNum(area_fraction, box_width, particle_radius):
    
    box_area = box_width ** 2
    all_particle_area = box_area * area_fraction
    one_particle_area = np.pi * (particle_radius**2)
    part_num = all_particle_area / one_particle_area
    
    return int(part_num) # This will always round down!

# Compute number of particles in dense/gas phases
def computePhaseNum(R4, R1, particle_radius, box_width, part_num):
    
    # Find area of each phase
    area_box = box_width**2
    area_big = np.pi * (R4**2)
    area_small = np.pi * (R1**2)
    area_dense = area_big - area_small
    one_particle_area = np.pi * (particle_radius**2)
    
    # Find: particles in dense phase/gas
    part_dense = int(area_dense / one_particle_area)
    part_gas = part_num - part_dense
    
    return int(part_dense), int(part_gas)

# Compute density of gas phase
def computeGasDensity(R4, R1, particle_radius, box_width, part_gas):
    
    # Find area of each phase
    area_box = box_width**2
    area_big = np.pi * (R4**2)
    area_small = np.pi * (R1**2)
    area_dense = area_big - area_small
    area_gas = area_box - area_dense
    one_particle_area = np.pi * (particle_radius**2)
    dense_gas = (float(part_gas) * one_particle_area) / area_gas

    return dense_gas

# NOTE: ORIENTATIONS ARE READ IN AS A 4-VECTOR QUATERNION IN HOOMD
#       MUST BE FLOAT or DOUBLE
# Output orientation for particle on outer edge
def orientIn(xyz, pe):
    x = xyz[0]
    y = xyz[1]
    z = 0
    
    # Point your vector inward
    if x < 0:
        x = abs(x)
    else:
        x *= -1
    if y < 0:
        y = abs(y)
    else:
        y *= -1

    # Normalize and set magnitude to active force
    mag = np.sqrt((x**2)+(y**2))
    x /= mag
    y /= mag
    x *= pe
    y *= pe
    
    orientation = (x, y, z)
    return orientation

# Output orientation for particle on inner edge
def orientOut(xyz, pe):
    x = xyz[0]
    y = xyz[1]
    z = 0
    
    mag = np.sqrt((x**2)+(y**2))
    x /= mag
    y /= mag
    x *= pe
    y *= pe
    
    orientation = (x, y, z)
    return orientation

# Output random orientation for gas and unaligned dense phase
def orientRand(pe):
    
    angle = np.random.rand() * 2 * np.pi
    x = np.cos(angle) * pe
    y = np.cos(angle) * pe
    z = 0
    
    orientation = (x, y, z)
    return orientation

# Immutable parameters:
area_fraction = 0.60
particle_radius = 0.50

leave = 1
while leave:
    # User inputs:
    R1 = float(raw_input("R1: "))                   # smallest circle
    R2 = float(raw_input("R2: "))                   # end outward alignment
    R3 = float(raw_input("R3: "))                   # start inward alignment
    R4 = float(raw_input("R4: "))                   # largest circle
    sep_x = float(raw_input("HCP spacing: " ))      # HCP spacing
    pe_a = float(raw_input("Activity A: " ))        # activity of species A
    pe_b = float(raw_input("Activity B: " ))        # activity of species B
    x_a = float(raw_input("Particle fraction A: " ))# particle fraction
    box_width = float(raw_input("Box width: "))
    box_lims = float(box_width / 2)

    # Some quick error checking for radii input
    if (R1 < R2 < R3 < R4) == 0:
        print "ERROR: Incorrect radii input"
        continue

    # Number of: total particles, dense phase particles, gas particles
    part_num =\
        computePartNum(area_fraction, box_width, particle_radius)
    part_dense, part_gas =\
        computePhaseNum(R4, R1, particle_radius, box_width, part_num)
    dense_gas =\
        computeGasDensity(R4, R1, particle_radius, box_width, part_gas)


    # Error checking for particle computations
    if (part_dense > part_num):
        print "ERROR: Too many particles in dense phase"
        continue

    else:
        leave = 0

print ""
print "ESTIMATED QUANTITIES:"
print "Total number of particles:", part_num
print "Dense phase particles:", part_dense
print "Gas phase particles:", part_gas
print "Gas phase density:", dense_gas
print ""

# Parameters have been set, now initialize the system:
hoomd.context.initialize()
snap = hoomd.data.make_snapshot(N = part_num,
                                box = hoomd.data.boxdim(L=box_width,
                                                        dimensions = 2),
                                particle_types = ['A','B'])

# Set up placement locations for dense phase:

start_y = -box_lims         # start placement in lower left corner of box
sep_y = computePythag(sep_x)
ith = 0                         # particle index
n = 0                           # column index (horizontal coord)
m = 0                           # row index (vertical coord)
row = 2                         # determines offset for x-coordinate
offset = 0.45                   # x offset between rows
activity_a = []                 # instantiate activity list (holds orientation)
activity_b = []                 # instantiate activity list for type b
part = np.zeros((3))            # temporarily holds (x,y,z) for each particle

# Place dense phase particles:
# Outer loop for y location
while (start_y + m) < box_lims:

    part[1] = start_y + m   # gives current y location
    n = 0                   # reset the x location
    # Ensures adjacent rows are offset in x dimension
    if row % 2 == 0:
        start_x = -box_lims
    else:
        start_x = -box_lims + offset

    # Inner loop for x location
    while (start_x + n) < box_lims:
        part[0] = start_x + n                       # gives current x location
        distance = distTo(part[0], part[1], 0, 0)   # compute distance
        
        if (R1 <= distance <= R4):                      # particle is in dense phase
            if (R1 <= distance <= R2):                  # particle is on inner edge
                activity_b.append(orientOut(part, pe_a))
                snap.particles.typeid[ith] = 1
            elif (R3 <= distance <= R4):                # particle is on outer edge
                activity_a.append(orientIn(part, pe_a))
                snap.particles.typeid[ith] = 0
            else:                                       # particle is in unaligned bulk
                random = np.random.rand()
                if random > x_a:                        # x_a gives particle bias
                    activity_b.append(orientRand(pe_b))
                    snap.particles.typeid[ith] = 1
                else:
                    activity_a.append(orientRand(pe_a))
                    snap.particles.typeid[ith] = 0
            
            # Place position for particle and increment position
            snap.particles.position[ith] = part
            
            ith += 1
            n += sep_x

        else:                                   # increment position
            n += sep_x

    row += 1                # increment x offset
    m += sep_y              # increment y position

part_dense = ith - 1
part_gas = part_num - part_dense
dense_gas =\
    computeGasDensity(R4, R1, particle_radius, box_width, part_gas)

'''
Evidently, the thickness of the orientationally aligned edge greatly effects
these computations for extreme values of x_a.  This script requires greater
care in particle placement in the dense phase should it be used in extreme
cases (x_a < 0.2 || x_a > 0.8).
'''

dense_A = len(activity_a)               # how many dense phase particles are type A
dense_B = len(activity_b)               # how many dense phase particles are type A
gas_A = int(x_a * part_num) - dense_A   # compute number of A in gas phase
gas_B = int(part_num -\
            dense_A -\
            dense_B -\
            gas_A)                      # leftover is gas_B

print ""
print "EXACT QUANTITIES:"
print "Total number of particles:", part_num
print "Dense phase particles:", part_dense
print "Dense phase A:", dense_A
print "Dense phase B:", dense_B
print "Gas phase particles:", part_gas
print "Gas phase A:", gas_A
print "Gas phase B:", gas_B
print "Gas phase density:", dense_gas
print ""

# Make a mesh for particle placement:
sigma = 2 * particle_radius                             # the particle diameter
grid_width = 2 * sigma                                  # gives width of each grid indx
spacer = int(box_width / (grid_width))                  # gives number of grid indices
mesh = np.zeros((spacer, spacer), dtype = np.ndarray)   # this will give grid occupation
occ = 20                                                # max occupancy (overestimated)
pos = np.zeros((occ, 2))                                # each index can hold occ coords

for iii in range(0, spacer):                            # occupancy list for each grid
    for jjj in range(0, spacer):
        mesh[iii][jjj] =  np.zeros_like(pos)

# Place remaining particles in gas phase:
while (ith < part_num):
    part = np.zeros(3, dtype=np.float64)
    part[0] = np.random.rand() * float(box_width)
    part[1] = np.random.rand() * float(box_width)
    x = part[0] - box_lims
    y = part[1] - box_lims

    # First make sure it isn't in the dense phase
    if (R1 - 0.5) < distTo(x, y, 0, 0) < (R4 + 0.5):
        continue
    
    # Get relevant indices and coordinates:
    loc_x = int(part[0] / grid_width)               # x index
    loc_y = int(part[1] / grid_width)               # y index
    lft_wrap = loc_x - 1                            # boxes to left
    rgt_wrap = loc_x + 1                            # boxes to right
    dwn_wrap = loc_y - 1                            # boxes below
    top_wrap = loc_y + 1                            # boxes above
    
    # Make sure indices are okay
    if rgt_wrap == spacer:
        rgt_wrap = 0
    if top_wrap == spacer:
        top_wrap = 0

    h = np.zeros((3), dtype = np.int)
    h[0], h[1], h[2] = lft_wrap, loc_x, rgt_wrap
    v = np.zeros((3), dtype = np.int)
    v[0], v[1], v[2] = dwn_wrap, loc_y, top_wrap

    # Check the grid you are placing the particle in and adjacent grids
    for iii in range(0, len(h)):
        for jjj in range(0, len(v)):
            for mmm in range(0, occ):
                if mesh[h[iii]][v[jjj]][mmm][0] != 0:
                    r = distTo(mesh[h[iii]][v[jjj]][mmm][0],
                               mesh[h[iii]][v[jjj]][mmm][1],
                               part[0],
                               part[1])
                    if r > 90:
                        r -= box_width
                        r = abs(r)
                    if r < 1.0:         # if too close, break from mmm loop
                        break
                else:                   # if no particle, use passable r
                    r = 2
            if r < 1.0:                 # break from jjj loop
                break
        if r < 1.0:                     # break from iii loop
            break
    if r < 1.0:                         # restart while statement
        continue

    # This particle has no violations
    for nnn in range(0, occ):
        if mesh[loc_x][loc_y][nnn][0] == 0:             # find unoccupied spot
            mesh[loc_x][loc_y][nnn][0] = part[0]        # add new x coord
            mesh[loc_x][loc_y][nnn][1] = part[1]        # add new y coord
            break

    snap.particles.position[ith] = (x, y, 0)            # place random pos
    ith += 1

# Randomly assign types and orientations to gas phase
for iii in range((part_dense + 1), part_num):
    if iii <= (part_dense + 1 + gas_A):
        snap.particles.typeid[iii] = 0
        activity_a.append(orientRand(pe_a))
    else:
        snap.particles.typeid[iii] = 1
        activity_b.append(orientRand(pe_b))

# All particles have now been placed, run simulation
system = hoomd.init.read_snapshot(snap)

all = hoomd.group.all()
gA = hoomd.group.type(type = 'A', update=True)
gB = hoomd.group.type(type = 'B', update=True)

# Pair potentials
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=2**(1/6), nlist=nl)
lj.set_params(mode='shift')
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('A', 'B', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0)

tau = 1
my_dt = 0.000001 * tau

hoomd.md.integrate.mode_standard(dt=my_dt)
hoomd.md.integrate.brownian(group=all, kT=1.0, seed=1345)

hoomd.md.force.active(group=gA,
                      seed=2034,
                      f_lst=activity_a,
                      rotation_diff=3.0,
                      orientation_link=False,
                      orientation_reverse_link=True)

hoomd.md.force.active(group=gB,
                      seed=9348,
                      f_lst=activity_b,
                      rotation_diff=3.0,
                      orientation_link=False,
                      orientation_reverse_link=True)

hoomd.dump.gsd('bubble_pea_150_peb_500_xa50.gsd',
               period = 10000,
               group = all,
               overwrite = True,
               dynamic = ['attribute', 'property', 'momentum'])

hoomd.run(1000000)


######################################################################################
#                                                                                    #
#                       IMPORTANT CHANGES TO SCRIPT                                  #
#                                                                                    #
# Means of re-initilizing a system:                                                  #
#   write out initial positions and orientations to restart file                     #
#                                                                                    #
#                                                                                    #
#                                                                                    #
#                                                                                    #
#                                                                                    #
#                                                                                    #
#                                                                                    #
#                                                                                    #
#                                                                                    #
#                                                                                    #
#                                                                                    #
#                                                                                    #
#                                                                                    #
#                                                                                    #
######################################################################################

