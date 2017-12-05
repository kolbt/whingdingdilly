'''
I want to compute the voronoi diagram to find nearest neighbors:
    1. Get your points (from a GSD snapshot, ONLY DENSE PHASE)
    2. Draw circles of maximum radius 2 sigma
    3. If circle contains more than 3 points, delete
        (you now have all triangles)
    4. Draw line segment from middle of each side to convergence of perpendicular
        bisectors of each triangle (one point in each triangle)
    5. Hypothetically, you now have a vornoi diagram...
'''

import sys

hoomd_path = str(sys.argv[4])
gsd_path = str(sys.argv[5])

# need to extract values from filename (pa, pb, xa) for naming
part_perc_a = int(sys.argv[3])
part_frac_a = float(part_perc_a) / 100.0
pe_a = int(sys.argv[1])
pe_b = int(sys.argv[2])

sys.path.append(hoomd_path)

import hoomd
from hoomd import md
from hoomd import deprecated
from mypy import load_bar

sys.path.append(gsd_path)
import gsd
from gsd import hoomd
from gsd import pygsd
import numpy as np

def computeDistance(point1, point2):
    distance = np.sqrt((point2[0]-point1[0])**2+(point2[1]-point1[1])**2)
    return distance

# This finds a circle defined by any 3 points (based off of this overflow discussion):
# https://stackoverflow.com/questions/28910718/give-3-points-and-a-plot-circle

def computeCircle(part1, part2, part3):
    x, y, z = part1[0]+(part1[1]*1j), part2[0]+(part2[1]*1j), part3[0]+(part3[1]*1j)
    w = z-x
    w /= y-x
    c = (x-y) * (w-abs(w)**2) / 2j / w.imag - x
    #print '(x%+.3f)^2+(y%+.3f)^2 = %.3f^2' % (c.real, c.imag, abs(c+x))
    return (-c.real, -c.imag), abs(c+x)

def lineCheck(part1, part2, part3):
    slope1 = (part1[1] - part2[1]) * (part1[0] - part3[0])
    slope2 = (part1[1] - part3[1]) * (part1[0] - part2[0])
    if slope1 == slope2:
        return 0
    else:
        return 1

# Reading in the GSD file
myfile = "pa" + str(pe_a) + "_pb" + str(pe_b) + "_xa" + str(part_perc_a) + ".gsd"
f = hoomd.open(name=myfile, mode='rb')
dumps = f.__len__()
position_array = np.zeros((1), dtype=np.ndarray)        # array of position arrays
box_data = np.zeros((1), dtype=np.ndarray)              # box dimensions

with hoomd.open(name=myfile, mode='rb') as t:           # open for reading
    snap = t[100]                                       # snap 100th snapshot
    box_data = snap.configuration.box                   # get box dimensions
    position_array[0] = snap.particles.position         # store all particle positions

part_num = len(position_array[0])
points = np.zeros((part_num), dtype=np.ndarray)
points = position_array[0]
keepers = []
diameter = 1

##########################################################################################

### Compute neighbor lists

# create mesh
float_side = box_data[0] / 2.0
side = float((int(box_data[0])+1)) / 2
box_width = 1
while side % box_width != 0:                            # grid boxes of width 5 sigma
    side += 1                                           # make sure side length is divisible

spacer = int(side * 2 / (box_width * diameter))         # spacing between indiviudal boxes
mesh = np.zeros((spacer, spacer), dtype = np.ndarray)   # array of each grid
occ = 100                                               # max occupancy of a grid box
test_occ = np.zeros((occ, 4))                           # occupation test index
for j in range(0, spacer):
    for k in range(0, spacer):
        mesh[j][k] =  np.zeros_like(test_occ)

final_index = mesh.shape[0] - 1                         # final index in mesh

# put particles in appropriate mesh grid
for iii in range(part_num):
    # get the index of the mesh the particle belongs in
    loc_x = int((points[iii][0] + float_side) / (box_width * diameter))
    loc_y = int((points[iii][1] + float_side) / (box_width * diameter))
    
    # place the particle in the first unoccupied space in this quadrant list
    for s in range(1, occ):                             # index 0 holds total number in grid
        if mesh[loc_x][loc_y][s][3] == 0:               # test occupancy of list
            mesh[loc_x][loc_y][s][0] = points[iii][0]   # x coord
            mesh[loc_x][loc_y][s][1] = points[iii][1]   # y coord
            mesh[loc_x][loc_y][s][2] = iii              # particle id
            mesh[loc_x][loc_y][s][3] = 1                # switch flag to occupied
            break

# count the number of particles in each grid
for iii in range(spacer):
    for jjj in range(spacer):
        count = 0
        for hhh in range(1,occ):
            if mesh[iii][jjj][hhh][0] != 0:
                count += 1
            else:
                mesh[iii][jjj][0][0] = count    # number of particles: first index of each box
                break

### Compute circles
load_bar.printLoadBar(0, part_num, prefix = "Progress:", suffix = "Complete")
for iii in range(0, part_num):
    load_bar.printLoadBar(iii+1, part_num, prefix = "Progress:", suffix = "Complete")
    
    loc_x = int((points[iii][0] + float_side) / (box_width * diameter))
    loc_y = int((points[iii][1] + float_side) / (box_width * diameter))
    # indices of surrounding grids
    lft = loc_x - 1
    rgt = loc_x + 1
    bot = loc_y - 1
    top = loc_y + 1
    # compute total number of particles / build a list of surrounding particles
    near = []
    for jjj in range(lft, rgt+1):
        for hhh in range(bot, top+1):
            if jjj >= spacer:
                jjj -= spacer
            if hhh >= spacer:
                hhh -= spacer
            for mmm in range(1, occ):
                if mesh[jjj][hhh][mmm][0] == points[iii][0]:    # can't use itself
                    continue
                elif mesh[jjj][hhh][mmm][0] != 0:               # finds an entry
                    near.append(mesh[jjj][hhh][mmm][2])         # append particle id to nearest list
                else:                                           # doesn't find an entry
                    break

    for jjj in range(len(near)):
        for hhh in range(len(near)):
            
            # steps to avoid drawing the same circle twice
            if int(near[jjj]) <= iii:        # particle id MUST be greater than ref
                continue
            elif int(near[hhh]) <= iii:      # particle id MUST be greater than ref
                continue
            elif hhh <= jjj:                 # inner loop id MUST be greater than outter
                continue
            
            # good to go! Compute the circle
            elif lineCheck(points[iii],
                         points[int(near[jjj])],
                         points[int(near[hhh])]):
                center, radius = computeCircle(points[iii],
                                               points[int(near[jjj])],
                                               points[int(near[hhh])])
                circle_check = 0
                if radius < 1.0:
                    # make sure the circle doesn't contain a point
                    circle_check = 1
                    for mmm in range(len(near)):
                        if int(near[mmm]) == iii:   # first circle-draw point
                            continue
                        elif mmm == jjj:            # second circle-draw point
                            continue
                        elif mmm == hhh:            # third circle-draw point
                            continue
                        elif computeDistance(center, points[int(near[mmm])]) <= radius:
                            circle_check = 0
                            break
                if circle_check:
                    keepers.append((center, radius))

### Plotting
import matplotlib.pyplot as plt

# get points in xs and ys for scatter
xs = np.zeros((part_num), dtype=np.float32)
ys = np.zeros((part_num), dtype=np.float32)
for iii in range(0, part_num):
    xs[iii] = points[iii][0]
    ys[iii] = points[iii][1]

# plot circles
print(len(keepers))
fig, ax = plt.subplots()
load_bar.printLoadBar(0, len(keepers), prefix = "Progress:", suffix = "Complete")
for iii in range(0, len(keepers)):
    load_bar.printLoadBar(iii+1, len(keepers), prefix = "Progress:", suffix = "Complete")
    circle = plt.Circle(keepers[iii][0], keepers[iii][1], linewidth=0.15, fill=False)
    ax.add_artist(circle)
# plot points
plt.scatter(xs, ys, s=0.15)
ax.set_aspect('equal')
ax.set_xlim(-side,side)
ax.set_ylim(-side,side)
plt.savefig('vornoi_diagram.png', dpi=1000)
plt.close()

'''
IMPROVEMENTS:
    -feed in ONLY the dense phase
    -radius checking shouldn't be necessary
        loc_x = int((center[0] + float_side) / (box_width * diameter))
        loc_y = int((center[1] + float_side) / (box_width * diameter))
        then look for > 3 points within dist <= radius 
'''
