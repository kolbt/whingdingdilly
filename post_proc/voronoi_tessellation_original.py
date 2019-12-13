'''
#                           This is an 80 character line                       #
What does this file do?

1.) Read in .gsd file of particle positions
2.) Feed data into Freud for Voronoi tessellation
3.) Output image of the tessellation for each tstep

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

import freud
from freud import parallel
from freud import box
from freud import voronoi

import numpy as np
from scipy import stats

import matplotlib
import matplotlib.pyplot as plt

### PLOTTER FUNCTION ###
def draw_voronoi(box, points, cells, nlist=None, color_by_sides=True):
    ax = plt.gca()
    # Draw Voronoi cells
    patches = [plt.Polygon(cell[:, :2]) for cell in cells]
    patch_collection = matplotlib.collections.PatchCollection(patches, edgecolors='black', alpha=0.4)
    cmap = plt.cm.Set1
    
    if color_by_sides:
        colors = [len(cell) for cell in voro.polytopes]
    else:
        colors = np.random.permutation(np.arange(len(patches)))

    cmap = plt.cm.get_cmap('Set1', np.unique(colors).size)
    bounds = np.array(range(min(colors), max(colors)+2))

    patch_collection.set_array(np.array(colors))
    patch_collection.set_cmap(cmap)
    patch_collection.set_clim(bounds[0], bounds[-1])
    ax.add_collection(patch_collection)

    # Draw points
    plt.scatter(points[:,0], points[:,1], c=colors)
    plt.title('Voronoi Diagram')
    plt.xlim((-box.Lx/2, box.Lx/2))
    plt.ylim((-box.Ly/2, box.Ly/2))
    
    # Set equal aspect and draw box
    ax.set_aspect('equal', 'datalim')
    box_patch = plt.Rectangle([-box.Lx/2, -box.Ly/2], box.Lx, box.Ly, alpha=1, fill=None)
    ax.add_patch(box_patch)
    # Draw neighbor lines
    if nlist is not None:
        bonds = np.asarray([points[j] - points[i] for i, j in zip(nlist.index_i, nlist.index_j)])
        box.wrap(bonds)
        line_data = np.asarray([[points[nlist.index_i[i]], points[nlist.index_i[i]]+bonds[i]] for i in range(len(nlist.index_i))])
        line_data = line_data[:, :, :2]
        line_collection = matplotlib.collections.LineCollection(line_data, alpha=0.3)
        ax.add_collection(line_collection)
    # Show colorbar for number of sides
    if color_by_sides:
        cb = plt.colorbar(patch_collection, ax=ax, ticks=bounds, boundaries=bounds)
        cb.set_ticks(cb.formatter.locs + 0.5)
        cb.set_ticklabels((cb.formatter.locs - 0.5).astype('int'))
        cb.set_label("Number of sides", fontsize=12)
    plt.show()
### END PLOTTER FUNCTION ###

def slowSort(array):
    """Sort an array the slow (but certain) way"""
    cpy = np.copy(array)
    ind = np.arange(0, len(array))
    for i in xrange(len(cpy)):
        for j in xrange(len(cpy)):
            if cpy[i] > cpy[j] and i < j:
                # Swap the copy array values
                tmp = cpy[i]
                cpy[i] = cpy[j]
                cpy[j] = tmp
                # Swap the corresponding indices
                tmp = ind[i]
                ind[i] = ind[j]
                ind[j] = tmp
    return ind

def indSort(arr1, arr2):
    """Take sorted index array, use to sort array"""
    # arr1 is array to sort
    # arr2 is index array
    cpy = np.copy(arr1)
    for i in xrange(len(arr1)):
        arr1[i] = cpy[arr2[i]]

def chkSort(array):
    """Make sure sort actually did its job"""
    for i in xrange(len(array)-2):
        if array[i] > array[i+1]:
            print("{} is not greater than {} for indices=({},{})").format(array[i+1], array[i], i, i+1)
            return False
    return True

try:
    # This is for the fine timescale data
    gsd_file = "pa" + str(pe_a) + \
                 "_pb" + str(pe_b) + \
                 "_xa" + str(part_perc_a) + \
                 ".gsd"
    # File to write all data to
    out_file = "voronoi_pa" + str(pe_a) + \
               "_pb" + str(pe_b) + \
               "_xa" + str(part_perc_a)
    g = hoomd.open(name=gsd_file, mode='rb') # open gsd file with hoomd
except:
    try:
        eps = str(sys.argv[6])
    except:
        eps = 1
    gsd_file = "pa" + str(pe_a) + \
                 "_pb" + str(pe_b) + \
                 "_xa" + str(part_perc_a) + \
                 "_ep" + str(eps) + \
                 ".gsd"
    # File to write all data to
    out_file = "voronoi_pa" + str(pe_a) + \
               "_pb" + str(pe_b) + \
               "_xa" + str(part_perc_a) + \
               "_ep" + str(eps)
    g = hoomd.open(name=gsd_file, mode='rb')  # open gsd file with hoomd

dumps = int(g.__len__())               # get number of timesteps dumped

start = 0                       # gives first frame to read
end = dumps                     # gives last frame to read

positions = np.zeros((end), dtype=np.ndarray)       # array of positions
types = np.zeros((end), dtype=np.ndarray)           # particle types
box_data = np.zeros((1), dtype=np.ndarray)          # box dimensions
timesteps = np.zeros((end), dtype=np.float64)       # timesteps

# Get relevant data from short.gsd file
with hoomd.open(name=gsd_file, mode='rb') as t:
    snap = t[0]
    box_data = snap.configuration.box
    for iii in range(start, end):
        snap = t[iii]                               # snapshot of frame
        types[iii] = snap.particles.typeid          # get types
        positions[iii] = snap.particles.position    # get positions
        timesteps[iii] = snap.configuration.step    # get timestep

# Get index order for chronological timestep sorting
newInd = slowSort(timesteps)
# Use these indexes to reorder other arrays
indSort(timesteps, newInd)
indSort(positions, newInd)
indSort(types, newInd)

if chkSort(timesteps):
    print("Array succesfully sorted")
else:
    print("Array not sorted")
timesteps -= timesteps[0]

# Get number of each type of particle
partNum = len(types[start])
part_A = int(partNum * part_frac_a)
part_B = partNum - part_A

# Feed data into freud analysis software
l_box = box_data[0]
h_box = l_box / 2.0
a_box = l_box * l_box
f_box = box.Box(Lx = l_box, Ly = l_box, is2D = True)    # make freud box
voro = voronoi.Voronoi(box=f_box, buff=(l_box / 2.0))

for j in range(end - 2, end):
#for j in range(start, end):
    # Easier accessors
    pos = positions[j]
    # Zero out the z component
    pos[:,-1] = 0.0
    # pos = np.delete(pos, 2, 1)
    typ = types[j]
    tst = timesteps[j]
    
    # Compute the Voronoi tessellation
    voro.compute(box=f_box, positions=pos, buff=h_box)
    draw_voronoi(f_box, pos, voro.polytopes)
