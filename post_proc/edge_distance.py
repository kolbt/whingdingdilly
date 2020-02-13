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

from shapely.ops import cascaded_union, polygonize
from scipy.spatial import Delaunay
from descartes import PolygonPatch
import shapely.geometry as geometry
from matplotlib.collections import LineCollection

import numpy as np
import math

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

def computeTauPerTstep(epsilon, mindt=0.000001):
    '''Read in epsilon, output tauBrownian per timestep'''
    if epsilon != 1.:
        mindt=0.00001
    kBT = 1.0
    tstepPerTau = float(epsilon / (kBT * mindt))
    return 1. / tstepPerTau

def plot_polygon(polygon):
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
#    ax.set_xlim([-h_box, h_box])
#    ax.set_ylim([-h_box, h_box])
    patch = PolygonPatch(polygon, fc='#999999',
                         ec='#000000', fill=True,
                         zorder=-1)
    ax.add_patch(patch)
    return fig

def alpha_shape(points, alpha):
    """
    Compute the alpha shape (concave hull) of a set
    of points.
    @param points: Iterable container of points.
    @param alpha: alpha value to influence the
        gooeyness of the border. Smaller numbers
        don't fall inward as much as larger numbers.
        Too large, and you lose everything!
    """
    if len(points) < 4:
        # When you have a triangle, there is no sense
        # in computing an alpha shape.
        return geometry.MultiPoint(list(points)).convex_hull
    def add_edge(edges, edge_points, coords, i, j):
        """
        Add a line between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add( (i, j) )
        edge_points.append(coords[ [i, j] ])

    coords = np.array([point for point in points])
    tri = Delaunay(coords)
    edges = set()
    edge_points = []
    # loop over triangles:
    # ia, ib, ic = indices of corner points of the
    # triangle
    for ia, ib, ic in tri.vertices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]
        # Lengths of sides of triangle
        a = math.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
        b = math.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
        c = math.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)
        # Semiperimeter of triangle
        s = (a + b + c)/2.0
        # Area of triangle by Heron's formula
        area = math.sqrt(s*(s-a)*(s-b)*(s-c))
        circum_r = a*b*c/(4.0*area)
        # Here's the radius filter.
        #print circum_r
        if circum_r < 1.0/alpha:
            add_edge(edges, edge_points, coords, ia, ib)
            add_edge(edges, edge_points, coords, ib, ic)
            add_edge(edges, edge_points, coords, ic, ia)
    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    return cascaded_union(triangles), edge_points

# Get infile and open
inFile = str(sys.argv[1])
if inFile[0:7] == "cluster":
    add = 'cluster_'
else:
    add = ''
    
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

# Outfile to write data to
outFile = add + 'edge_pa' + str(peA) +\
          '_pb' + str(peB) +\
          '_xa' + str(parFrac) +\
          '_phi' + str(intPhi) +\
          '_ep' + '{0:.3f}'.format(eps)

start = 0                   # first frame to process
dumps = int(f.__len__())    # get number of timesteps dumped
end = dumps                 # final frame to process
start = dumps - 1

box_data = np.zeros((1), dtype=np.ndarray)  # box dimension holder
r_cut = 2**(1./6.)                          # potential cutoff
tauPerDT = computeTauPerTstep(epsilon=eps)  # brownian time per timestep

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
        system = freud.AABBQuery(f_box, f_box.wrap(pos))
#        # Compute neighbor list for entire system
#        lc = freud.locality.AABBQuery(box=f_box, points=pos)
#        nlist = lc.query(pos, dict(num_neighbors=6, exclude_ii=True)).toNeighborList()
#        nlist.filter_r(r_max=1.0, r_min=0.2)
#        neigh = nlist.neighbor_counts
#        print(min(neigh))
#        print(max(neigh))
#        x = pos[:, 0]
#        y = pos[:, 1]
#        plt.scatter(x, y, c=list(neigh), s=0.5, edgecolor='none')
#        ax = plt.gca()
#        ax.axis('equal')
#        plt.show()
        
        # Compute neighbor list for only largest cluster
        my_clust.compute(system, neighbors={'r_max': 1.0})
        ids = my_clust.cluster_idx              # get id of each cluster
        c_props.compute(system, ids)            # find cluster properties
        clust_size = c_props.sizes              # find cluster sizes

        # Try and grab edge of largest cluster
        lClust = max(clust_size)
        # Get the id of the largest cluster
        for k in range(0, len(clust_size)):
            if clust_size[k] == lClust:
                lcID = ids[k]

        # Get the positions of all particles in LC
        lcPos = []
        for k in range(0, len(ids)):
            if ids[k] == lcID:
                lcPos.append(pos[k])
                
        # Feed LC into neighbor computation, extract low-neighbor points
        lc = freud.locality.AABBQuery(box=f_box, points=f_box.wrap(lcPos))
        nlist = lc.query(lcPos, dict(num_neighbors=6, exclude_ii=True)).toNeighborList()
        nlist.filter_r(r_max=1.0, r_min=0.2)
        neigh = nlist.neighbor_counts
        
        # YOU SHOULD BE ABLE TO USE SHAPELY ON YOUR INITIAL NEIGHBORLIST
        
        # Let's throw out any point that has >3 neighbors
        edPos = []
        for k in range(0, len(neigh)):
            if neigh[k] <= 4:
                edPos.append(lcPos[k])
                
        # Remove z-component of edPos
        x = list(list(zip(*edPos))[0])
        y = list(list(zip(*edPos))[1])
        edgePos = list(zip(x, y))
                
        # Let's use shapely and convex hulls
        concave_hull, edge_points = alpha_shape(edgePos, alpha=0.15)
        plot_polygon(concave_hull)
#        plot_polygon(concave_hull.buffer(1.5))
        plt.scatter(x, y, s=0.5, c='#FF6103', edgecolors='none')
        ax = plt.gca()
        ax.axis('equal')
        plt.show()
        
#        # THIS WORKS VERY WELL, BUT I CANNOT ACCESS EDGE LENGTH
#        # Let's try outsourcing this to a module
#        from polylidar import extractPlanesAndPolygons, extractPolygons, Delaunator
#        from polylidarutil import plot_polygons, get_estimated_lmax
#        lmax = 2.0
#        kwargs = dict(alpha=0.0, lmax=lmax)
#        print(np.asarray(lcPos).shape)
#        delaunay, planes, polygons = extractPlanesAndPolygons(np.asarray(lcPos), **kwargs)
#        print(polygons[0].length)
#        fig = plt.figure()
#        ax = plt.gca()
#        x = list(list(zip(*lcPos))[0])
#        y = list(list(zip(*lcPos))[1])
#        z = list(list(zip(*lcPos))[2])
#        ax.scatter(x, y, s=0.2, edgecolor='none')
#        plot_polygons(polygons, delaunay, np.asarray(lcPos), ax=ax)
#        ax.axis('equal')
#        plt.show()
        

#        # Compute the neighbor list for this subset cluster
#        lc = freud.locality.AABBQuery(box=f_box, points=lcPos)
##        nlist = lc.query(lcPos, dict(r_max=1.)).toNeighborList()
#        nlist = lc.query(lcPos, dict(num_neighbors=6, exclude_ii=True)).toNeighborList()
#        nlist.filter_r(r_max=1.0, r_min=0.2)
#        neigh = nlist.neighbor_counts
#
##        # Plot to look at nearest neighbor overlayed on position
##        x = list(list(zip(*lcPos))[0])
##        y = list(list(zip(*lcPos))[1])
##        z = list(list(zip(*lcPos))[2])
##        plt.scatter(x, y, c=list(neigh), s=0.5, edgecolor='none')
##        ax = plt.gca()
##        ax.axis('equal')
##        plt.show()
#
##        # This give me the average lattice spacing
##        distances = []
##        for (i, j) in nlist:
##            distances.append(np.linalg.norm(lcPos[i] - lcPos[j]))
##        lat = sum(distances) / len(distances)
##        print(lat)
#
#        # Let's throw out any point that has >3 neighbors
#        edPos = []
#        for k in range(0, len(neigh)):
#            if neigh[k] <= 3:
#                edPos.append(lcPos[k])
#
##        # Recluster this data
##        ed = freud.locality.AABBQuery(box=f_box, points=edPos)
##        edlist = ed.query(edPos, dict(num_neighbors=2, exclude_ii=True)).toNeighborList()
##        nlist.filter_r(r_max=1.5, r_min=0.2)
##        edneigh = nlist.neighbor_counts
##
##        # This give me the average lattice spacing
##        distances = []
#        x = list(list(zip(*edPos))[0])
#        y = list(list(zip(*edPos))[1])
#        z = list(list(zip(*edPos))[2])
##        plt.scatter(x, y, s=0.5, edgecolor='none')
##        for (i, j) in edlist:
##            distances.append(np.linalg.norm(edPos[i] - edPos[j]))
##            plt.plot([edPos[i][0], edPos[j][0]], [edPos[i][1], edPos[j][1]], lw=0.5)
##        print(len(x))
##        print(len(distances))
##        surface_distance = sum(distances)
##        print(surface_distance)
##        ax = plt.gca()
##        ax.axis('equal')
##        plt.show()
#
#        # Let's try this as a voronoi compute
#        voro = freud.locality.Voronoi()
#        voro.compute((f_box, f_box.wrap(edPos)))
#        plt.figure()
#        ax = plt.gca()
#        voro.plot(ax=ax)
#        ax.scatter(x, y, s=0.5)
#        plt.show()
