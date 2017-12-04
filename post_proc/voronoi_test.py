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

# going to make a test situation of 5 points
part_num = 5
point1 = (0, 0)
point2 = (1, 0)
point3 = (0, 1)
point4 = (-1, 0)
point5 = (5, 5)

dumps = 1

points = np.zeros((part_num), dtype=np.ndarray)
points = [point1, point2, point3, point4, point5]
keepers = []

# Find close points...
neighbor_list = np.empty(part_num, dtype=np.object)
neighbor_list[:] = [], [], [], [], []
for iii in range(0, part_num):
    for jjj in range(0, part_num):
        if iii != jjj:
            if computeDistance(points[iii], points[jjj]) <= 10.0:
                neighbor_list[iii].append(jjj)

# Draw circles
for iii in range(0, part_num):
    for jjj in range(0, len(neighbor_list[iii])):
        for hhh in range(0, len(neighbor_list[iii])):
            if hhh > jjj:
                if lineCheck(points[iii], points[neighbor_list[iii][jjj]], points[neighbor_list[iii][hhh]]):
                    center, radius = computeCircle(points[iii],
                                                   points[neighbor_list[iii][jjj]],
                                                   points[neighbor_list[iii][hhh]])
                    keepers.append((center, radius))

# Check if any points are in a circle
iii = 0
while iii < len(keepers):
    for jjj in range(0, part_num):
        if computeDistance(keepers[iii][0], points[jjj]) < keepers[iii][1]:
            del keepers[iii]
            iii -= 1
            break
    iii += 1

# Plotting
import matplotlib.pyplot as plt

# get points in xs and ys for scatter
xs = np.zeros((part_num), dtype=np.float32)
ys = np.zeros((part_num), dtype=np.float32)
for iii in range(0, part_num):
    xs[iii] = points[iii][0]
    ys[iii] = points[iii][1]

# plot circles
fig, ax = plt.subplots()
for iii in range(0, len(keepers)):
    circle = plt.Circle(keepers[iii][0], keepers[iii][1], fill=False)
    ax.add_artist(circle)
# plot points
plt.scatter(xs, ys)

ax.set_aspect('equal')
ax.set_xlim(-6,6)
ax.set_ylim(-6,6)
plt.show()






