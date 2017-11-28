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


