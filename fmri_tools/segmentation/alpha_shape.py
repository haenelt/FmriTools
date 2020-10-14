# -*- coding: utf-8 -*-

# python standard library inputs
import math

# external inputs
import numpy as np
import shapely.geometry as geometry
from scipy.spatial import Delaunay
from shapely.ops import cascaded_union, polygonize


def alpha_shape(points, alpha):    
    """ Alpha shape

    This function computes the alpha shape (concave hull) of a set of points. It 
    is taken and slightly changed from [1]. Delaunay triangles are computed 
    which establish a connection between each point and nearby points and then 
    we remove some of the triangles that are too far from their neighbors. This 
    removal part is the key. By identifying candidates for removal we are saying 
    that these points are too far from their connected points so don't use that 
    connection as part of the boundary.    

    Parameters
    ----------
    points : list
        Iterable container of points.
    alpha : float
        Alpha value to influence the gooeyness of the border. Smaller numbers 
        don't fall inward as much as larger numbers. Too large, and you lose 
        everything!

    Returns
    -------
    cascaded_union(triangles) : poly
        Concave hull coordinates.
    edge_points : list
        Start and end points of all edges.

    References
    -------
    .. [1] http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/ 
    (accessed 22-10-2018)

    Notes
    -------
    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 12-10-2020

    """
    
    # when you have one triangle, there is no sense in computing an alpha shape
    if len(points) < 4:
        return geometry.MultiPoint(list(points)).convex_hull
    
    # add a line between point i and j if not alread added 
    def add_edge(edges, edge_points, coords, i, j):
        if (i, j) in edges or (j, i) in edges:
            return 
        edges.add( (i, j) )
        edge_points.append(coords[ [i, j] ])

    # loop over triangles (ia, ib, ic are indices of triangle corner points)
    coords = np.array(points)
    tri = Delaunay(coords)
    edges = set()
    edge_points = []
    for ia, ib, ic in tri.vertices:
        
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]
        
        # lengths of sides of triangle
        a = math.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
        b = math.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
        c = math.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)
        
        # semiperimeter of triangle
        s = (a + b + c)/2.0
        
        # area of triangle by Heron's formula
        area = math.sqrt(s*(s-a)*(s-b)*(s-c))
        circum_r = a*b*c/(4.0*area)
        
        # radius filter.
        if circum_r < 1.0/alpha:
            add_edge(edges, edge_points, coords, ia, ib)
            add_edge(edges, edge_points, coords, ib, ic)
            add_edge(edges, edge_points, coords, ic, ia)
            
    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    
    return cascaded_union(triangles), edge_points
