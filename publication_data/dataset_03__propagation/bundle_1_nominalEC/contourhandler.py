"""
This module administrates the contours
8 June 2018.
"""

import numpy as np
from collections import OrderedDict
import csv




def load_contours(name):
    """ Open a cotours csv file and read it """
    contours = OrderedDict()
    # Load the contours
    with open(name, "r") as f:
        fr = csv.reader(f)
        frl = list(fr)
        for line in frl:
            try:
                values = [float(item) for item in line]
            except ValueError:
                key = line[0]
                contours[key] = []
            else:
                contours[key].append(values)
    return contours

def read_poly(file_name):
    """
    Simple poly-file reader, that creates a python dictionary 
    with information about vertices, edges and holes.
    It assumes that vertices have no attributes or boundary markers.
    It assumes that edges have no boundary markers.
    No regional attributes or area constraints are parsed.
    Taken and adapted from: ...
    """

    output = {'vertices': None, 'holes': None, 'segments': None}
    
    # open file and store lines in a list
    file = open(file_name, 'r')
    lines = file.readlines()
    file.close()
    lines = [x.strip('\n').split() for x in lines]
    
    # Store vertices
    vertices= []
    N_vertices, dimension, attr, bdry_markers = [int(x) for x in lines[0]]
    # We assume attr = bdrt_markers = 0
    for k in range(N_vertices):
        label, x, y = [items for items in lines[k+1]]
        vertices.append([float(x), float(y)])
    output['vertices']=np.array(vertices)
    
    # Store segments
    segments = []
    N_segments, bdry_markers = [int(x) for x in lines[N_vertices+1]]
    for k in range(N_segments):
        label, pointer_1, pointer_2 = [items for items in lines[N_vertices+k+2]]
        segments.append([int(pointer_1)-1, int(pointer_2)-1])
    output['segments'] = np.array(segments)
    
    # Store holes
    N_holes = int(lines[N_segments+N_vertices+2][0])
    holes = []
    for k in range(N_holes):
        label, x, y = [items for items in lines[N_segments + N_vertices + 3 + k]]
        holes.append([float(x), float(y)])
    
    output['holes'] = np.array(holes)
    
    return output

def move_element(dic, key, newpos):
    """ Move the element with key 'key' of a dictionary 'dic' to the 
    position 'newpos' """
    newdic = OrderedDict()
    for i, (k, v) in enumerate(dic.items()):
        if k != key:
            newdic[k] = v
        if i == newpos:
            newdic[key] = dic[key]
    return newdic
    
def fill_random(c, n):
    """ Fill a contour given by the array c with a maximum of n random 
    points """

    # Contour's enclosing box
    left, right = c[:, 0].min(), c[:, 0].max()
    bottom, top = c[:, 1].min(), c[:, 1].max()

    # Dimensions
    width = right - left
    height = top - bottom

    # Generate random points in the box
    x = np.random.uniform(low=left, high=right, size=n)
    y = np.random.uniform(low=bottom, high=top, size=n)
    points = np.array([x, y]).T

    # 'Planar' polygon object for the contour
    polygon = pl.Polygon(c.tolist())

    # Check if the points are inside the polygon
    good_points = []
    for p in points:
        if polygon.contains_point(p):
            good_points.append(p)

    return polygon, np.array(good_points)

def remove_points(points):
    """ Remove points from 'points' that fall inside fascicles """
    survivors = []
    # Iterate over input points
    for p in points:
        discard = False
        # Iterate over fascicles
        for k, fp in polygons.items():
            if "Fascicle" in k:
                # Discard it and end this sub-loop if the point is inside the
                # fascicle's polygon
                if fp.contains_point(p):
                    discard = True
                    break
        # If it survived the loop, save it
        if not discard:
            survivors.append(p)

    return np.array(survivors)

def reduce_points(c, n):
    """ Reduce the number of points of contours 'c' by a factor of n in 
    order to reduce the computational burden. Only for contours with 
    more than n points """
    for k, a in c.items():
        if len(a) > n:
            c[k] = a[::n]
    return c

def c2pslg(c):
    """ Convert the contours of a nerve stored in the dictionary 'c' to a 
    dictionary describing a Planar Straight Line Graph (PSLG) """

    # newc = {'vertices': None, 'holes': None, 'segments': None}
    newc = {}
    # Total number of points
    nt = np.array([len(a) for a in c.values()]).sum()
    # Counter of processed points
    npr = 0
    # Generate the information
    vertices = []
    holes = []
    segments = []

    # Iterate over all the contours
    for k, a in c.items():
        n = len(a)
        i = 0
        # Append vertices and segments
        for (x, y) in a:
            vertices.append([x, y])
            # Only if there's more than one point here, add segments
            if n > 1:
                if i < n-1:
                    segments.append([npr+i, npr+i+1])
                    i += 1
                else:
                    segments.append([npr+i, npr])
        # Update npr
        npr += n

        # Holes
        if "Fascicle" in k:
            # Add the centroid (hand-maded for speed-up)
            a = np.array(a)
            holes.append([a[:, 0].mean(), a[:, 1].mean()])

    newc['vertices'] = np.array(vertices)
    # Condition: for there to be a 'holes' entry, there must be holes
    if len(holes) > 0:
        newc['holes'] = np.array(holes)
    newc['segments'] = np.array(segments)
    return newc