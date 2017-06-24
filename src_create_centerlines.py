#!/usr/bin/env python

# Author:  Joachim Ungar <joachim.ungar@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2015 EOX IT Services GmbH
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies of this Software or works derived from this Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#-------------------------------------------------------------------------------

from shapely.geometry import shape, Polygon, LineString,\
    MultiLineString, Point, MultiPoint, mapping
from shapely.wkt import loads
import ogr
from scipy.spatial import Voronoi
import networkx as nx
from itertools import combinations
import numpy as np
from scipy.ndimage import filters
from math import *

from pdb import set_trace

# number of longest paths sent to get_least_curved_path()
# TODO: someday this will make something break
NTOP = 10

MIN_AREA = 4.0

# this is for the geometry cleaning
TINY = 1.e-9


import fiona
def write_shapes(geoms, filename):
    print "writing", filename
    schema = {'geometry': geoms[0].geom_type, 'properties': {'id': 'int'}}    
    with fiona.open(filename, 'w', 'ESRI Shapefile', schema) as c:
        for i, geom in enumerate(geoms):
            c.write({'geometry': mapping(geom), 'properties': {'id': i}})


def get_centerlines_from_geom(
    geometry, feature_name, segmentize_maxlen, max_points, simplification,
    smooth_sigma, morpho_dist, minbranchlen, smoothed=False):
    """
    Returns centerline (for Polygon) 
    or centerlines (for MultiPolygons)
    as LineString or MultiLineString geometries.
    """
    if geometry.geom_type not in ["MultiPolygon", "Polygon"]:
        raise TypeError(
            "Geom type must be Polygon or MultiPolygon, not %s" %\
                geometry.geom_type)

    if geometry.area < MIN_AREA:
        return None
    
    # clean the geometry
    geometry = geometry.buffer(TINY)

    # the first num_inner interior rings are OK
    #interiors = geometry.interiors
    #num_inner = len(interiors)

    # dilate and erode the feature to smooth it
    if morpho_dist and not smoothed:
        print "smoothing", morpho_dist
        geometry = geometry.buffer(-morpho_dist)
        geometry = geometry.buffer(morpho_dist)
        smoothed = True

    
    if geometry.geom_type == "MultiPolygon":
        # recursion so that code below operates on Polygon objects
        print "recursing", len(geometry)
        centerline_geoms = []
        for subgeom in geometry:
            geom = get_centerlines_from_geom(
                subgeom, feature_name, segmentize_maxlen, max_points,
                simplification, smooth_sigma, morpho_dist, minbranchlen,
                smoothed)
            if geom is not None:
                if geom.geom_type == "LineString":
                    centerline_geoms.append(geom)
                else:
                    centerline_geoms.extend([g for g in geom])
        try:
            out_centerlines = MultiLineString(centerline_geoms)
        except TypeError, e:
            print e
            out_centerlines = LineString(centerline_geoms)
        except Exception, e:
            print e
        return out_centerlines

    else:
        """
        # remove interior created by the smoother
        num_inner2 = len(geometry.interiors)
        print "interior rings", num_inner, num_inner2
        if num_inner2 > num_inner:
            interiors = geometry.interiors
            nbad = num_inner2 - num_inner
            print "interior rings created by smoother", nbad
            for i in range(1,nbad+1):
                geometry = geometry.union(Polygon(interiors[-i]))
                print len(interiors), len(geometry.interiors)
        """

        # make cuts to interior rings
        if len(geometry.interiors) > 0:
            print "cutting interior rings", len(geometry.interiors)
            interiors = geometry.interiors
            exterior = geometry.exterior
            for interior in interiors:
                centroid = interior.centroid
                dproj = exterior.project(centroid)
                pproj = exterior.interpolate(dproj)
                line = LineString([centroid, pproj])
                thickline = line.buffer(0.0001)
                geometry = geometry.difference(thickline)

        if geometry.geom_type == "Polygon":
            print "final geometry.interiors", len(geometry.interiors)

        print "done cutting"
        
        if geometry.geom_type == "MultiPolygon":
            # recursion so that code below operates on Polygon objects
            print "recursing - cuts created multiple polygons", len(geometry)
            centerline_geoms = []
            for subgeom in geometry:
                print "subgeom.area", subgeom.area
                geom = get_centerlines_from_geom(
                    subgeom, feature_name, segmentize_maxlen, max_points,
                    simplification, smooth_sigma, morpho_dist,
                    minbranchlen, smoothed)
                if geom is not None:
                    if geom.geom_type == "LineString":
                        centerline_geoms.append(geom)
                    else:
                        centerline_geoms.extend([g for g in geom])
            try:
                out_centerlines = MultiLineString(centerline_geoms)
            except TypeError, e:
                print e
                out_centerlines = LineString(centerline_geoms)
            except Exception, e:
                print e
            return out_centerlines

        else:
            print "not multipolygon"
            assert len(geometry.interiors) == 0, len(geometry.interiors)

        # convert Polygon to Linestring
        boundary = geometry.boundary
        print "boundary.length", boundary.length
        
        # convert to OGR object and segmentize
        ogr_boundary = ogr.CreateGeometryFromWkb(boundary.wkb)
        ogr_boundary.Segmentize(segmentize_maxlen)
        segmentized = loads(ogr_boundary.ExportToWkt())

        # get points from the polygon
        points = segmentized.coords

        tolerance = simplification
        while len(points) > max_points:
            # if geometry is too large, apply simplification
            # until geometry is simplified enough
            # (indicated by the "max_points" value)
            tolerance += simplification
            simplified = boundary.simplify(tolerance)
            points = simplified.coords

        # calculate Voronoi diagram
        print "calculate Voronoi diagram"
        vor = Voronoi(points)

        # the next three steps are the most processing intensive and probably
        # not the most efficient method to get the skeleton centerline

        # convert to networkx graph
        print "convert to networkx graph"
        graph = graph_from_voronoi(vor, geometry)

        # get end nodes from graph
        print "get end nodes from graph"
        end_nodes = get_end_nodes(graph)

        if len(end_nodes) < 2:
            print "not enough nodes"
            return None

        # get longest path SLOW
        print "get longest path SLOW"
        paths_sorted, path_dists = get_longest_paths(end_nodes, graph)
        # TODO: maybe can change this back
        #paths_sorted = get_longest_paths(end_nodes, graph)

        # get least curved path out of the NTOP longest
        print "get least curved path out of the NTOP longest"
        longest_paths = paths_sorted[:NTOP]
        best_path = get_least_curved_path(longest_paths, vor.vertices)
        centerline = LineString(vor.vertices[best_path])

        print "for path in paths_sorted"
        
        for i,path in enumerate(paths_sorted):
            if path != best_path:
                # get the branch geometries
                line = LineString(vor.vertices[path])
                branches = line.difference(centerline)
                # branches might have multiple segments
                if branches.type != "MultiLineString":
                    branches = [branches]
                
                for branch in branches:
                    nodes = set(path).difference(best_path)
                    nnodes = len(nodes)
                    # attach the branch if it is long enough
                    if nnodes > 1 and branch.length > minbranchlen \
                       and centerline.distance(branch) == 0:
                        if branch.length <=0:
                            raise Exception, "this shouldn't happen"
                        centerline = centerline.union(line)
        
        # smooth out geometry
        if smooth_sigma > 0.:
            centerline_smoothed = smooth_linestring(centerline, smooth_sigma)
        else:
            centerline_smoothed = centerline

        return centerline_smoothed


def smooth_linestring(linestring, smooth_sigma):
    """
    Uses a gauss filter to smooth out the LineString coordinates.
    """
    smooth_x = np.array(filters.gaussian_filter1d(
        linestring.xy[0],
        smooth_sigma)
        )
    smooth_y = np.array(filters.gaussian_filter1d(
        linestring.xy[1],
        smooth_sigma)
        )
    smoothed_coords = np.hstack((smooth_x, smooth_y))
    smoothed_coords = zip(smooth_x, smooth_y)
    linestring_smoothed = LineString(smoothed_coords)
    return linestring_smoothed


def get_longest_paths(nodes, graph):
    """ returns longest path of all possible paths between a list of nodes """
    paths = []
    distances = []
    possible_paths = list(combinations(nodes, r=2))
    for node1, node2 in possible_paths:
        try:
            path = nx.shortest_path(graph, node1, node2, "weight")
        except Exception,e:
            path = []
        if len(path) > 1:
            distance = get_path_distance(path, graph)
            paths.append(path)
            distances.append(distance)

    szdp =  sorted(zip(distances, paths), reverse=True)
    path_dists, paths_sorted = zip(*szdp)
    return paths_sorted, path_dists


def get_least_curved_path(paths, vertices):
    angle_sums = []
    for path in paths:
        path_angles = get_path_angles(path, vertices)
        angle_sum = abs(sum(path_angles))
        angle_sums.append(angle_sum)
    paths_sorted = [x for (y,x) in sorted(zip(angle_sums, paths))]

    return paths_sorted[0]


def get_path_angles(path, vertices):
    angles = []
    prior_line = None
    next_line = None
    for index, point in enumerate(path):
        if index > 0 and index < len(path)-1:
            prior_point = vertices[path[index-1]]
            current_point = vertices[point]
            next_point = vertices[path[index+1]]
            angles.append(
                get_angle(
                    (prior_point, current_point), (current_point, next_point)
                )
            )

    return angles


def get_angle(line1, line2):
    v1 = line1[0] - line1[1]
    v2 = line2[0] - line2[1]
    angle = np.math.atan2(np.linalg.det([v1,v2]),np.dot(v1,v2))
    return np.degrees(angle)

def get_path_distance(path, graph):
    """ return weighted path distance """
    distance = 0
    for i, w in enumerate(path):
        j = i + 1
        if j < len(path):
            distance += round(graph.edge[path[i]][path[j]]['weight'], 6)
    return distance


def get_end_nodes(graph):
    """ return list of nodes with just one neighbor node """
    nodelist = [
        i for i in graph.nodes_iter() if len(graph.neighbors(i))==1
    ]
    return nodelist


def graph_from_voronoi(vor, geometry):
    """
    Creates a networkx graph out of all Voronoi ridge vertices which are inside
    the original geometry.
    """
    graph = nx.Graph()
    for i in vor.ridge_vertices:
        if i[0]>-1 and i[1]>-1:
            point1 = Point(vor.vertices[i][0])
            point2 = Point(vor.vertices[i][1])
            # Eliminate all points outside our geometry.
            if point1.within(geometry) and point2.within(geometry):
                dist = point1.distance(point2)
                graph.add_nodes_from([i[0], i[1]])
                graph.add_edge(i[0], i[1], weight=dist)
    return graph


def multilinestring_from_voronoi(vor, geometry):
    """
    Creates a MultiLineString geometry out of all Voronoi ridge vertices which
    are inside the original geometry.
    """
    linestrings = []
    for i in vor.ridge_vertices:
        if i[0]>-1 and i[1]>-1:
            point1 = Point(vor.vertices[i][0])
            point2 = Point(vor.vertices[i][1])
            # Eliminate all points outside our geometry.
            if point1.within(geometry) and point2.within(geometry):
                linestring = LineString([point1, point2])
                linestrings.append(linestring)
    multilinestring = MultiLineString(linestrings)
    return multilinestring
