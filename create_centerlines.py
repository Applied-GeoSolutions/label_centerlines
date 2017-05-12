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

import os
import sys
import argparse
import fiona
import multiprocessing
from shapely.geometry import shape, mapping
from functools import partial

from src_create_centerlines import get_centerlines_from_geom

from pdb import set_trace

NUMPROC = 1
# TODO: this default only works with projected vector data sets
MINBRANCHLEN = 30

def worker(
    segmentize_maxlen,
    max_points,
    simplification,
    smooth_sigma,
    morpho_dist,
    minbranchlen,
    feature
    ):

    geom = shape(feature['geometry'])

    for name_field in ["name", "Name", "NAME"]:
        if name_field in feature["properties"]:
            feature_name = feature["properties"][name_field]
            break
        else:
            feature_name = None
    print "processing", feature_name

    #centerlines_geom = get_centerlines_from_geom(
    #    geom,
    #    segmentize_maxlen=segmentize_maxlen,
    #    max_points=max_points,
    #    simplification=simplification,
    #    smooth_sigma=smooth_sigma,
    #    morpho_dist=morpho_dist
    #)

    try:
        centerlines_geom = get_centerlines_from_geom(
            geom,
            segmentize_maxlen=segmentize_maxlen,
            max_points=max_points,
            simplification=simplification,
            smooth_sigma=smooth_sigma,
            morpho_dist=morpho_dist,
            minbranchlen=minbranchlen
            )
    except TypeError as e:
        print e
    except:
        raise

    if centerlines_geom:
        return (
            feature_name,
            {
                'properties': feature['properties'],
                'geometry': mapping(centerlines_geom)
            }
        )
    else:
        return (None, None)

def run(
    input_shp,
    output_file,
    segmentize_maxlen,
    max_points,
    simplification,
    smooth_sigma,
    morpho_dist,
    driver,
    numproc,
    minbranchlen
    ):

    extensions = {'ESRI Shapefile': '.shp', 'GeoJSON': '.geojson'}
    if os.path.splitext(output_file)[1] != extensions[driver]:
            output_file += extensions[driver]

    with fiona.open(input_shp, "r") as inp_polygons:
        out_schema = inp_polygons.schema.copy()
        out_schema['geometry'] = "LineString"
        if os.path.exists(output_file):
            os.remove(output_file)
        with fiona.open(
            output_file,
            "w",
            schema=out_schema,
            crs=inp_polygons.crs,
            driver=driver
            ) as out_centerlines:

            pool = multiprocessing.Pool(processes=numproc)
            func = partial(
                worker,
                segmentize_maxlen,
                max_points,
                simplification,
                smooth_sigma,
                morpho_dist,
                minbranchlen
            )

            try:
                feature_count = 0
                for feature_name, output in pool.imap_unordered(
                    func,
                    inp_polygons
                    ):
                    feature_count += 1
                    if output:
                        out_centerlines.write(output)
                        print "written feature %s: %s" %(
                            feature_count,
                            feature_name
                            )
                    else:
                        print "Invalid output for feature", feature_name
            except KeyboardInterrupt:
                print "Caught KeyboardInterrupt, terminating workers"
                pool.terminate()
            except Exception as e:
                if feature_name:
                    print ("%s: FAILED (%s)" %(feature_name, e))
                else:
                    print ("feature: FAILED (%s)" %(e))
                raise
            finally:
                pool.close()
                pool.join()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_shp",
        type=str,
        help="input polygons"
        )
    parser.add_argument(
        "output_file",
        type=str,
        help="output centerlines"
        )
    parser.add_argument(
        "--segmentize_maxlen",
        type=float,
        help="maximum length used when segmentizing polygon borders",
        default=0.5
        )
    parser.add_argument(
        "--max_points",
        type=int,
        help="number of points per geometry allowed before simplifying",
        default=3000
        )
    parser.add_argument(
        "--simplification",
        type=float,
        help="value which increases simplification when necessary",
        default=0.05
        )
    parser.add_argument(
        "--smooth",
        type=int,
        help="smoothness of the output centerlines",
        default=0
        )
    parser.add_argument(
        "--morpho_dist",
        type=float,
        help="distance for erosion and dilation",
        default=0.0
        )
    parser.add_argument(
        "--output_driver",
        type=str,
        help="write to 'ESRI Shapefile' or 'GeoJSON' (default)",
        default="ESRI Shapefile"
        #default="GeoJSON"
        )
    parser.add_argument(
        "--numproc",
        type=int,
        help="number of processors to use",
        default=NUMPROC
        )
    parsed = parser.parse_args(sys.argv[1:])

    run(
        parsed.input_shp,
        parsed.output_file,
        parsed.segmentize_maxlen,
        parsed.max_points,
        parsed.simplification,
        parsed.smooth,
        parsed.morpho_dist,
        parsed.output_driver,
        parsed.numproc,
        parsed.minbranchlen
    )
