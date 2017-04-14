#!/usr/bin/env python

import os
from itertools import product
from create_centerlines import run

from pdb import set_trace

def latin(*paramsets):
    outparams = set()
    for i,paramset1 in enumerate(paramsets):
        np1 = len(paramset1)
        for params1 in paramset1:
            params = []
            for j,paramset2 in enumerate(paramsets):
                if j != i:
                    params.append(paramset2[0])
                else:
                    params.append(params1)
            outparams.add(tuple(params))
    return tuple(outparams)

def simpletest():
    input_shp = "test/output_roads_2.shp"
    output_driver = "GeoJSON"

    max_points = 3000 # 3000
    simplification = 0.05 # 0.05
    smooth = 0. # 5.
    segmentize_maxlen = 0.5 # 0.5
    morpho_dist = 0. # 0.

    param = (
        max_points,
        simplification,
        smooth,
        segmentize_maxlen,
        morpho_dist
    )

    paramstr = "_".join([str(p) for p in param])
    output_file = os.path.splitext(input_shp)[0] + "_" + paramstr + ".geojson"
    print output_file
    
    run(
        input_shp,
        output_file,
        segmentize_maxlen,
        max_points,
        simplification,
        smooth,
        morpho_dist,
        output_driver
    )

                
def fulltest():

    input_shp = "test/output_roads_2.shp"
    output_driver = "GeoJSON"

    max_points = (3000, 6000, 1500)
    simplification = (0.05, 0.1, 0.025)
    smooth = (5, 10, 2.5)
    segmentize_maxlen = (0.5, 1.0, 0.25)
    morpho_dist = (0.0, 10.)

    params = latin(
        max_points,
        simplification,
        smooth,
        segmentize_maxlen,
        morpho_dist
    )

    for param in params:
        print param

        (max_points,
         simplification,
         smooth,
         segmentize_maxlen,
         morpho_dist) = param

        paramstr = "_".join([str(p) for p in param])
        output_file = os.path.splitext(input_shp)[0] + outstr + ".geojson"
        print output_file
        
        run(
            input_shp,
            output_file,
            segmentize_maxlen,
            max_points,
            simplification,
            smooth,
            morpho_dist,
            output_driver
        )
         
if __name__ == "__main__":

    #fulltest()
    simpletest()
