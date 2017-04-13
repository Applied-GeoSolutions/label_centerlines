#!/usr/bin/env python

from itertools import product
from create_centerlines import run

def main():

    input_shp = "test/output_roads_1.shp"
    output_driver = "GeoJSON"

    max_points = [3000, 6000, 1500]
    simplification = [0.05, 0.1, 0.025]
    smooth = [5, 10, 2.5]
    segmentize_maxlen = [0.5, 1.0, 0.25]
    morpho_dist = [0.0, 10.]
    
    params = product(
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

        outstr = "_".join([str(p) for p in param])

        output_file = "test/out_" + outstr + ".geojson"

        print output_file
        
        run(
            input_shp,
            output_file,
            segmentize_maxlen,
            max_points,
            simplification,
            smoothg,
            morpho_dist,
            output_driver
        )
         
if __name__ == "__main__":
    main()
