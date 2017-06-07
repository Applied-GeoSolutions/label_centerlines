#!/usr/bin/env python

import os
from itertools import product
from create_centerlines import run

from pdb import set_trace

START_DIR = "/projects/cmsindo/data/logging/automated_accuracy/test_postprocess"
timberdana_roads = "005/roads_mixed_gtet050_dissolve_multiparts_idsep.shp"
timberdana_trails = "005/skidtrails_mixed_lt050_dissolve_multiparts_idsep.shp"
rodamas_roads = "008/roads_mixed_gtet050_dissolve_multiparts_idsep.shp"
rodamas_trails = "008/skidtrails_mixed_lt050_dissolve_multiparts_idsep.shp"


def makelines(input_shp):

    output_driver = "GeoJSON"

    max_points = 3000 # 3000
    simplification = 0.05 # 0.05
    smooth = 0. # 5.
    segmentize_maxlen = 0.5 # 0.5
    morpho_dist = 1. # 0.
    numproc = 5 # 5
    minbranchlen = 30 # 30

    param = (
        max_points,
        simplification,
        smooth,
        segmentize_maxlen,
        morpho_dist,
        minbranchlen
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
        output_driver,
        numproc,
        minbranchlen
    )


if __name__ == "__main__":

    input_shp = os.path.join(START_DIR, rodamas_trails)
    makelines(input_shp)
