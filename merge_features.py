#!/usr/bin/env python

import os
import fiona
from shapely.geometry import shape, mapping, box, Point
from shapely.ops import cascaded_union


from pdb import set_trace

STARTDIR = "/home/rbraswell/repo/label_centerlines/test"
INFILE = "roads_mixed_gtet050.shp"
OUTFILE = "roads_mixed_gtet050_merged.shp"

CHECK_PT = Point((286787.855, 107588.807))


if __name__ == "__main__":

    features = fiona.open(os.path.join(STARTDIR, INFILE))
    schema = features.schema
    crs = features.crs
    driver = features.driver
    #out_schema = features.schema.copy()

    geoms = []
    
    for feature in features:

        fid = feature['id']
        geom = shape(feature['geometry'])

        assert geom.is_valid, "geom not valid"
        assert geom.geometryType() == "Polygon", "geom not polygon"

        assert geom.contains(CHECK_PT) is False, "oops"
        
        #geom = geom.buffer(1.e-9)
        geom = geom.buffer(0.1)
        geom = geom.buffer(-0.1)

        area = geom.area
        length = geom.length
        n_inner = len(geom.interiors)
        bbox = box(*geom.bounds)
        bbarea = bbox.area

        print int(fid), n_inner, area, length, length/area, bbarea

        geoms.append(geom)

    geom_all = cascaded_union(geoms)
    assert geom_all.contains(CHECK_PT) is False, "crap"

    #set_trace()

    out_schema = {'geometry': 'MultiPolygon', 'properties': {'id': 'int'}}

    outfile = os.path.join(STARTDIR, OUTFILE)

    out_features = fiona.open(outfile, 'w', schema=out_schema,
                              crs=crs, driver=driver)

    out_feature = {}
    out_feature['geometry'] = mapping(geom_all)
    out_feature['type'] = 'Feature'
    out_feature['id'] = '0'
    out_feature['properties'] = {'id': 0}
    
    out_features.write(out_feature)

    out_features.close()

