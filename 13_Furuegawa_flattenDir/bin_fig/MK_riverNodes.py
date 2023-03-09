
def reproj_Points(x_list, y_list, src_epsg, dst_epsg):
    """
    this function is reprojection of a point from src_epsg to dst_espg
    """
    import geopandas as gpd
    from shapely.geometry import Point
    x_src = x_list
    y_src = y_list
    gdf = gpd.GeoDataFrame( geometry=gpd.points_from_xy(x=x_src, y=y_src ) ).set_crs(src_epsg)
    gdf_pl = gdf.to_crs(dst_epsg)
    x_dst = []
    y_dst = []
    for p in list(gdf_pl.geometry):
        xx, yy = p.xy
        x_dst.append(xx[0])
        y_dst.append(yy[0])
    return x_dst, y_dst

def make_GeoJSON(x, y, id, outJsonName="river_node_lonlat.geojson"):
    from geojson import Point, Feature, FeatureCollection, dump
    # https://qiita.com/XPT60/items/80827565dcf0f4ef3e91
    ft_all = []
    for x_i, y_i, id_i in zip(x,y,id):
        ft = Feature(geometry=Point((x_i,y_i)),
                properties = {'name': "id={}".format(id_i),
#                    'marker_symbol': ""
                }
            )
        ft_all.append(ft)
    ft_collection = FeatureCollection(ft_all)
    with open(outJsonName, "w") as dst:
        dump(ft_collection, dst, indent=2)


def make_KmlMarkerSymbolFrmGoogle(datatype="streamConfigration"):
    if datatype == "streamConfigration":
        marker_url = "https://labs.google.com/ridefinder/images/mm_20_blue.png"
    elif datatype == "CVcenter":
        marker_url = "https://labs.google.com/ridefinder/images/mm_20_blue.png"
    elif datatype == "CVflux_x":
        marker_url = "https://maps.google.com/mapfiles/dir_90.png"
    elif datatype == "CVflux_y":
        marker_url = "https://maps.google.com/mapfiles/dir_0.png"
    else:
        marker_url = ""
    return marker_url
    

def make_Kml(x,y,id, marker_url="", kmlname="river_node_lonlat.kml"):
    import simplekml
    kml = simplekml.Kml()
    for x_i, y_i, id_i in zip(x,y,id):
        # namestr = 'id: ' + str(point[2]) + ' ; ' + 'to_id: ' + str(point[3])
        pnt = kml.newpoint(name="id={}".format(id_i), coords=[(x_i, y_i)])
        pnt.style.iconstyle.icon.href = marker_url
    kml.save(kmlname)
    # return True


def make_RiverNodes(inputfile="streamConfiguration_inRR.txt", srcEPSG="EPSG:6670", tgtEPSG="EPSG:6668", format="kml"):
    import os
    import pandas as pd
    import geopandas as gpd
    tdir = os.getcwd()
    stream_inRR = pd.read_csv(os.path.join(tdir, inputfile),skiprows=1)
    id = [i for i in stream_inRR["id"]]
    x  = [x for x in stream_inRR["x_2d"]]
    y  = [y for y in stream_inRR["y_2d"]]
    lon, lat = reproj_Points(x,y, srcEPSG, tgtEPSG)
    print("Processing:  {}".format(inputfile))
    print("    num. of river nodes = {}".format(len(x)))
    outfilename = "river_node"
    if format == "GeoJSON":
        print("..  get, ", outfilename+"_xy.GeoJSON", outfilename+"_lonlat.GeoJSON")
        make_GeoJSON(x,y,id, outfilename+"_xy.GeoJSON")
        make_GeoJSON(lon,lat, id, outfilename+"_lonlat.GeoJSON")
    elif format == "kml":
        print("..  get, ", outfilename+"_xy.kml", outfilename+"_lonlat.kml")
        kml_marker_url = make_KmlMarkerSymbolFrmGoogle(datatype="streamConfigration")
        make_Kml(x, y, id, marker_url=kml_marker_url, kmlname=outfilename+"_xy.kml")
        make_Kml(lon,lat,id, marker_url=kml_marker_url, kmlname=outfilename+"_lonlat.kml")
    else:
        pass


def make_CVcenter(inputfile="controlVolume_centerPoint_inRR.txt", srcEPSG="EPSG:6670", tgtEPSG="EPSG:6668", format="kml"):
    import os
    import pandas as pd
    import geopandas as gpd
    tdir = os.getcwd()
    data = pd.read_csv(os.path.join(tdir, inputfile),skiprows=1)
    x = [x for x in data["x"]]
    y = [y for y in data["y"]]
    id = [int(i)+1 for i in range(len(x))]
    lon, lat = reproj_Points(x,y, srcEPSG, tgtEPSG)
    print("Processing:  {}".format(inputfile))
    print("    num. of center point = {}".format(len(x)))
    outfilename = os.path.splitext(inputfile)[0]
    if format == "GeoJSON":
        print("..  get, ", outfilename+"_xy.GeoJSON", outfilename+"_lonlat.GeoJSON")
        make_GeoJSON(x,y,id, outfilename+"_xy.GeoJSON")
        make_GeoJSON(lon,lat, id, outfilename+"_lonlat.GeoJSON")
    elif format == "kml":
        print("..  get, ", outfilename+"_xy.kml", outfilename+"_lonlat.kml")
        kml_marker_url = make_KmlMarkerSymbolFrmGoogle(datatype="CVcenter")
        make_Kml(x, y, id, marker_url=kml_marker_url, kmlname=outfilename+"_xy.kml")
        make_Kml(lon,lat,id, marker_url=kml_marker_url, kmlname=outfilename+"_lonlat.kml")   
    else:
        pass    

def make_CVflux_x(inputfile="fluxPoint_xBoudary_inRR.txt", srcEPSG="EPSG:6670", tgtEPSG="EPSG:6668", format="kml"):
    import os
    import pandas as pd
    import geopandas as gpd
    tdir = os.getcwd()
    data = pd.read_csv(os.path.join(tdir, inputfile),skiprows=1)
    x = [x for x in data["x"]]
    y = [y for y in data["y"]]
    id = [int(i)+1 for i in range(len(x))]
    lon, lat = reproj_Points(x,y, srcEPSG, tgtEPSG)
    print("Processing:  {}".format(inputfile))
    print("    num. of point on edge in x dir. = {}".format(len(x)))
    outfilename = os.path.splitext(inputfile)[0]
    if format == "GeoJSON":
        print("..  get, ", outfilename+"_xy.GeoJSON", outfilename+"_lonlat.GeoJSON")
        make_GeoJSON(x,y,id, outfilename+"_xy.GeoJSON")
        make_GeoJSON(lon,lat, id, outfilename+"_lonlat.GeoJSON")
    elif format == "kml":
        print("..  get, ", outfilename+"_xy.kml", outfilename+"_lonlat.kml")
        kml_marker_url = make_KmlMarkerSymbolFrmGoogle(datatype="CVflux_x")
        make_Kml(x, y, id, marker_url=kml_marker_url, kmlname=outfilename+"_xy.kml")
        make_Kml(lon,lat,id, marker_url=kml_marker_url, kmlname=outfilename+"_lonlat.kml")   
    else:
        pass

def make_CVflux_y(inputfile="fluxPoint_yBoudary_inRR.txt", srcEPSG="EPSG:6670", tgtEPSG="EPSG:6668", format="kml"):
    import os
    import pandas as pd
    import geopandas as gpd
    tdir = os.getcwd()
    data = pd.read_csv(os.path.join(tdir, inputfile),skiprows=1)
    x = [x for x in data["x"]]
    y = [y for y in data["y"]]
    id = [int(i)+1 for i in range(len(x))]
    lon, lat = reproj_Points(x,y, srcEPSG, tgtEPSG)
    print("Processing:  {}".format(inputfile))
    print("    num. of point on edge in y dir = {}".format(len(x)))
    outfilename = os.path.splitext(inputfile)[0]
    if format == "GeoJSON":
        print("..  get, ", outfilename+"_xy.GeoJSON", outfilename+"_lonlat.GeoJSON")
        make_GeoJSON(x,y,id, outfilename+"_xy.GeoJSON")
        make_GeoJSON(lon,lat, id, outfilename+"_lonlat.GeoJSON")
    elif format == "kml":
        print("..  get, ", outfilename+"_xy.kml", outfilename+"_lonlat.kml")
        kml_marker_url = make_KmlMarkerSymbolFrmGoogle(datatype="CVflux_y")
        make_Kml(x, y, id, marker_url=kml_marker_url, kmlname=outfilename+"_xy.kml")
        make_Kml(lon,lat,id, marker_url=kml_marker_url, kmlname=outfilename+"_lonlat.kml")   
    else:
        pass

import argparse
parser = argparse.ArgumentParser(\
    description="""outputs nodes as kml and GeoJSON formats in two coordinate system.
        1st argument is given the EPSG code of the source x-y coodinate system, which is the same system as the DFSS calculation.
        2nd argument is given the EPSG code of the target coordinate system. It is assumed to use lon-lat system.
        Three kind of output files can be chosen: 1) river nodes (default), 2) boundary nodes on Control Volume, 
        3) center nodes in Control Volume.

    """
    )
parser.add_argument("sourceEPGScode", \
                    type=str, \
                    # default="", \
                    help = "EPSG code of x-y coordinate system of DFSS'S output files; ex) JGD2011 CS2 ->  'EPSG:6670'",
                    )
parser.add_argument("targetEPGScode", \
                    type=str, \
                    default="EPSG:6668", \
                    help = "EPSG code of lon-lat coordinate system; ex) JGD2011 lonlat ->  'EPSG:6668'",
                    )
parser.add_argument("-k", "--kind",\
                    choices = ["river", "edge", "center", "full"],
                    default="river",
                    help = "to choose kind of output files. except river, it takes a lot of calculation time",
                    )
parser.add_argument("-f", "--format",\
                    choices = ["kml", "GeoJSON"],
                    default="kml",
                    help = "to choose output format: kml [default] or GeoJSON",
                    )
parser.add_argument("--version", version="%(prog)s 1.0.0",
                    action="version",
                    default=False)
args = parser.parse_args()

try:
    srcEPSG = int(args.sourceEPGScode.split(":")[1])
    tgtEPSG = int(args.targetEPGScode.split(":")[1])
except Exception as e:
    import sys, traceback
    print(traceback.format_exc())
    sys.exit(-1)

if args.kind == "river":
    make_RiverNodes(srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
elif args.kind == "edge":
    make_CVflux_x(inputfile="fluxPoint_x_inRR.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVflux_x(inputfile="fluxPoint_x_inDF.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVflux_x(inputfile="fluxPoint_xBoudary_inRR.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVflux_x(inputfile="fluxPoint_xBoudary_inDF.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVflux_y(inputfile="fluxPoint_y_inRR.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVflux_y(inputfile="fluxPoint_y_inDF.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVflux_y(inputfile="fluxPoint_yBoudary_inRR.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVflux_y(inputfile="fluxPoint_yBoudary_inDF.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)                
elif args.kind == "center":
    make_CVcenter(inputfile="controlVolume_centerPoint_inRR.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVcenter(inputfile="controlVolume_centerPoint_inDF.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
elif args.kind == "full":
    make_RiverNodes(srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVflux_x(inputfile="fluxPoint_x_inRR.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVflux_x(inputfile="fluxPoint_x_inDF.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVflux_x(inputfile="fluxPoint_xBoudary_inRR.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVflux_x(inputfile="fluxPoint_xBoudary_inDF.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVflux_y(inputfile="fluxPoint_y_inRR.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVflux_y(inputfile="fluxPoint_y_inDF.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVflux_y(inputfile="fluxPoint_yBoudary_inRR.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVflux_y(inputfile="fluxPoint_yBoudary_inDF.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)                
    make_CVcenter(inputfile="controlVolume_centerPoint_inRR.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
    make_CVcenter(inputfile="controlVolume_centerPoint_inDF.txt", srcEPSG=srcEPSG, tgtEPSG=tgtEPSG, format=args.format)
else:
    import sys
    parser.print_help(sys.stderr)
    sys.exit(1)

