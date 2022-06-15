
# coding: utf-8

# In[1]:


# Script to take polygons of training data and convert to points by location of pixel centroid per a single naip tile
import os
import h5py
import gdal
import pandas as pd
import geopandas as gpd
import numpy as np

from collections import OrderedDict
import numpy
import fiona
import geopandas as gpd
from shapely.geometry import shape, Point, mapping
from datetime import datetime


# In[11]:

hdf5_dir = os.path.abspath("Z:/SRER/Martha/RFdata/plots_hyperspec")
in_polygon = os.path.abspath("Z:/SRER/Martha/RFdata/3mbuff.shp")
out_dir = os.path.abspath("Z:/SRER/Martha/RFdata")

h5_files = []
for root,dirs,files in os.walk(hdf5_dir):
    for file in files:
        if file.endswith(".h5"):
            fpath = os.path.join(root,file)
            h5_files.append(fpath)


# In[3]:


def getRasterData(rasterloc):
    #Read in reflectance hdf5 file 
    hdf5_file = h5py.File(rasterloc,"r")

    #Get the site name
    file_attrs_string = str(list(hdf5_file.items()))
    file_attrs_string_split = file_attrs_string.split("\"")
    sitename = file_attrs_string_split[1]
    
    #Extract the reflectance & wavelength datasets
    refl = hdf5_file[sitename]["Reflectance"]
    reflData = refl["Reflectance_Data"]
    
    wavelengths = refl["Metadata"]["Spectral_Data"]["Wavelength"].value
 
    hdf5_file.close        
    
    return wavelengths, reflData


# In[4]:


def getRasterData(rasterloc):
    #Read in reflectance hdf5 file 
    hdf5_file = h5py.File(rasterloc,"r")

    #Get the site name
    file_attrs_string = str(list(hdf5_file.items()))
    file_attrs_string_split = file_attrs_string.split("\"")
    sitename = file_attrs_string_split[1]
    
    #Extract the reflectance & wavelength datasets
    refl = hdf5_file[sitename]["Reflectance"]
    reflData = refl["Reflectance_Data"]
    
    wavelengths = refl["Metadata"]["Spectral_Data"]["Wavelength"].value
 
    hdf5_file.close        
    
    return wavelengths, reflData


# In[5]:


def getRasterInfo(rasterloc):
    #Read in reflectance hdf5 file 
    try:
        hdf5_file = h5py.File(rasterloc,"r")
    except:
        raise Exception(f"h5py file open failure for {rasterloc}")

    #Get the site name
    file_attrs_string = str(list(hdf5_file.items()))
    file_attrs_string_split = file_attrs_string.split("\"")
    sitename = file_attrs_string_split[1]
    
    #Extract the reflectance & wavelength datasets
    refl = hdf5_file[sitename]["Reflectance"]
    reflData = refl["Reflectance_Data"]
    #reflRaw = refl["Reflectance_Data"].value
    
    #Create dictionary containing relevant metadata information
    metadata = {}
    #Map_Info sample : "UTM,  1.000,  1.000,  505834.000,  3515363.000,  1.0000000000e+000,  1.0000000000e+000,  12,  North,  WGS-84,  units=Meters, 0"
    map_info = str(refl["Metadata"]["Coordinate_System"]["Map_Info"].value).split(",")
    metadata["resx"] = float(map_info[5])
    metadata["resy"] = float(map_info[6])
    #metadata["wavelength"] = refl["Metadata"]["Spectral_Data"]["Wavelength"].value

    #Extract no data value & scale factor
    metadata["data ignore value"] = float(reflData.attrs["Data_Ignore_Value"])
    metadata["reflectance scale factor"] = float(reflData.attrs["Scale_Factor"])
    #metadata["interleave"] = reflData.attrs["Interleave"]

    #Extract spatial extent from attributes
    #    [  516087.   516921.  3516971.  3532372.]
    # TODO recalcualte extent so it"s units independent
    spatial_extent = reflData.attrs["Spatial_Extent_meters"]
    bbox = {}
    bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"] = map(float, spatial_extent )
    metadata["extent"] = bbox

    #Extract bad band windows
    metadata["bad band window1"] = (refl.attrs["Band_Window_1_Nanometers"])
    metadata["bad band window2"] = (refl.attrs["Band_Window_2_Nanometers"])

    #Extract projection information
    metadata["projection"] = refl["Metadata"]["Coordinate_System"]["Proj4"].value
    metadata["epsg"] = int(refl["Metadata"]["Coordinate_System"]["EPSG Code"].value)

    hdf5_file.close        
    
    return metadata


# In[6]:


def getSnappedLoc(x, y, raster_extent, px_x_size, px_y_size):
    
    diffx = abs(x - raster_extent["xmin"])
    diffy = abs(raster_extent["ymax"]- y)

    number_of_pixels_to_x = int(diffx / px_x_size)
    number_of_pixels_to_y = int(diffy / px_y_size)

    pixel_xdiff = float("{0:.11f}".format( diffx % px_x_size ))  # get modulo pixel difference to float precision of 11 decimals
    pixel_ydiff = float("{0:.11f}".format( diffy % px_y_size ))  # get modulo pixel difference to float precision of 11 decimals

    #snapped pixel locations
    if abs(pixel_xdiff) > abs(px_x_size / 2):
        snapped_x = number_of_pixels_to_x + 1
    else:
        snapped_x = number_of_pixels_to_x

    if abs(pixel_ydiff) > abs(px_y_size / 2):
        snapped_y = number_of_pixels_to_y + 1
    else:
        snapped_y = number_of_pixels_to_y

    if snapped_x % px_x_size != raster_extent["xmin"] % px_x_size:
        print(snapped_x % pix_xsize)
        raise ValueError("BAD PIXEL VALUE FOR ULX - ", snapped_ulx)
    if snapped_y % px_y_size != raster_extent["ymax"] % px_y_size:
        print(snapped_y % pix_ysize)
        raise ValueError("BAD PIXEL VALUE FOR ULY - ", snapped_uly)

    return snapped_x, snapped_y

def getSnappedBbox(vector_bbox, raster_metadata):
    ulx = vector_bbox["xmin"]
    uly = vector_bbox["ymax"]
    lrx = vector_bbox["xmax"]
    lry = vector_bbox["ymin"]
    
    #print("ULX: ", ulx, "ULY: ", uly)
    #print("LRX: ", lrx, "LRY: ", lry)
    """ Returns set of upper-uly snapped pixel locations in set as (x, y)"""
    pix_xsize = raster_metadata["resx"]  # ras_aff.a
    pix_ysize = raster_metadata["resy"]  # ras_aff.e

    raster_xmin = raster_metadata["extent"]["xmin"]
    raster_ymax = raster_metadata["extent"]["ymax"]
    
    #print("STARTING UPPER LEFT...")
    ulx, uly = getSnappedLoc(ulx, uly, raster_metadata["extent"], pix_xsize, pix_ysize)
    #print("STARTING LOWER RIGHTLEFT...")
    lrx, lry = getSnappedLoc(lrx, lry, raster_metadata["extent"], pix_xsize, pix_ysize)

    raster_xmin = raster_metadata["extent"]["xmin"]
    raster_ymax = raster_metadata["extent"]["ymax"]

    crs_xmin = (ulx * pix_xsize) + raster_xmin
    crs_xmax = (lrx * pix_xsize) + raster_xmin
    crs_ymin = raster_ymax - (lry * pix_ysize)
    crs_ymax = raster_ymax - (uly * pix_ysize)

    crs_bbox = {"xmin": crs_xmin, "ymin": crs_ymin, "xmax": crs_xmax, "ymax": crs_ymax}
    
    # RELATIVE BBOX NEEDS NEEDS TO REFERENCE FROM THE BOTTOM-LEFT INSTEAD OF TOP-LEFT. SWAP Y COODINATES
    relative_bbox = {"xmin":ulx, "ymin":uly, "xmax":lrx, "ymax":lry}
    #print("RELATIVE BBOX: ", relative_bbox)
    
    return relative_bbox, crs_bbox
  


# In[7]:


def findRaster(raster_files, polybbox):
    for raster_path in raster_files:
        rasterinfo = getRasterInfo(raster_path)
        rasterbbox = rasterinfo["extent"]

        if polybbox["xmin"] > rasterbbox["xmin"] and polybbox["xmax"] < rasterbbox["xmax"]:
            if polybbox["ymin"] > rasterbbox["ymin"] and polybbox["ymax"] < rasterbbox["ymax"]:
                return raster_path, rasterinfo
            
    #print("Unable to locate raster file for poly. Returning null values")
    return False


# In[8]:


def splitGeom(geom):
    x = geom.centroid.x
    y = geom.centroid.y
    return pd.Series([x, y], index=[["x", "y"]])


# In[23]:


def createShapefilePoints(polys, raster_list, out_dir, csv=False, shp=False, overwrite=False):
    start = datetime.now()
    
    rows_list = []
    
    with fiona.open(polys, "r") as poly_vector:
        print("NUM POLY FEATURES: ", len(poly_vector))
        
        # get crs
        out_crs = poly_vector.crs

        count = 0
        
        wavelengths = getRasterData(raster_list[0]) # just sample raster
        wavelength_names = []
        for name in wavelengths[0]:
            attrib_name = "b_" + str(name)
            wavelength_names.append(attrib_name)
        
        point_index = 0
        
        # BEGIN INTERATING POLYGONS
        for polygon in poly_vector:

            # MARTHA ONLY. Skip sub-site. E.g. don"t use "D1mesq2-2", only "D1mesq2"
            siteId = polygon["properties"]["siteID"]
            if "-" in siteId:
                continue

            print("Starting feature ", siteId)
            geometry = shape(polygon["geometry"])

            #(507897.81164529826, 3523674.995766839, 507904.882699674, 3523682.0668212175)
            geom_b = geometry.bounds

            poly_bbox = {}
            poly_bbox["xmin"], poly_bbox["ymin"], poly_bbox["xmax"], poly_bbox["ymax"] =  map(float, geom_b)

            ras_file, ras_info = findRaster(raster_list, poly_bbox)

            if not ras_info:
                continue

            raster_extent = ras_info["extent"]

            # geometry bounds lower-left and upper-right
            #ll = geom_b[0:2] # lower-left
            #ur = geom_b[2:4] # upper-right
            snapped_bboxes = getSnappedBbox(poly_bbox, ras_info)
            relative_bbox = snapped_bboxes[0]
            crs_bbox = snapped_bboxes[1]

            # outshape of the polygon bbox in pixels
            outshape_x = abs(relative_bbox["xmax"] - relative_bbox["xmin"])
            outshape_y = abs(relative_bbox["ymax"] - relative_bbox["ymin"])
            
            #print("RELATIVE BBOX: ", relative_bbox)
            if outshape_x == 0 or outshape_y == 0:
                raise ValueError("Problem with snapped bounding box ",
                                 outshape_x, outshape_y)

            raster_pixel_x_size = abs(ras_info["resx"])
            raster_pixel_y_size = abs(ras_info["resy"])

            half_x_size = raster_pixel_x_size / 2
            half_y_size = raster_pixel_y_size / 2

            polygon_internal_points = []
            
            raster_info = getRasterData(ras_file)
            scale_fact = ras_info["reflectance scale factor"]
            reflectance_data_ds = raster_info[1]

            #print("REFLECTANCE_DATA_DS: ", reflectance_data_ds)
            subset_reflectance = reflectance_data_ds[relative_bbox["ymin"]:relative_bbox["ymax"],relative_bbox["xmin"]:relative_bbox["xmax"],:] 

            #print("REFELCTANCE_DATA_DS: ", reflectance_data_ds)
            #print("\n SUBSET_REFLECTANCE: ", subset_reflectance)
            num_y, num_x, num_z = subset_reflectance.shape                                  # obtain the dimensions of the subsetted data
            num_rows = num_x*num_y                                                          # and the total number of z rows in the subset
         
            #print(ras_file)
            # iterate through row of pixels
            for x in range(outshape_x):
                pointx = (crs_bbox["xmin"] + half_x_size) + (x * raster_pixel_x_size)
                # iterate through column of row
                for y in range(outshape_y):
                    pointy = (crs_bbox["ymin"] + half_y_size) + (y * raster_pixel_y_size)
                    point = Point(pointx, pointy)
                    # check for intersection
                    if point.within(geometry):
                        count += 1
                        
                        point_dict = OrderedDict()
                        point_dict["geometry"] = point
                        
                        # SET POINT ATTRIBUTES BASE ON POLYGON ATTRIBUTES
                        for k,v in polygon["properties"].items():
                            point_dict[k] = v
                        
                        point_dict["raster"] = os.path.basename(ras_file)
                        
                        wavelength_data = subset_reflectance[y][x]
                        
                        for i in range(len(wavelength_data)):         # process each wavelength/z bucket
                            if wavelength_data[i] != ras_info["data ignore value"]:
                                value = wavelength_data[i]/scale_fact
                            else:
                                value = -9999
                            #print(wavelength_names[i], "---", value)
                            point_dict[wavelength_names[i]] = float("{0:.4f}".format(value))
                        
                        rows_list.append(point_dict)
                        
                        point_index += 1
                        
                        if point_index % 500 == 0:
                            print("\t- {} points created".format(point_index))
                
            if count == 0:
                raise ValueError("PROBLEM. NO POINTS CREATED FOR FEATURE - ", feature)

    # Create list of columns in the correct order. Passing list of dictionaries doesn"t preserve column order otherwise
    cols = []
    for k,v, in rows_list[0].items():
        cols.append(k)
        
    points_gdf = gpd.GeoDataFrame(rows_list, columns=cols)
    points_gdf.crs = out_crs
    
    #for i in wavelength_names:
    #    points_gdf[i] = pd.to_numeric(points_gdf[i])
    
    print("\nFinished creating points for %s. %d points created." % (polys, count))
    
    ofile_name = os.path.splitext(os.path.basename(polys))[0] + "_points"
    
    if csv:
        ocsv = ofile_name + ".csv"
        ocsv_path = os.path.join(out_dir, ocsv)
        print("\nWriting points to csv file: %s ..." % ocsv_path)
        # coordinates exist in column "geometry". break into x, y
        points_csv_gdf = points_gdf.copy()
        # initialize x and y columns
        points_csv_gdf["x"] = None
        points_csv_gdf["y"] = None
        points_csv_gdf[["x","y"]] = points_csv_gdf["geometry"].apply(splitGeom)
        del points_csv_gdf["geometry"]
        points_csv_gdf.to_csv(ocsv_path)
        
    if shp:
        oshp = ofile_name + ".shp"
        oshp_path = os.path.join(out_dir, oshp)
        print("\nWriting points to shapefile: %s ..." % oshp_path)
        points_shp_gdf = points_gdf.copy()
        points_shp_gdf["FID"] = points_shp_gdf.index
        points_shp_gdf.to_file(oshp_path, driver="ESRI Shapefile")
        
    return points_gdf


# In[24]:


df = createShapefilePoints(in_polygon, h5_files, out_dir, csv=True, shp=True, overwrite=True)
        
print("\n\nFINISHED")

