import csv
import argparse
import sys
import os
import pdb
# Import GDAL, NumPy, and matplotlib
from osgeo import gdal, gdal_array, ogr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import mean_squared_error
import seaborn as sns
import logging
from tifffile import imwrite
import psutil
from datetime import datetime
import gc
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
t0 = datetime.now()

##### Parse arguments
parser = argparse.ArgumentParser(description='Reflectance Analyzer', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--profile',            type=str, default="om",                 help='run profile to be used')
parser.add_argument('--regress_small',      action="store_true",                    help='use small regression for testing')
parser.add_argument('--regress_verbose',    type=int, default=0,                    help='regression verbosity... see sklearn.ensemble.RandomForestRegressor')

arggroup = parser.add_mutually_exclusive_group()
arggroup.add_argument('--bands',        type=int, nargs='*',                        help="band numbers selected for x data")
arggroup.add_argument('--bandnames',    type=str, nargs='*',                        help="band names selected for x data")
args = parser.parse_args()
logging.info(f'***** Beginning Data Preparation Phase for profile "{args.profile}" *****')

profiles = {
        'cfi': {
        'xfile'         : "RFdataCFI.tif",
        'yfile'         : "CFI.tif",
        'output_tif'    : 'CFI.tif',
        'bandnames'     : [
            "soilage",
            "sand",
            "foliarN",
            "ref711",
            "slope",            
        ],
    },
    'om': {
        'xfile'         : "RFdataOM.tif",
        'yfile'         : "om.tif",
        'output_tif'    : 'om.tif',
        'bandnames'     : [
            "soilage",
            "sand",
            "landform",
            "X8moprecip",
            "slope",            
        ],
    },
    'doc': {
        'xfile'         : "RFdataDOC.tif",
        'yfile'         : "doc.tif",
        'output_tif'    : 'doc.tif',
        'bandnames'     : [
            "CHM",
            "AGbiomass",
            "foliarN",
            "ref706",
            "veg_class",
            "ref711",
            "silt",
            "ref716",
            "ref386",
            "ref2454",
        ],
    },
    'tdn': {
        'xfile'         : "RFdataTDN.tif",
        'yfile'         : "tdn.tif",
        'output_tif'    : 'tdn.tif',
        'bandnames'     : [
            "foliarN",
            "ref711",
        ],
    },
    'nh4': {
        'xfile'         : "RFdataNH4.tif",
        'yfile'         : "nh4.tif",
        'output_tif'    : 'nh4.tif',
        'bandnames'     : [
            "foliarN",
            "canopyht",
            "ref611",
            "ref541",
            "ref546",
            "ref551",
            "ref556",
            "ref531",
            "ref526",
            "ref1347",
            "slope",
        ],
    },
    'no3': {
        'xfile'         : "RFdataNO3.tif",
        'yfile'         : "no3.tif",
        'output_tif'    : 'no3.tif',
        'bandnames'     : [
            "ref711",
            "ref2474",
            "canopyht",
            "ref1813",
            "AGbiomass",
            "ref1357",
            "ref2509",
            "slope",
        ],
    },
    'ag': {
        'xfile'         : "RFdataAG.tif",
        'yfile'         : "ag.tif",
        'output_tif'    : 'ag.tif',
        'bandnames'     : [
            "slope",
            "silt",
            "sand",
            "landform",
        ],
    },
    'bg': {
        'xfile'         : "RFdataBG.tif",
        'yfile'         : "bg.tif",
        'output_tif'    : 'bg.tif',
        'bandnames'     : [
            "slope",
            "soilage",
            "PM",
            "ph",
            "landform",
        ],
    },
    'cb': {
        'xfile'         : "RFdataCB.tif",
        'yfile'         : "cb.tif",
        'output_tif'    : 'cb.tif',
        'bandnames'     : [
            "slope",
            "silt",
            "rock",
            "soil_age",
            "landform",
            "ref696",
            "ref2489",
            "ref1352",
            "ref1387",
            "ref691",
            "ref1132",
            "canopyht",
            "ref1127",
        ],
    },
    'xyl': {
        'xfile'         : "RFdataXYL.tif",
        'yfile'         : "xyl.tif",
        'output_tif'    : 'xyl.tif',
        'bandnames'     : [
            "slope",
            "ref1127",
            "silt",
            "ref1132",
            "ref1142",
        ],
    },
    'nag': {
        'xfile'         : "RFdataNAG.tif",
        'yfile'         : "nag.tif",
        'output_tif'    : 'nag.tif',
        'bandnames'     : [
            "silt",
            "ref2509",
            "ref2504",
            "precip8mo",
            "precip3mo",
            "ref1933",
            "elevation",
            "ref1893",
        ],
    },
    'lap': {
        'xfile'         : "RFdataLAP.tif",
        'yfile'         : "lap.tif",
        'output_tif'    : 'lap.tif',
        'bandnames'     : [
            "soil_age",
            "ref2489",
            "sand",
            "ref1863",
            "ref1843",
            "ref1913",
        ],
    },
    'phos': {
        'xfile'         : "RFdataPHOS.tif",
        'yfile'         : "phos.tif",
        'output_tif'    : 'phos.tif',
        'bandnames'     : [
            "precip3mo",
            "silt",
            "ref2509",
            "ref2504",
            "rock",
            "slope",
        ],
    },
    'mbc': {
        'xfile'         : "RFdataMBC.tif",
        'yfile'         : "mbc.tif",
        'output_tif'    : 'mbc.tif',
        'bandnames'     : [
            "soil_age",
            "AGbiomass",
            "silt",
            "slope",
            "ref691",
            "ref1833",
            "ref701",
            "ref386",
            "aspect",
            "ref1893",
            "ref571",
            "elevation",
        ],
    },
    'mbn': {
        'xfile'         : "RFdataMBN.tif",
        'yfile'         : "mbn.tif",
        'output_tif'    : 'mbn.tif',
        'bandnames'     : [
            "ref1377",
            "ref706",
            "precip8mo",
        ],
    },
   'dna': {
        'xfile'         : "RFdataDNA.tif",
        'yfile'         : "dna.tif",
        'output_tif'    : 'dna.tif',
        'bandnames'     : [
            "ref2509",
            "soil_age",
            "canopyht",
            "ref386",
            "silt",
            "ref391",
            "ref691",
        ],
    },
   'swc': {
        'xfile'         : "RFdataSWC.tif",
        'yfile'         : "swc.tif",
        'output_tif'    : 'swc.tif',
        'bandnames'     : [
            "ref2509",
            "soil_age",
            "landform",
            "elevation",
            "ref2494",
            "foliarN",
        ],
    },
   'ph': {
        'xfile'         : "RFdatapH.tif",
        'yfile'         : "ph.tif",
        'output_tif'    : 'ph.tif',
        'bandnames'     : [
            "BD",
            "X8moprecip",
            "X3moprecip",
            "clay",
            "elevation",
            "silt",
            "PM",
            "soilage",
            "X1yrprecip",
        ],
    },
    'type' : {
        'xfile'         : "RFdatatype.tif",
        'yfile'         : "type.tif",
        'output_tif'    : 'type.tif',
        'bandnames'     : [
            "CHM",
            "FoliarN",
            "aspect",
            "b344",
            "NEON_AGB",
        ],
    },
}

if args.profile not in profiles:
    logging.info(f'Defined profiles are: {profiles.keys()}')
    raise Exception(f'Specified profile "{args.profile}" not found in defined profiles')
profile = profiles[args.profile]

#define the bands in the x-dataframe
bandnames = profile['bandnames']

#define the factor variables
onehots = [
    "soilage",
    "landform",
    "PM",
    ]


sel_bands = []
if args.bands:
    sel_bands = args.bands
elif args.bandnames:
    for arg in args.bandnames:
        if arg in bandnames:
            sel_bands.append(bandnames.index(arg))
        else:
            raise ValueError(f'--bandnames name "{arg}" is not valid') 
 
# Tell GDAL to throw Python exceptions, allocate buffer space, and register all drivers
#gdal.UseExceptions()
bufmem = ((1000 ** 3) * 20) 
logging.info(f'Allocating {bufmem/(1000 ** 3):2.1f}G memory for gdal buffers')
gdal.SetCacheMax(bufmem)
gdal.AllRegister()

# Read in our image
logging.info(f"gdal.opening xfile {profile['xfile']}")
img_ds = gdal.Open(profile['xfile'], gdal.GA_ReadOnly)
ysize = img_ds.RasterYSize
xsize = img_ds.RasterXSize
bands_in_xfile = img_ds.RasterCount 
raster_dtype = gdal_array.GDALTypeCodeToNumericTypeCode(img_ds.GetRasterBand(1).DataType)

if bands_in_xfile != len(bandnames):
    raise ValueError(f'In "{profile["xfile"]}" there are {bands_in_xfile} bands which does not match defined bandnames {len(bandnames)}... cannot continue')

logging.info("Beginning band extraction from xfile...")
band_data_list = []
all_bandnames = []
for b in range(bands_in_xfile):
    if not sel_bands or b in sel_bands:
        bandname = bandnames[b]
        if bandname in onehots:
            logging.info(f'Processing {b+1}/{bands_in_xfile} band "{bandnames[b]}" as onehot...')  
            df = pd.DataFrame(data = img_ds.GetRasterBand(b + 1).ReadAsArray())
            logging.info("Calculating dummies and appending bands to band_data_list...")
            df[df < 0] = np.NaN
            nrow,ncol = df.shape
            dummies = pd.get_dummies(df.values.flatten())
            del df
            ndumrow,ndumcol = dummies.shape
            dummies = dummies.values.reshape(nrow, ncol, ndumcol)
            for i in range(ndumcol):
                band_data_list.append(dummies[:, :, i])
                all_bandnames.append(f'{bandname}-{i}')
            del dummies
            gc.collect()
        else:
            logging.info(f'Processing {b+1}/{bands_in_xfile} band "{bandnames[b]}" via GetRasterBand().ReadAsArray() and appending to band_data_list...')
            band_data_list.append(img_ds.GetRasterBand(b + 1).ReadAsArray())
            all_bandnames.append(bandname)
del img_ds
gc.collect()
logging.info("Completed band extraction.  Closed gdal xfile and creating img...")
img = np.zeros((ysize, xsize, len(band_data_list)), raster_dtype)
for band_data_idx, band_data in enumerate(band_data_list):
    img[:, :, band_data_idx] = band_data.copy('K')
del band_data_list
gc.collect()

yfile = f"./train_pts/{profile['yfile']}"
logging.info(f"img created.  gdal.opening yfile {yfile} and obtaining roi...")
roi_ds = gdal.Open(yfile, gdal.GA_ReadOnly)
roi = roi_ds.GetRasterBand(1).ReadAsArray().astype(np.float32)
del roi_ds
gc.collect()
logging.info("roi created.  Closed gdal yfile")
logging.info(f"***** Data Preparation Phase complete in {datetime.now()-t0}")

logging.info("***** Beginning Analysis Phase *****")
t1 = datetime.now()
# Find how many non-zero entries we have -- i.e. how many training data samples?
n_samples = (roi > 0).sum()
logging.info(f'We have {n_samples} samples')

# What are our classification labels?
labels = np.unique(roi[roi > 0])
logging.info('The training data include {n} classes: {classes}'.format(n=labels.size, classes=labels))
# We will need a "X" matrix containing our features, and a "y" array containing our labels
#     These will have n_samples rows
#     In other languages we would need to allocate these and them loop to fill them, but NumPy can be faster

X = img[roi > 0, :]  
y = roi[roi > 0]

logging.info(f'Our X matrix is sized: {X.shape}')
logging.info(f'Our y array is sized: {y.shape}')

# Initialize our model with "n_estimators" # of trees
if args.regress_small:
    max_depth = 2
    n_estimators = 20    
else:
    max_depth = None
    n_estimators = 1000 

rf = RandomForestRegressor(n_estimators=n_estimators, max_depth=max_depth, oob_score=True, criterion='mse', random_state=42, max_features='auto', n_jobs=-2, verbose=args.regress_verbose)
#if doing a classification
#rf = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, oob_score=True, criterion='mse', random_state=42, max_features='auto', n_jobs=-2, verbose=args.regress_verbose)

# Fit our model to training data
logging.info(f"Fitting model to training data via rf.fit with parameters {rf.get_params()}...")
rf = rf.fit(X, y)

logging.info(f'Our OOB prediction of accuracy is: {rf.oob_score_ * 100}')
logging.info(f'Our R2 value is: {rf.score(X,y) * 100}')
for bandname, imp in zip(all_bandnames, rf.feature_importances_):
    logging.info(f'Band {bandname} importance: {imp}')

# Take our full image and reshape into long 2d array (nrow * ncol, nband) for classification
new_shape = (img.shape[0] * img.shape[1], img.shape[2])
img_array = img[:, :, :].reshape(new_shape)
logging.info(f'Reshaped img from {img.shape} to {img_array.shape}')

logging.info("Predicting classes via rf.predict...")
# Now predict for each pixel
class_prediction = rf.predict(img_array)
#logging.info(f'Our mse is: {mean_squared_error(y, class_prediction)}')

logging.info("Reshaping classification map...")
# Reshape our classification map
class_prediction = class_prediction.reshape(img[:, :, 0].shape)

output_tif = f"./RFresults/{profile['output_tif']}"
logging.info(f"writing {output_tif}...")
#write the prediction across SRER as a .tif
imwrite(output_tif, class_prediction) 
logging.info(f"{output_tif} written")
logging.info(f"***** Analysis Phase complete in {datetime.now()-t1}")
logging.info(f"***** All Phases complete in {datetime.now()-t0}")

logging.info("Presenting plot...")
# Visualize
plt.subplot(111)
plt.imshow(class_prediction, interpolation='none')
plt.colorbar();
plt.show()