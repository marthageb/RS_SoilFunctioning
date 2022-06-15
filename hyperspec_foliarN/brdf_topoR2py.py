#take NEON .h5 reflectance data, perform topographic and BRDF correction based on .h5 reflectance and meta-data, write an output .h5 files with results

#modified code from python code written by Adam Chlus: https://github.com/EnSpec/HyTools-sandbox and R code written by Aaron Kamoske: https://github.com/akamoske/hypRspec
#TOPOGRAPHIC CORRECTION BASED ON: 
#   Soenen, S.A., Peddle, D.R., and Coburn, C.A., 2005. SCS+C: A Modified Sun-Canopy-Sensor Topographic Correction 
#     in Forested Terrain. IEEE Transactions on Geoscience and Remote Sensing, 43(9): 2148-2159.

# BRDF correction is based on the following papers:
#    Colgan, M.S., Baldeck, C.A., Feret, J.B., and Asner, G.P., 2012. Mapping savanna tree species at ecosystem scales 
#       using support vector machine classification and BRDF correction on airborne hyperspectral and LiDAR data.
#       Remote Sensing, 4(11): 3462-3480.  
#
#   Collings, S., Caccetta, P., Campbell, N., and Wu, X., 2010. Techniques for BRDF correction of hyperspectral mosaics. 
#       IEEE Transactions on Geoscience and Remote Sensing, 48(10): 3733-3746.
#
#   Schlapfer, D., Richter, R., and Feingersh, T., 2015. Operational BRDF effects correction for wide-field-of-view 
#       optical scanners (BREFCOR). IEEE Transactions on Geoscience and Remote Sensing, 53(4): 1855-1864.
#
#   Wanner, W., Li, X., and Strahler, A.H., 1995. On the derivation of kernels for kernel-driven models of 
#       bidirectional reflectance. Journal of Geophysical Research: Atmospheres, 100(D10): 21077-21089.
# 
#   Weyermann, J., Kneubuhler, M., Schlapfer, D., and Schaepman, M.E., 2015. Minimizing Reflectance Anisotropy 
#       Effects in Airborne Spectroscopy Data Using Ross-Li Model Inversion With Continuous Field Land Cover Stratification. 
#       IEEE Transactions on Geoscience and Remote Sensing, 53(11): 5814-5823.


#python brdf_topoR2py.py --ross thick --li sparse --hy_file NEON_D14_SRER_DP1_20170824_162407_reflectance.h5 NEON_D14_SRER_DP1_20170824_163106_reflectance.h5 NEON_D14_SRER_DP1_20170825_165121_reflectance.h5 NEON_D14_SRER_DP1_20170825_170642_reflectance.h5 NEON_D14_SRER_DP1_20170829_181113_reflectance.h5
#Center hdf files: NEON_D14_SRER_DP1_20170824_171440_reflectance.h5 NEON_D14_SRER_DP1_20170824_172149_reflectance.h5 NEON_D14_SRER_DP1_20170824_172924_reflectance.h5   NEON_D14_SRER_DP1_20170829_174731_reflectance.h5 NEON_D14_SRER_DP1_20170824_173621_reflectance.h5 NEON_D14_SRER_DP1_20170824_175127_reflectance.h5 NEON_D14_SRER_DP1_20170824_175910_reflectance.h5 NEON_D14_SRER_DP1_20170824_180609_reflectance.h5 NEON_D14_SRER_DP1_20170824_181344_reflectance.h5 NEON_D14_SRER_DP1_20170824_182059_reflectance.h5 NEON_D14_SRER_DP1_20170824_182820_reflectance.h5 NEON_D14_SRER_DP1_20170824_183531_reflectance.h5 NEON_D14_SRER_DP1_20170824_200520_reflectance.h5
#Northwest hdf files: NEON_D14_SRER_DP1_20170824_201237_reflectance.h5 NEON_D14_SRER_DP1_20170824_203449_reflectance.h5 NEON_D14_SRER_DP1_20170824_204230_reflectance.h5 NEON_D14_SRER_DP1_20170824_204947_reflectance.h5 NEON_D14_SRER_DP1_20170824_205924_reflectance.h5 NEON_D14_SRER_DP1_20170824_210427_reflectance.h5 NEON_D14_SRER_DP1_20170824_210936_reflectance.h5 NEON_D14_SRER_DP1_20170824_211426_reflectance.h5 NEON_D14_SRER_DP1_20170824_211937_reflectance.h5 NEON_D14_SRER_DP1_20170824_212429_reflectance.h5 NEON_D14_SRER_DP1_20170824_213002_reflectance.h5 NEON_D14_SRER_DP1_20170824_213511_reflectance.h5 NEON_D14_SRER_DP1_20170825_162024_reflectance.h5 NEON_D14_SRER_DP1_20170825_162519_reflectance.h5 NEON_D14_SRER_DP1_20170825_163020_reflectance.h5
#south hdf files: NEON_D14_SRER_DP1_20170829_184602_reflectance.h5 NEON_D14_SRER_DP1_20170829_190359_reflectance.h5 NEON_D14_SRER_DP1_20170829_182805_reflectance.h5 NEON_D14_SRER_DP1_20170829_181113_reflectance.h5 NEON_D14_SRER_DP1_20170829_175512_reflectance.h5 NEON_D14_SRER_DP1_20170829_174013_reflectance.h5 NEON_D14_SRER_DP1_20170829_172625_reflectance.h5
#MAIN hdf files that cover sample locs: NEON_D14_SRER_DP1_20170824_161655_reflectance.h5 NEON_D14_SRER_DP1_20170824_162407_reflectance.h5 NEON_D14_SRER_DP1_20170824_163106_reflectance.h5 NEON_D14_SRER_DP1_20170824_170703_reflectance.h5 NEON_D14_SRER_DP1_20170824_174415_reflectance.h5 NEON_D14_SRER_DP1_20170824_195807_reflectance.h5 NEON_D14_SRER_DP1_20170824_202008_reflectance.h5 NEON_D14_SRER_DP1_20170824_202731_reflectance.h5 NEON_D14_SRER_DP1_20170825_165121_reflectance.h5 NEON_D14_SRER_DP1_20170825_165904_reflectance.h5 NEON_D14_SRER_DP1_20170829_161226_reflectance.h5 NEON_D14_SRER_DP1_20170829_180310_reflectance.h5 NEON_D14_SRER_DP1_20170829_181113_reflectance.h5
#northeast hdf files: NEON_D14_SRER_DP1_20170824_164049_reflectance.h5 NEON_D14_SRER_DP1_20170824_164453_reflectance.h5 NEON_D14_SRER_DP1_20170824_165157_reflectance.h5 NEON_D14_SRER_DP1_20170824_165936_reflectance.h5 NEON_D14_SRER_DP1_20170824_171440_reflectance.h5 NEON_D14_SRER_DP1_20170824_172149_reflectance.h5 NEON_D14_SRER_DP1_20170824_172924_reflectance.h5 NEON_D14_SRER_DP1_20170825_170642_reflectance.h5 NEON_D14_SRER_DP1_20170825_171420_reflectance.h5 NEON_D14_SRER_DP1_20170825_172210_reflectance.h5 NEON_D14_SRER_DP1_20170825_172943_reflectance.h5 NEON_D14_SRER_DP1_20170825_173703_reflectance.h5 NEON_D14_SRER_DP1_20170825_174439_reflectance.h5 NEON_D14_SRER_DP1_20170825_175314_reflectance.h5 NEON_D14_SRER_DP1_20170829_174731_reflectance.h5
#OTHER hdf files that cover sample locs: NEON_D14_SRER_DP1_20170824_201237_reflectance.h5 NEON_D14_SRER_DP1_20170824_214058_reflectance.h5 NEON_D14_SRER_DP1_20170829_181113_reflectance.h5


import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from sklearn import linear_model
import gc
import pdb
import psutil
from datetime import datetime

parser = argparse.ArgumentParser()

#will need to change the 'site name' for various sites
parser.add_argument('--hy_file', type=str, nargs = '*', help='hdf5 file containing hyperspectral imagery and associated metadata')
parser.add_argument("--ndvi_mask", help="NDVI mask created with the .ndvi_mask function", action='store_true')
parser.add_argument("--brightness_mask", help=" brightness mask created with the .brightness_mask function", action='store_true')
parser.add_argument('--reflectance_path', type=str, default="SRER/Reflectance/Reflectance_Data", help="hdf5 path to reflectance data")
parser.add_argument('--wavelength_path', type=str, default="SRER/Reflectance/Metadata/Spectral_Data/Wavelength", help="hdf5 path to wavelength metadata")
parser.add_argument('--solar_az_path', type=str, default="SRER/Reflectance/Metadata/Logs/Solar_Azimuth_Angle", help="hdf5 path to solar azimuth data")
parser.add_argument('--solar_zn_path', type=str, default="SRER/Reflectance/Metadata/Logs/Solar_Zenith_Angle", help="hdf5 path to solar zenith data")
parser.add_argument('--slope_path', type=str, default="SRER/Reflectance/Metadata/Ancillary_Imagery/Slope", help="hdf5 path to slope data")
parser.add_argument('--aspect_path', type=str, default="SRER/Reflectance/Metadata/Ancillary_Imagery/Aspect", help="hdf5 path to aspect data")
parser.add_argument('--sensor_az_path', type=str, default="SRER/Reflectance/Metadata/to-sensor_Azimuth_Angle", help="hdf5 path to sensor azimuth data")
parser.add_argument('--sensor_zn_path', type=str, default="SRER/Reflectance/Metadata/to-sensor_Zenith_Angle", help="hdf5 path to sensor zenith data")
parser.add_argument('--coordinate_path', type=str, default="SRER/Reflectance/Metadata/Coordinate_System/Map_Info", help="hdf5 path to hdf5 path to coordinate data")
parser.add_argument("--ross", help=".ross kernel types-- thick or thin", type =str)
parser.add_argument("--li", help="Li kernel types-- dense or sparse", type =str)

args = parser.parse_args()
if len(args.hy_file) == 1:
    image = args.hy_file[0]
for image in args.hy_file:
    outfile = f'output/{image}'

    #load .h5 file
    f = h5py.File(image, 'r')

    # lets load the reflectance data
    refl_info = f[args.reflectance_path]

    # lets read in the wavelength info
    wavelengths = f[args.wavelength_path] 

    # lets save the dimensions of the dataset for future use
    n_rows = refl_info.shape[0]
    n_cols = refl_info.shape[1]
    n_bands = refl_info.shape[2]
    print(n_bands)

    # save the scale factor and the data ignore value
    scaleFactor = refl_info.attrs['Scale_Factor']
    noDataVal = refl_info.attrs['Data_Ignore_Value']
    #scaleFactor = 10000.0
    #noDataVal = -9999.0
    
    #---------------------------------------------------------------------------------------------------
    # lets save all the information that will not need to be repeated for each band
    #---------------------------------------------------------------------------------------------------

    # save several pieces of data as numpy arrays
    sensor_az = f.get(args.sensor_az_path).value
    sensor_zn = f.get(args.sensor_zn_path).value
    solar_az = f.get(args.solar_az_path).value
    solar_zn = f.get(args.solar_zn_path).value
    slope = f.get(args.slope_path).value
    aspect = f.get(args.aspect_path).value


    #  need to remove the no data values from all of these
    slope[slope == int(noDataVal)] = np.nan
    aspect[aspect == int(noDataVal)] = np.nan
    #solar_az[solar_az == int(noDataVal)] = np.nan #Are these ever going to = noDataVal?
    #solar_zn[solar_zn == int(noDataVal)] = np.nan #Are these ever going to = noDataVal?
    sensor_az[sensor_az == int(noDataVal)] = np.nan
    sensor_zn[sensor_zn == int(noDataVal)] = np.nan

    # now we need to convert these to radians
    slope = np.deg2rad(slope)
    aspect = np.deg2rad(aspect)
    solar_az = np.deg2rad(solar_az)
    solar_zn = np.deg2rad(solar_zn)
    sensor_az = np.deg2rad(sensor_az)
    sensor_zn = np.deg2rad(sensor_zn)
    
    
    #---------------------------------------------------------------------------------------------------
    # create the NDVI mask
    #---------------------------------------------------------------------------------------------------
    #define the red and nir band
    red_nm = 668
    nir_nm = 828
    #define the ndvi_threshold
    ndvi_threshold = -1
    
    #find the index for the closest band to red and nir 
    wl = f.get(args.wavelength_path).value
    nir_index = np.argmin(np.abs(wl-nir_nm))
    red_index = np.argmin(np.abs(wl-red_nm))
    
    #load the bands
    red = refl_info[:,:,red_index]
    nir  = refl_info[:,:,nir_index]
    
    #calculate ndvi
    ndvi = (nir - red) / (nir + red)
    
    #create the mask
    ndvi_mask = (ndvi > ndvi_threshold) & (nir != noDataVal)

    #---------------------------------------------------------------------------------------------------
    # calculate the topographic correction coefficients
    #---------------------------------------------------------------------------------------------------

    print(f"{datetime.now()} calculating topographic correction variables.")

    # Generate the cosine i
    # the cosine of the incidence angle (i ), defined as the angle between the normal to the pixel 
    # surface and the solar zenith direction
    rel_az = aspect - solar_az
    cosine_i = np.cos(solar_zn) * np.cos(slope) + np.sin(solar_zn) * np.sin(slope) * np.cos(rel_az)

    # apply the NDVI mask to the cosine_i
    cosine_i[~ndvi_mask] = np.nan 


    # Eq 11. Soenen et al., IEEE TGARS 2005
    # cos(alpha)* cos(theta)
    # alpha -- slope (slope), theta -- solar zenith angle (solar_zn)
    c1 = np.cos(solar_zn) * np.cos(slope)

    #---------------------------------------------------------------------------------------------------
    # calculate the ross Volumetric Scattering Kernel 
    #---------------------------------------------------------------------------------------------------

    print(f"{datetime.now()} calculating brdf correction variables.")

    # calculate the relative azimuth 
    relative_az = sensor_az - solar_az

    # first we need to use Eq 2. from Schlapfer et al. IEEE-TGARS 2015 which uses the inverse cosine
    phase = np.arccos(np.cos(solar_zn) * np.cos(sensor_zn) + np.sin(solar_zn) * np.sin(sensor_zn) * np.cos(relative_az))

    if args.ross == 'thick':
        # for the Thick Ross Kernel - Eq 7. Wanner et al. JGRA 1995
        kvol = ((np.pi/2 - phase) * np.cos(phase) + np.sin(phase)) / (np.cos(sensor_zn) + np.cos(solar_zn)) - np.pi/4
    elif args.ross == 'thin':
        # for the Thin ross Kernal - Eq 13. Wanner et al. JGRA 1995
        kvol = ((np.pi/2 - phase) * np.cos(phase) + np.sin(phase)) / (np.cos(sensor_zn) * np.cos(solar_zn)) - np.pi/2

    #apply the NDVI mask to kvol
    kvol[~ndvi_mask] = np.nan
    
    #---------------------------------------------------------------------------------------------------
    # calculate the ross Volumetric Scattering Kernel at nadir (sensor_zn = 0)
    #---------------------------------------------------------------------------------------------------

    # first we need to use Eq 2. from Schlapfer et al. IEEE-TGARS 2015 which uses the inverse cosine
    phase_nad = np.arccos(np.cos(solar_zn) * np.cos(0) + np.sin(solar_zn) * np.sin(0) * np.cos(relative_az))

    if args.ross == 'thick':
        # for the Thick args.ross Kernel - Eq 13. Wanner et al. JGRA 1995
        kvol_nad = ((np.pi/2 - phase_nad) * np.cos(phase_nad) + np.sin(phase_nad)) / (np.cos(0) + np.cos(solar_zn)) - np.pi/4
    elif args.ross == 'thin':
        # for the Thin args.ross Kernal - Eq 13. Wanner et al. JGRA 1995
        kvol_nad = ((np.pi/2 - phase_nad) * np.cos(phase_nad) + np.sin(phase_nad)) / (np.cos(0) * np.cos(solar_zn)) - np.pi/2
    
    #apply the NDVI mask to kvol_nad
    kvol_nad[~ndvi_mask] = np.nan
    #---------------------------------------------------------------------------------------------------
    # calculate the Li Geometric Scattering Kernel

    # we will use the constants from Colgan et al. 2012 in Remote Sensing
    # h/b = 2 ; b/r = 10
    #---------------------------------------------------------------------------------------------------

    # first we need to implement Eq. 37,52. Wanner et al. JGRA 1995
    solar_zn_at = np.arctan(10 * np.tan(solar_zn))
    sensor_zn_at = np.arctan(10 * np.tan(sensor_zn))

    # next we need to use Eq 50. Wanner et al. JGRA 1995
    d = np.sqrt((np.tan(solar_zn_at) ** 2) + (np.tan(sensor_zn_at) ** 2) - (2 * np.tan(solar_zn_at) * np.tan(sensor_zn_at) * np.cos(relative_az)))

    # next we need to use Eq 49. Wanner et al. JGRA 1995
    # we will need to restraint these values between -1 and 1 to not return NaN and Pi values
    # for more info see this stack exchange thread 
    # https://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
    t_num = 2 * np.sqrt(d**2 + (np.tan(solar_zn_at) * np.tan(sensor_zn_at) * np.sin(relative_az)) ** 2)
    t_denom = (1 / np.cos(solar_zn_at)) + (1 / np.cos(sensor_zn_at))
    t = np.arccos(np.clip(t_num/t_denom, -1, 1))

    # next we need to use Eq 33,48. Wanner et al. JGRA 1995
    o = (1 / np.pi) * (t - np.sin(t) * np.cos(t)) * t_denom

    # next we need to use Eq 51. Wanner et al. JGRA 1995
    cos_phase = np.cos(solar_zn_at) * np.cos(sensor_zn_at) + np.sin(solar_zn_at) * np.sin(sensor_zn_at) * np.cos(relative_az)

    if args.li == 'sparse':
        # for the Sparse Li Kernel - Eq 32. Wanner et al. JGRA 1995
        kgeom = o - (1 / np.cos(solar_zn_at)) - (1 / np.cos(sensor_zn_at)) + 0.5 * (1 + cos_phase) * (1 / np.cos(sensor_zn_at))
    elif args.li == 'dense':
        # for the Dense Li Kernel - Eq 47. Wanner et al. JGRA 1995
        kgeom = (((1 + cos_phase) * (1 / np.cos(sensor_zn_at))) / (t_denom - o)) - 2

    #apply the NDVI mask to kgeom
    kgeom[~ndvi_mask] = np.nan
    
    #---------------------------------------------------------------------------------------------------
    # calculate the Li Geometric Scattering Kernel at nadir (sensor.zn = 0)

    # we will use the constants from Colgan et al. 2012 in Remote Sensing
    # h/b = 2 ; b/r = 10
    #---------------------------------------------------------------------------------------------------

    # first we need to implement Eq. 37,52. Wanner et al. JGRA 1995
    solar_zn_at = np.arctan(10 * np.tan(solar_zn))
    sensor_zn_at = np.arctan(10 * np.tan(0))

    # use Eq 50. Wanner et al. JGRA 1995
    d = np.sqrt((np.tan(solar_zn_at) ** 2) + (np.tan(sensor_zn_at) ** 2) - (2 * np.tan(solar_zn_at) * np.tan(sensor_zn_at) * np.cos(relative_az)))
    
    # next we need to use Eq 49. Wanner et al. JGRA 1995
    # we will need to restraint these values between -1 and 1 to not return NaN and Pi values
    # for more info see this stack exchange thread 
    # https://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
    t_num = np.sqrt(d**2 + (np.tan(solar_zn_at) * np.tan(sensor_zn_at) * np.sin(relative_az)) ** 2)
    t_denom = (1 / np.cos(solar_zn_at)) + (1 / np.cos(sensor_zn_at))
    t = np.arccos(np.clip(2*(t_num/t_denom), -1, 1))

    # next we need to use Eq 33,48. Wanner et al. JGRA 1995
    o = (1 / np.pi) * (t - np.sin(t) * np.cos(t)) * t_denom
   
   # next we need to use Eq 51. Wanner et al. JGRA 1995
    cos_phase = np.cos(solar_zn_at) * np.cos(sensor_zn_at) + np.sin(solar_zn_at) * np.sin(sensor_zn_at) * np.cos(relative_az)
    if args.li == 'sparse':
        # for the Sparse Li Kernel - Eq 32. Wanner et al. JGRA 1995
        kgeom_nad = o - (1 / np.cos(solar_zn_at)) - (1 / np.cos(sensor_zn_at)) + 0.5 * (1 + cos_phase) * (1 / np.cos(sensor_zn_at))
    elif args.li == 'dense':
        # for the Dense Li Kernel - Eq 47. Wanner et al. JGRA 1995
        kgeom_nad = (((1 + cos_phase) * (1 / np.cos(sensor_zn_at))) / (t_denom - o)) - 2
    
    #apply the NDVI mask to kgeom_nad
    kgeom_nad[~ndvi_mask] = np.nan

    #---------------------------------------------------------------------------------------------------
    # clear up some memory since the only things we need from above are the kernels
    #---------------------------------------------------------------------------------------------------
    del red_nm
    del nir_nm
    del ndvi_threshold
    del wl
    del nir_index
    del red_index
    del red
    del nir
    del ndvi   
    del cos_phase
    del d
    del o
    del phase_nad
    del relative_az
    del sensor_az
    del sensor_zn
    del t
    del t_num
    del t_denom
    del sensor_zn_at
    del solar_zn_at
    del phase
    del rel_az
    del slope
    del aspect
    del solar_az
    del solar_zn
    gc.collect()

    #---------------------------------------------------------------------------------------------------
    # create vegetation masks for the kernels we will use and then rearrange the data
    # for regression
    #---------------------------------------------------------------------------------------------------
    #apply the NDVI mask to the kernels 
    #ross_mask = kvol[ndvi_mask]
    #li_mask = kgeom[ndvi_mask]
    #ross_mask_nad = kvol_nad[ndvi_mask]
    #li_mask_nad = kgeom_nad[ndvi_mask]
    ross_mask = kvol
    li_mask = kgeom
    # transform the data into the appropriate shape for regression
    #ross = np.expand_dims(ross_mask, axis = 1)
    #li = np.expand_dims(li_mask, axis = 1)
    #x_brdf = np.concatenate([ross, li, np.ones(li.shape)], axis = 1)

    #---------------------------------------------------------------------------------------------------
    # load in and apply corrections the reflectance imagery
    #---------------------------------------------------------------------------------------------------
    print(f'{datetime.now()} reading in reflectance data')
    #read the reflectance data
    reflClean = f.get(args.reflectance_path).value.astype(float)

    #memory clean up
    f.close()
    del f
    gc.collect()

    #Apply no data value
    arr_size = reflClean.shape
    print(f'{datetime.now()} applying nan values to reflectance bands')
    #print('% No Data: ',np.round(np.count_nonzero(reflClean==noDataVal)*100/(arr_size[0]*arr_size[1]*arr_size[2]),1))
    #nodata_ind = np.where(reflClean==noDataVal)
    #reflClean[nodata_ind]=np.nan 
    reflClean[reflClean == noDataVal] = np.nan
    
    
    #memory clean up
    #del nodata_ind
    #gc.collect()

    #Apply scale factor
    print(f'{datetime.now()} applying the scale factor')
    reflArray = reflClean/scaleFactor

    #memory clean up
    del reflClean
    gc.collect()

    #also could apply the noDataVal and scaleFactor band x band in the loop rather than doing all bands at once. 

    #---------------------------------------------------------------------------------------------------
    # apply the topographic correction to each reflectance band
    #---------------------------------------------------------------------------------------------------

    #apply the topgraphic and brdf corrections band x band
    print(f'{datetime.now()} applying topo correction')
    #create an empty lists to store the topographic coefficients
    #topo_coeff = []
    #topo_applied = []
    #brdf_coeff = []
    #brdf_applied = []
    topo_coeff = np.ndarray(shape=(426,1), dtype = np.float16)
    topo_applied = np.ndarray(shape=(n_bands, n_rows, n_cols), dtype = np.float16)
    brdf_coeff = np.ndarray(shape = (426,3), dtype = np.float16)
    brdf_applied = np.ndarray(shape= (n_bands, n_rows, n_cols), dtype = np.float16)


    #apply the topgraphic and brdf corrections band x band
    for i in range(n_bands):
        #select one band form the corrected reletance Array
        reflMatrix = reflArray[:,:,i]
        #apply the NDVI mask
        reflMatrix[~ndvi_mask] = np.nan
        #ID cells with missing data in the reflectance band
        missing = np.isnan(reflMatrix)
        #remove the 'missing val cells' from the x (cosine_i) and y (reflMatrix) datasets. This also converts the data from 2D from 1D array needed for the linear regression
        y = reflMatrix[~missing]
        x = cosine_i[~missing]
        # run the regression -- Eq 7. Soenen et al., IEEE TGARS 2005
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        # Eq 8. Soenen et al., IEEE TGARS 2005
        c = intercept/slope
        # Set a large number if slope is zero
        if not np.isfinite(c):
            c = 100000.0
        # save the coefficients
        #topo_coeff.append(c)
        print(type(c))
        topo_coeff[i,:] = c
        
        # find correction factor - Eq 11. Soenen et al., IEEE TGARS 2005
        cor_fact = (c1 + c) / (cosine_i + c)
        # apply the correction factor
        topo_cor = reflMatrix * cor_fact
        #write the corrected topgraphic correction values to an output
        #topo_applied.append(topo_cor)
        topo_applied[i,:,:]=(topo_cor)
        
        
        #---------------------------------------------------------------------------------------------------
        # apply the brdf correction to the topo corrected band
        #---------------------------------------------------------------------------------------------------

        print(f'{datetime.now()} applying brdf correction to band {i}')

        # apply the brightness mask to this topo corrected band
        #topo_matrix <- ifelse(args.brightness_mask, topo.cor, NA)

        #ID cells with missing data in the topographically corrected band
        missing = np.isnan(topo_cor)
        #remove the 'missing val cells' from the x (topo corrected band) and y (reflMatrix) datasets. This also converts the data from 2D from 1D array needed for the linear regression
        y = reflMatrix[~missing]
        x_kvol = kvol[~missing]
        x_kgeom = kgeom[~missing]
        #combine both IVs into a single dataframe
        x_brdf = x_kvol, x_kgeom
        #transpose x variables to match y variables
        x = np.transpose(x_brdf)
        
        # run the regression now
        reg = linear_model.LinearRegression()
        reg.fit(x,y)
        #get and save the regression coefficients and y intercept value
        coef = reg.coef_
        yint = reg.intercept_
        
        #brdf_coeff.append([yint, coef[0], coef[1]])
        brdf_coeff[i,:] = ([yint, coef[0], coef[1]])
        # memory management
        del y
        gc.collect()

        # apply the coefficients to the band - eq 5. Weyermann et al. IEEE-TGARS 2015
        brdf =  yint + (kgeom * coef[1]) + (kvol * coef[0])
        brdf_nad = yint + (kgeom_nad * coef[1]) + (kvol_nad * coef[0])

        # find the correction factor: eq 4. Weyermann et al. IEEE-TGARS 2015
        brdf_cor = brdf_nad / brdf

        # apply the correction factor to the band
        band_brdf = topo_cor * brdf_cor
        
        # first we need to rescale the data so we don't have floating points
        band_brdf_scale = np.round(scaleFactor * band_brdf)

        # next we can replace the NA values with the data ignore value
        band_brdf_scale[np.isnan(band_brdf_scale)] = noDataVal

        #write the corrected topgraphic correction values to an output
        #brdf_applied.append(band_brdf_scale)
        brdf_applied[i,:,:] =(band_brdf_scale)

        print(f'band {i} is finished')

    #---------------------------------------------------------------------------------------------------
    # clear up some memory
    #---------------------------------------------------------------------------------------------------
    del ndvi_mask
    del refl_info
    del kvol_nad
    del kgeom_nad   
    del ross_mask
    del li_mask
    del arr_size
    del reflArray
    del reflMatrix
    del missing
    del x
    del c
    del cor_fact
    del topo_cor
    del x_kvol
    del x_kgeom
    del x_brdf
    del reg
    del coef
    del yint
    del brdf
    del brdf_nad
    del brdf_cor
    del band_brdf
    del band_brdf_scale
    gc.collect()

    #print(json.dumps(dict(psutil.virtual_memory()._asdict()),indent=2)

    #create a new hdf output file  
    output = h5py.File(outfile,'w')

    #write output brdf datasets and move the "bands" axis
    print(f'{datetime.now()} writing brdf  coeff results')
    output.create_dataset("BRDF/Coeffs", data = brdf_coeff)
    print(f'{datetime.now()} flipping the brdf dataset axis')
    print(brdf_applied.shape)
    brdf_applied = np.moveaxis(brdf_applied, 0, 2)
    print(f'{datetime.now()} done flipping brdf axis. now writing the dataset')
    output.create_dataset("BRDF/Correction", data = brdf_applied, compression = "szip")
    #memory clean up
    del brdf_applied
    del brdf_coeff
    gc.collect()

    #write the topo datasets and move the axis
    #print(f'{datetime.now()} writing topo coeff results')
    #output.create_dataset("Topo/Coeffs", data = topo_coeff)
    #print(f'{datetime.now()} flipping the topo dataset axis')
    #topo_applied = np.moveaxis(topo_applied, 0, 2)
    #print(f'{datetime.now()} done flipping the topo axis. now writing the dataset')
    #output.create_dataset("Topo/Correction", data = topo_applied, compression = "szip")
    print(f'{datetime.now()} done')
    output.close()