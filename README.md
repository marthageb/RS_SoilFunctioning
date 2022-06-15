# RS_SoilFunctioning
code/data to run ML on RS data (hyperspectral and other data types) to predict collective soil functioning across space

Corresponding Author for Code:   
  Martha Farella, farellam@iu.edu

License:
    MIT
---------------------------------------------
# Files and Descriptions
## hyperspec_foliarN - use hyperspectral data to predict foliar N across SRER
  brdf_topoR2py.py:
    apply NDVI mask, BRDF, and topographic corrections to individual flightlines. NEON flightlines accessed from the NEON data portal (https://data.neonscience.org/)
	hyperspec_info.py:
    add metadata fields to the .h5 topo/BRDF corrected reflectance data (output(s) from brdf_topoR2py.py)
	RGBhyperspec.R:
    convert .h5 file to RGB raster stack for original and BRDF corrected imagery mosaic all of the RGB flightlines together
	extract_pixel_hyperspec.R:
    extract corrected reflectance values for ENVI x/y cords (plant sample/plot locations)
	refagg.R:
    aggregate all of the extracted reflectance values together by taking the mean for each plot; combine reflectance values with plot averaged foliar chemistry values
	PCA_plot_hyperspec.R:
    select 10 pixels/plot, run PCA and calculate 95% confidence ellipse to ID outliers; reselect pixels until all 10 fall within the confidence ellipse. (iterate between ‘extract_pixel_hyperspec.R’ and ‘PCA_plot_hyperspec.R’)
	PLSRloop.R:
    perform brightness normalization on the corrected reflectance values; run PLSR in a loop storing model coefficients and performance metrics with each run; evaluate performance from all loop iterations and retain robust models; test model residuals for spatial auto correlation
	analyze_reflect.py:
    apply PLSR coefficients to the landscape (.h5 flightlines) for foliar N prediction. *Make sure brdf_topoR2py.py and hyperspec_info.py has been run on these flightlines before applying coefficients*
	analyze_overlap.py:
    calculate average values for all overlapping flightline areas
	
## RF_samplelocs - run RF analysis on data from sample locations. Perform variable reduction for downstream analysis
	3mbuff.R:
    create 3m2 buffer around soil sample locations, export as polygon feature class .shp file
	extract_hyperspec.py:
    extract hyperspectral reflectance readings for plot area
	ref_prepro.R:
    take outputs from extract_hyperspec.py and format for downstream analysis
	raster_extents.R:
    clip/extend each raster layer that will be used in RF analysis to align with the Foliar N data product. 
	extract_raster_vals.R:
    extract values for all RF raster layers at soil sample locations
	RFdataprepro.R:
    get data compiled/structured for Random Forest Analysis
	RFloop.R:
    run Random Forest on the data in a loop and output loop results (R2 and mse); VSURF package for variable reduction/importance analysis

## RF_SRER - run RF on the reduced variables to predict soil functionality (and individual components therein) across SRER
	RFstack.R:
    1. take NEON mosaicked hyperspec reflectance data, extract a single reflectance band from each tile, mosaic all tiles together into a single .tif;
    2. get extracted reflectance bands to “fit” the other raster layers (if needed); then, combine all predictor layers for each soil variable into a single raster stack and write results
	soil_pts_to_rast.R:
    convert soil lab measurements into raster with points/values corresponding to soil sample locations across SRER
	RF1hot.py:
    runs RF algorithms on rasterstack (output from RFstack.R) to predict desired soil variable across SRER. Before running RF, transforms categorical variables to OneHot variables
	RFoutput.R:
    take .tif output from RF1hot.py and transform raster so it has the same cs, extent, spatial coverage as all other raster layers. extract raster vals for soil sample locs and evalaute actual vs predicted vals
	
## data
	plotlvlN.csv:
    foliar N values for sample plots
	allsamples.shp:
    point locations for all samples
	soildata.csv:
    soil biogeochemical measurements
