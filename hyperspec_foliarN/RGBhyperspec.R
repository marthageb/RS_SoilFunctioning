###Code adopted from NEON tutorial: https://www.neonscience.org/intro-hsi-r-series###
###Take .h5 file, extract relevant information such as extent, no data vals, scale factor###
###Extract RGB bands from .h5 file, create raster stack and write raster .tiff file###

library(rhdf5)
library(raster)
library(rgdal)
library(gdalUtils)

setwd("Z:/SRER/Martha/hyperspec/flightlines")

f <- "./output/NDVI_thresh_0.1/NEON_D14_SRER_DP1_20170829_175512_reflectance.h5"

# define coordinate reference system from the EPSG code provided in the HDF5 file
myEPSG <- h5read(f,"/SRER/Reflectance/Metadata/Coordinate_System/EPSG Code" )
myCRS <- crs(paste0("+init=epsg:",myEPSG))

# get the Reflectance_Data attributes
#reflInfo <- h5readAttributes(f,"/SRER/Reflectance/Reflectance_Data" ) # for orig reflectance data
reflInfo <- h5readAttributes(f,"/BRDF/Correction" ) # for BRDF corrected reflectance data


# Grab the UTM coordinates of the spatial extent
xMin <- reflInfo$Spatial_Extent_meters[1]
xMax <- reflInfo$Spatial_Extent_meters[2]
yMin <- reflInfo$Spatial_Extent_meters[3]
yMax <- reflInfo$Spatial_Extent_meters[4]

# define the extent (left, right, top, bottom)
rasExt <- extent(xMin,xMax,yMin,yMax)

# view the extent to make sure that it looks right
rasExt

# Define the no data and scale factor 
myNoDataValue <- -9999
Scale_Factor <- 10000


#define the function to extract one band from the reflectance data
# file: the hdf file
# band: the band you want to process
# returns: a matrix containing the reflectance data for the specific band
band2Raster <- function(file, band, noDataValue, extent, CRS){
  # first, read in the raster
  #out <- h5read(file,"/SRER/Reflectance/Reflectance_Data",index=list(band,NULL,NULL)) # for the original reflectance data
  out <- h5read(file,"/Topo/Correction",index=list(band,NULL,NULL)) # for the BRDF corrected reflectance data
  
  # Convert from array to matrix
  out <- (out[1,,])
  # transpose data to fix flipped row and column order 
  # depending upon how your data are formatted you might not have to perform this
  # step.
  out <- t(out)
  # assign data ignore values to NA
  # note, you might chose to assign values of 15000 to NA
  out[out == myNoDataValue] <- NA
  
  # turn the out object into a raster
  outr <- raster(out,crs=CRS)
  
  # assign the extents to the raster
  extent(outr) <- extent
  
  # return the raster object
  return(outr)
}

# create a list of the bands we want in our stack
rgb <- list(58,34,19)

# lapply tells R to apply the function to each element in the list
rgb_rast <- lapply(rgb,FUN=band2Raster, file = f,
                   noDataValue=myNoDataValue, 
                   extent=rasExt,
                   CRS=myCRS)

# check out the properties or rgb_rast
# note that it displays properties of 3 rasters.
rgb_rast

# create a raster stack from our list of rasters
rgbStack <- stack(rgb_rast)

# Create a list of band names
bandNames <- paste("Band_",unlist(rgb),sep="")

# set the rasterStack's names equal to the list of bandNames created above
names(rgbStack) <- bandNames

# check properties of the raster list - note the band names
rgbStack

# scale the data as specified in the reflInfo$Scale Factor
#rgbStack <- rgbStack/as.integer(reflInfo$Scale_Factor)

# create a 3 band RGB image
#plotRGB(rgbStack, r=1,g=2,b=3, stretch = "lin")

# write out final raster    
writeRaster(rgbStack, file="./RGB/BRDF/SRER212429brdf.tif", format="GTiff")

####################################################################################
##############################MOSAIC RASTERS TOGETHER###############################
####################################################################################
r1 <- stack("./RGB/BRDF/SRER201237BRDF.tif")
r2 <- stack("./RGB/BRDF/SRER202008BRDF.tif")
r3 <- stack("./RGB/BRDF/SRER202731BRDF.tif")
r4 <- stack("./RGB/BRDF/SRER203449BRDF.tif")
r5 <- stack("./RGB/BRDF/SRER204230BRDF.tif")
r6 <- stack("./RGB/BRDF/SRER204947BRDF.tif")
r7 <- stack("./RGB/BRDF/SRER205924BRDF.tif")
r8 <- stack("./RGB/BRDF/SRER210427BRDF.tif")
r9 <- stack("./RGB/BRDF/SRER210936BRDF.tif")
r10 <- stack("./RGB/BRDF/SRER211426BRDF.tif")
r11 <- stack("./RGB/BRDF/SRER211937BRDF.tif")
r12 <- stack("./RGB/BRDF/SRER212429BRDF.tif")
r13 <- stack("./RGB/BRDF/SRER213002BRDF.tif")
r14 <- stack("./RGB/BRDF/SRER213511BRDF.tif")
r15 <- stack("./RGB/BRDF/SRER162024BRDF.tif")
r16 <- stack("./RGB/BRDF/SRER162519BRDF.tif")
r17 <- stack("./RGB/BRDF/SRER163020BRDF.tif")

#assign NA val
NAvalue(r1) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r2) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r3) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r4) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r5) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r6) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r7) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r8) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r9) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r10) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r11) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r12) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r13) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r14) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r15) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r16) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999
NAvalue(r17) <- -0.9999 #for BRDF correction: -0.9999; for orig image: -0.9999

#create a list of the rasters, define the function and NA.action
x <- list(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17)
x$fun <- mean
x$na.omit <- TRUE

#mosaic the rasters together
m <- do.call(raster::mosaic, x)
m

#name the bands
names(m) <- c("band_58", "band_34", "band_19")

#setNAvals
m[m < 0] <- NA

#see what the mosaiced raster looks like 
#plotRGB(m, r=1,g=2,b=3, stretch = "lin")

#for some reason writing to '.tif' messes up the file 
writeRaster(m, "./RGB/mosaic/brdfNWmosaic.grd", format="raster", overwrite=TRUE)






brdfrast <- brick("./RGB/comparison/SRERmerge_brdf.tif")
plotRGB(brdfrast,
        r=1,g=2,b=3,
        stretch = "lin")

origrast <- brick("./RGB/comparison/SRERmerge_orig.tif")
plotRGB(origrast,
        r=1,g=2,b=3,
        stretch = "lin")




plots <- readOGR("Z:/SRER/Martha/arcGIS/plots", "plots5x5")
plots$plot <- plots$siteID
plots$siteID <- rep("SRER", 50)
writeOGR(plots, "Z:/SRER/Martha/arcGIS/plots", "plots", driver="ESRI Shapefile")
???