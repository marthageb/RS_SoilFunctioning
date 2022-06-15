getwd()
setwd("Z:/SRER/Martha/RFdata")
library(raster)
library(sp) # used to create a SpatialPoint object
library(rgdal)
library(sf)
library(maptools)
library(rgeos)

###FIRST NEED TO CLIP/EXTEND EACH RASTER SO IT ALLIGNS WITH THE FOLIAR N DATA PRODUCT

#load raters
setwd("Z:/SRER/Martha/RFdata")

###VEGETATION PARAMETERS
#FoliarN
r1 <- raster("Z:/SRER/Martha/hyperspec/flightlines/FoliarN/FoliarN.tif")
#NEON's aboveground biomass prediction
r2 <- raster("C:/Users/marthag/Documents/ArcGIS/Projects/SRER_AGB/agbALL.tif")
#r2 <- raster("Z:/SRER/Martha/RFdata/same_extent/NEON_AGB.tif")
#canopy ht derived from NEON LiDAR data
r3 <- raster("Z:/SRER/NEONAOP/NEON_2017_Data/L3/CanopyModel2017Mosaic.tif") #Kyle's CHM
#r3 <- raster("Z:/SRER/Martha/RFdata/same_extent/CHM.tif")

###TOPOGRAPHIC PARAMETERS
#elevation
r4 <- raster("Z:/SRER/Martha/arcGIS/DSMmosaic/allSRER.tif")
#r4 <- raster("Z:/SRER/Martha/RFdata/same_extent/elevation.tif")
#slope
r5 <- raster("C:/Users/marthag/Documents/ArcGIS/Projects/SRERslope/slopeALL.tif")  
#r5 <- raster("Z:/SRER/Martha/RFdata/same_extent/slope.tif")
#aspect
r6 <- raster("C:/Users/marthag/Documents/ArcGIS/Projects/SRERaspect/aspectALL.tif")
#r6 <- raster("Z:/SRER/Martha/RFdata/same_extent/aspect.tif")
#landform
r7 <- raster("Z:/SRER/Martha/RFdata/same_extent/landform.tif")
  
###CLIMATE
#1yr precip
r8 <- raster("C:/Users/marthag/Documents/ArcGIS/Projects/SRERprecip/c1yrprecip.tif")
#8mo precip
r9 <- raster("C:/Users/marthag/Documents/ArcGIS/Projects/SRERprecip/c8moprecip.tif") 
#3mo precip
r10 <- raster("C:/Users/marthag/Documents/ArcGIS/Projects/SRERprecip/c3moprecip.tif") #3mo precip

###SOIL TEXTURE/TYPE
#sand % from POLARIS
r11 <- raster("Z:/SRER/Martha/arcGIS/POLARIS/sandavg.tif")
#silt % from POLARIS
r12 <- raster("Z:/SRER/Martha/arcGIS/POLARIS/siltavg.tif")
#clay % from POLARIS
r13 <- raster("Z:/SRER/Martha/arcGIS/POLARIS/clayavg.tif")
#bulk density from POLARIS
r14 <- raster("Z:/SRER/Martha/arcGIS/POLARIS/BDavg.tif")
#rock type/parent material
r15 <- raster("Z:/SRER/Martha/RFdata/same_extent/PM.tif")
#soil age
r16 <- raster("Z:/SRER/Martha/RFdata/same_extent/soilage.tif")

#test raster layers to see if they 'match up'
#https://gis.stackexchange.com/questions/217082/handling-multiple-extent-problem-to-create-raster-stack-in-r
files <- list()
files[[1]] <- r1
files[[2]] <- r2
files[[3]] <- r3
files[[4]] <- r4
files[[5]] <- r5
files[[6]] <- r6
files[[7]] <- r7
files[[8]] <- r8
files[[9]] <- r9
files[[10]] <- r10
files[[11]] <- r11
files[[12]] <- r12
files[[13]] <- r13
files[[14]] <- r14
files[[15]] <- r15
files[[16]] <- r16

results <- list()

###test rasters
library(testthat)

test <- list()

for(i in 1:length(files)){
  test[[i]] <- capture_warnings(compareRaster(r7,files[[i]], res=T, orig=T, stopiffalse=F, showwarning=T))
}

test


##need to rebuild the CHM raster because when the "crop" function is applied (line 64) there is an offset in the extent and it causes stack/mask issue
#https://stackoverflow.com/questions/53440229/make-raster-stack-with-different-extent/53440900#53440900
#r7 -- landform used as the 'master raster'

#rc <- crop(r3,extent(r7))
#r3.1 <- raster(vals=values(rc), ext=extent(r7), crs=crs(r7), nrows=dim(r7)[1], ncols=dim(r7)[2])

##the cell resolution in the precipitation datasets is larger than in the foliarN dataset, so need to resample
r8.1 <- resample(r8,r7)
r9.1 <- resample(r9, r7)
r10.1 <- resample(r10, r7)
#r8.1 <- raster("Z:/SRER/Martha/RFdata/same_extent/1yrprecip.tif")
#r9.1 <- raster("Z:/SRER/Martha/RFdata/same_extent/8moprecip.tif")
#r10.1 <- raster("Z:/SRER/Martha/RFdata/same_extent/3moprecip.tif")

#reproject the soil texture rasters so that resolution, extent, and crs are matching
beginCluster()
r1.1 <- projectRaster(r1, r7)
endCluster()

beginCluster()
r11.1 <-projectRaster(r11, r7)
endCluster()

beginCluster()
r12.1 <- projectRaster(r12, r7)
endCluster()

beginCluster()
r13.1 <- projectRaster(r13, r7)
endCluster()

beginCluster()
r14.1 <- projectRaster(r14, r7)
endCluster()

#https://gis.stackexchange.com/questions/217082/handling-multiple-extent-problem-to-create-raster-stack-in-r
files <- list()
files[[1]] <- r1.1
files[[2]] <- r2
files[[3]] <- r3
files[[4]] <- r4
files[[5]] <- r5
files[[6]] <- r6
files[[7]] <- r7
files[[8]] <- r8.1
files[[9]] <- r9.1
files[[10]] <- r10.1
files[[11]] <- r11.1
files[[12]] <- r12.1
files[[13]] <- r13.1
files[[14]] <- r14.1
files[[15]] <- r15
files[[16]] <- r16

results <- list()

  ###test rasters
  library(testthat)
  
  test <- list()
  
  for(i in 1:length(files)){
    test[[i]] <- capture_warnings(compareRaster(r7,files[[i]], res=T, orig=T, stopiffalse=F, showwarning=T))
  }
  
  test

##Define function and run to extent/crop raster layers to foliar N layer
beginCluster()
for(i in 1:length(files)) {
  e <- extent(r7)
  r <-files[[i]] # raster(files[i])
  print(names(r))
  rc <- crop(r, e)
  if(sum(as.matrix(extent(rc))!=as.matrix(e)) == 0){ # edited
    rc <- mask(rc, r7) # You can't mask with extent, only with a Raster layer, RStack or RBrick
  }else{
    rc <- extend(rc,r7)
    rc<- mask(rc, r7)
  }
  
  # commented for reproducible example      
  results[[i]] <- rc # rw <- writeRaster(rc, outfiles[i], overwrite=TRUE)
  # print(outfiles[i])
  
}
endCluster()


env_data<- stack(results)

plot(env_data)
plot(r5)

#####
#writeout each raster layer that has been edited to have the same extent as the foliar N layer
sub1 <- subset(env_data,subset=14) #subset=raster layer number
sub1
writeRaster(sub1, "./same_extent/BD.tif", overwrite=TRUE)
#####

############soil texture polygon conversions###############
#load soil texture polygon feature class (only has map unit key for each soil type... no texture info)
text <- st_read("soiltexture.shp")
text <- as(text, "Spatial")
#read soil texture csv file (has soil texture info)
info <- read.csv("SRERtexture.csv")
#combine sand, silt, clay, BD info with MUKEY polys
texture <- merge(text, info,  by='MapUnitKey')
texture
#load raster you want the "fit" this one to
r5 <- raster("./same_extent/NEON_AGB.tif") #NEON_AGB

#create blank raster with dimensions and extents same as FoliarN raster
r <- raster()
dim(r) <- dim(r5)
extent(r) <- extent(r5)
crs(r) <- crs(r5)

#rasterize soil texture polygon featureclasses
mukey <- rasterize(texture, r, 'MUKEY') #soil map unit key raster
sand <- rasterize(texture, r, 'per_sand')
silt <- rasterize(texture, r, 'per_silt')
clay <- rasterize(texture, r, 'per_clay')
bulk_den <- rasterize(texture, r, 'BD')

#double check to make sure it is the same extent, orgin, resolution, and crs as the Foliar N raster 
capture_warnings(compareRaster(clay,r5, res=T, orig=T, stopiffalse=F, showwarning=T))

writeRaster(mukey, "./same_extent/soil_mukey.tif")
writeRaster(sand, "./same_extent/per_sand.tif", overwrite=TRUE)
writeRaster(silt, "./same_extent/per_silt.tif", overwrite=TRUE)
writeRaster(clay, "./same_extent/per_clay.tif", overwrite=TRUE)
writeRaster(bulk_den, "./same_extent/bulk_den.tif", overwrite=TRUE)


