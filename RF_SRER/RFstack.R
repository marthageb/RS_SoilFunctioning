library(testthat)
library(raster)
library(gdalUtils)
library(rhdf5)
library(rgdal)
library(raster)
library(maps)
library(parallel)
library(lme4)


#############################################################
#########Extracting single bands from NEON hdf files#########
#############################################################
#setwd("Z:/SRER/Martha/RFdata")
alltiles <- list.files("D:/martha/hyperspec/flightlines/output/noNDVI_thick_sparse", pattern="h5", full.names = TRUE)
#define the crs
myCRS <- as.character("+proj=utm +zone=12 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
#get wavelengths
f <- "D:/martha/hyperspec/a1_20170829_180310_reflectance.h5"
wavelengths<- h5read(f,"/SRER/Reflectance/Metadata/Spectral_Data/Wavelength")
#convert to matrix
wavelengths <- as.matrix((wavelengths))
#figure out which index matches the band of interest, then define it
b_int <- 31




#function to extract a single band from each mosaiced file
bandextraction <- function(x){
  #load results.hdf5 file where each tile is a seperate result
  tile1 <- x
  #get file name by extracting characters after "thick_sparse"
  name <- sub(".*thick_sparse","", tile1)
  print(name)
  og <- paste("D:/martha/hyperspec/flightlines", name, sep="")
  #look at hdf5 file format
  #h5ls(tile1,all=T) #there are 9 sub folders in results each with about 50 tiles
  # get spatialInfo using the h5readAttributes function 
  spInfo <- h5readAttributes(tile1,"/SRER/Reflectance/Metadata/Coordinate_System/Map_Info")
  # get attributes for the Reflectance dataset
  reflInfo <- h5readAttributes(og,"/SRER/Reflectance/Reflectance_Data") #if using the raw/original NEON data
  
  #get the wavelenght vals
  wavelengths<- h5read(tile1,"/SRER/Reflectance/Metadata/Spectral_Data/Wavelength")
  
  nRows <- reflInfo$Dimensions[1]
  nCols <- reflInfo$Dimensions[2]
  nBands <- reflInfo$Dimensions[3]
  
  
  oneband<- h5read(tile1,"/BRDF/Correction",index=list(18, 1:nCols,1:nRows))
  # what type of object is b34?
  #class(oneband)
  
  # convert from array to matrix
  oneband <- oneband[1,,]
  # check it
  #class(oneband)
  #divide by scalefator
  oneband <- oneband / 10000
  # there is a no data value in our raster - let's define it
  #myNoDataValue <- as.numeric(reflInfo$`Data_Ignore_Value`)
  myNoDataValue <- -1
  
  # set all values = -9999 to NA
  oneband[oneband == myNoDataValue] <- NA
  
  # We need to transpose x and y values in order for our 
  # final image to plot properly
  oneband<-t(oneband)
  
  
  # Populate the raster image extent value. 
  # get the map info, split out elements
  mapInfo<-h5read(tile1,"/SRER/Reflectance/Metadata/Coordinate_System/Map_Info")
  
  # Extract each element of the map info information 
  # so we can extract the lower left hand corner coordinates.
  mapInfo<-unlist(strsplit(mapInfo, ","))
  
  # view the attributes in the map dataset
  #mapInfo
  
  # Create the projection in as object
  #crsinfo <- h5read(tile1,"/SRER/Reflectance/Metadata/Coordinate_System")
  #myCRS <- crsinfo$Proj4
  #myCRS
  
  #str(myCRS)
  
  #myCRS <- as.character("+proj=utm +zone=12 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  
  # define final raster with projection info 
  # if UTM is all caps it might cause an error!
  onebandr <- raster(oneband, crs=myCRS)
  
  onebandr
  
  # grab resolution of raster as an object
  res <- as.numeric(mapInfo[2])
  
  # Grab the UTM coordinates of the upper left hand corner of the raster
  
  #grab the left side x coordinate (xMin)
  xMin <- as.numeric(mapInfo[4]) 
  #grab the top corner coordinate (yMax)
  yMax <- as.numeric(mapInfo[5])
  
  
  # Calculate the lower right hand corner to define the full extent of the 
  # raster. To do this we need the number of columns and rows in the raster
  # and the resolution of the raster.
  
  # note that you need to multiple the columns and rows by the resolution of 
  # the data to calculate the proper extent!
  xMax <- (xMin + (ncol(oneband))*res)
  yMin <- (yMax - (nrow(oneband))*res) 
  
  #xMax
  #yMin
  
  # define the extent (left, right, top, bottom)
  rasExt <- extent(xMin,xMax,yMin,yMax)
  rasExt
  
  # assign the spatial extent to the raster
  extent(onebandr) <- rasExt
  
  # look at raster attributes
  onebandr
  
  #close the h5files
  #h5closeAll()
  
}

#Lapply "rast_files" function over full list of files: takes less than a minute for each subfolder
#rast_tmp <- lapply(alltiles, bandextraction)


#do the lapply function on parallel processors to speed up time
numCores <- detectCores()
cl <- makeCluster(numCores)
#load package needed for each process
clusterEvalQ(cl, {
  library(rhdf5)
  library(rgdal)
  library(raster)
  myCRS <- as.character("+proj=utm +zone=12 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
})

rast_tmp <- parLapply(cl, alltiles, bandextraction)

#best practice to avoid processes running in the background 
stopCluster(cl)

rast_tmp$fun <- mean

#mosaic all of the rasters together
mos <- do.call(raster::mosaic, rast_tmp)

getwd()
raster::writeRaster(mos,"D:/martha/RFdata/same_extent/unfit_ref_bands/b18.tif")


###################################################################################
############Get extracted reflectance band to 'fit' other raster layers############
###################################################################################
setwd("D:/martha/RFdata/")

##first get re'flectance bands to "fit" other raster layers##
mast <- raster("./same_extent/CHM.tif")
b1 <- raster("./same_extent/unfit_ref_bands/b67.tif")
b2 <- raster("./same_extent/unfit_ref_bands/b39.tif")
b3 <- raster("./same_extent/unfit_ref_bands/b200.tif")
b4 <- raster("./same_extent/unfit_ref_bands/b3.tif")
b5 <- raster("./same_extent/unfit_ref_bands/b423.tif")


b6 <- raster("./same_extent/unfit_ref_bands/b151.tif")
b7 <- raster("./same_extent/unfit_ref_bands/b150.tif")
b8 <- raster("./same_extent/unfit_ref_bands/b153.tif")


#https://gis.stackexchange.com/questions/217082/handling-multiple-extent-problem-to-create-raster-stack-in-r
files <- list()
files[[1]] <- b1
files[[2]] <- b2
files[[3]] <- b3
files[[4]] <- b4
files[[5]] <- b5
files[[6]] <- b6
files[[7]] <- b7
files[[8]] <- b8


files[[9]] <- 
files[[10]] <- 
files[[11]] <- 


###test rasters
test <- list()

for(i in 1:length(files)){
  test[[i]] <- capture_warnings(compareRaster(r4,files[[i]], res=T, orig=T, stopiffalse=F, showwarning=T))
}
test

##Define function and run to extend/crop raster layers to foliar N layer
results <- list()
for(i in 1:length(files)) {
  e <- extent(mast)
  r <-files[[i]] # raster(files[i])
  rc <- crop(r, e)
  if(sum(as.matrix(extent(rc))!=as.matrix(e)) == 0){ # edited
    rc <- mask(rc, mast) # You can't mask with extent, only with a Raster layer, RStack or RBrick
  }else{
    rc <- extend(rc,mast)
    rc<- mask(rc, mast)
  }
    results[[i]] <- rc
  # print(outfiles[i])
  
}

env_data <- stack(results[1:1])

#####
#writeout each raster layer that has been edited to have the same extent as the foliar N layer
sub1 <- subset(env_data,subset=1) #subset=raster layer number
writeRaster(sub1, "./same_extent/b67.tif")

###################################################################################
########################Create raster stack for RF#################################
###################################################################################

#####then, stack and write raster#########
r1 <- raster("./same_extent/soilage.tif") 
r2 <- raster("./same_extent/sand.tif")
r3 <- raster("./same_extent/FoliarN.tif")
r4 <- raster("./same_extent/b67.tif")
r5 <- raster("./same_extent/slope.tif") 
r6 <- raster("./same_extent/silt.tif") 
r7 <- raster("./same_extent/PM.tif") 
r8 <- raster("./same_extent/soilage.tif") 
r9 <- raster("./same_extent/1yrprecip.tif") 


r6 <- raster("./same_extent/b3.tif")
r7 <- raster("./same_extent/b63.tif") 


r8 <- raster("./same_extent/b303.tif")
r9 <- raster("./same_extent/b202.tif")
r10 <- raster("./same_extent/b63.tif")
r11 <- raster("./same_extent/b151.tif")
r12 <- raster("./same_extent/CHM.tif")
r13 <- raster("./same_extent/b150.tif")



#create a raster brick of all raster layers
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

rfdata<- stack(files)
#write raster stack
writeRaster(rfdata, "RFdataCFI.tif", overwrite=TRUE)


#####CREATE one HOT VARIABLES FOR CATEGORICAL VARS#####
library(RStoolbox)
#load raster layer w categorical vars
cats <- raster("./same_extent/PM.tif")
cats
sc <- unsuperClass(cats, nClasses = 4)

## one-hot encode 
sc_oneHot <- oneHotEncode(sc$map, classes = c(1,2,3,4))

## check results
sc_oneHot

getwd()
writeRaster(sc_onehot, "./same_extent/PM1hot.tif")
