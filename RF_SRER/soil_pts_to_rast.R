getwd()
setwd("D:/martha/RFdata")

#load .shp file with soil sample location and corresponding lab measurements
samples <- readOGR("soildata.shp")

#load a raster dataset to build a blank one from. (one you want to "fit" the soil samples to)
mast <- raster("./same_extent/CHM.tif")

#create blank raster with dimensions and extents same as FoliarN raster
r <- raster()
dim(r) <- dim(mast)
extent(r) <- extent(mast)
crs(r) <- crs(mast)

#turn soil sample .shp into a raster based on values for one lab measurement
rast <- rasterize(samples, r, field="pH", fun=mean, filename="./train_pts/ph.tif" )
