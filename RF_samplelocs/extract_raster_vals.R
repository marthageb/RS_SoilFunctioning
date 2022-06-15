library(raster)
library(sp) # used to create a SpatialPoint object
library(rgdal)
library(sf)
library(maptools)

#############################
#load "environmental variables" (raster layers) used in RF
## all layers have been previously modified to have the same extent, orgin, resolution, and crs as Foliar N raster layer using "raster layer extents.R"
setwd("Z:/SRER/Martha/RFdata")

###VEGETATION PARAMETERS
#FoliarN
r1 <- raster("./same_extent/FoliarN.tif")
#NEON's aboveground biomass prediction
r2 <- raster("./same_extent/NEON_AGB.tif")
#canopy ht derived from NEON LiDAR data
r3 <- raster("./same_extent/CHM.tif")

###TOPOGRAPHIC PARAMETERS
#elevation
r4 <- raster("./same_extent/elevation.tif")
#slope
r5 <- raster("./same_extent/slope.tif")
#aspect
r6 <- raster("./same_extent/aspect.tif")
#landform
r7 <- raster("./same_extent/landform.tif")

###CLIMATE
#1yr precip
r8 <- raster("./same_extent/1yrprecip.tif")
#8mo precip
r9 <- raster("./same_extent/8moprecip.tif") 
#3mo precip
r10 <- raster("./same_extent/3moprecip.tif") #3mo precip

###SOIL TEXTURE/TYPE
#sand % from POLARIS
r11 <- raster("./same_extent/sand.tif")
#silt % from POLARIS
r12 <- raster("./same_extent/silt.tif")
#clay % from POLARIS
r13 <- raster("./same_extent/clay.tif")
#bulk density from POLARIS
r14 <- raster("./same_extent/BD.tif")
#rock type/parent material
r15 <- raster("./same_extent/PM.tif")
#soil age
r16 <- raster("./same_extent/soilage.tif")
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
files[[14]] <- r14
files[[15]] <- r15
files[[16]] <- r16

rfdata<- stack(files)
#write raster stack
writeRaster(rfdata, "RFdata.tif", overwrite=TRUE)

#load polygon plot (or sample location points) feature class
plots <- st_read("allsamples.shp") #single pixel
plots <- st_read("3mbuff_points.shp") #3x3m buffer area
#convert plots to Spatial* object (which will be needed to get x/y coords from the extract function)
plots2 <- as(plots, "Spatial")

#plot(r5)
#plot(plots, add=TRUE)

#extract values from raster stack based on plot locations
data1 <- raster::extract(rfdata, plots, df=TRUE, cellnumbers=T) 

#get coordinates for each of plot pixels
xandy <- as.data.frame(coordinates(rfdata)[data1[,2],])

#cbind coordinates to extracted raster stack values
data2 <- cbind(xandy, data1)
data2 <- cbind(plots$id,data2)
#data2 <- cbind(plots$uniqueID, data2)
names(data2)[names(data2) == "plots$id"] <- "sampleID"

######For 3x3 buffs will need to calculate mode/mean######
#mode function
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

landform <- aggregate(landform ~ sampleID, data2, FUN = Mode )
PM <- aggregate(PM ~ sampleID, data2, FUN = Mode)
soilage <- aggregate(soilage ~ sampleID, data2, FUN = Mode)

catvar <- merge(landform, PM, by="sampleID")
catvar <- merge(catvar, soilage, by="sampleID")

contvar <- data2[,c(1,6:9,11,10,13:19)]
contvar <- aggregate(contvar[,2:14], by=list(contvar$sampleID), FUN = mean, na.rm=TRUE)
names(contvar)[names(contvar) == "Group.1"] <- "sampleID"

finalrf <- merge(catvar, contvar, by="sampleID") #3x3 buffer area


#rename columns
#single pixels
colnames(data2) <- c("sampleID","x", "y", "ID", "cell", "geomorph", "vegclass","NEON_AGB", 
                     "precip1yr", "precip8mo", "precip3mo", "elevation", "aspect", "slope", 
                     "CHM", "sand", "silt", "clay", "BD")
#3x3 buffer area
finalrf <- finalrf[c("sampleID","FoliarN", "NEON_AGB", "CHM", "elevation", "aspect", "slope", "landform", "X1yrprecip", "X8moprecip", "X3moprecip", "sand", "silt", "clay", "BD", "PM", "soilage")]

write.csv(data2, file="rasterdata.csv") #single pixel
write.csv(finalrf, file="./3mbuff/rasterdata.csv") #3x3 buffer area






