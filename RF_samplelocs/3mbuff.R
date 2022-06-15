#################################################################################
###########create 3m rectangle buffer around plant sample locations##############
#################################################################################
#####code adopted from: https://www.neonscience.org/field-data-polygons-centroids


library(sp)
library(rgdal)

setwd("Z:/SRER/Martha/RFdata")
#load soil sample locations
soillocs <- readOGR(".", "allsamples")
samples <- as.data.frame(soillocs)

#set the radius for the plots
radius <- 1.5
#define the plot edges based upon the plot radius. 
##"coords.x1" is the adjEasting coordinate
##"coords.x2" is the adjNorthing coordinate
yplus <- samples$coords.x2 + radius
xplus <- samples$coords.x1 + radius
yminus <- samples$coords.x2 - radius
xminus <- samples$coords.x1 - radius

# calculate polygon coordinates for each plot centroid. 
square=cbind(xminus,yplus,  # NW corner
             xplus, yplus,  # NE corner
             xplus,yminus,  # SE corner
             xminus,yminus, # SW corner
             xminus,yplus)  # NW corner again - close ploygon

# Extract the sample ID information
ID=samples$sampleID
# make siteID vector
siteID <- rep("SRER", times=nrow(samples))


# create spatial polygons from coordinates, set the crs to "soillocs" spatialpoints dataframe 
polys <- SpatialPolygons(mapply(function(poly, id) 
{
  xy <- matrix(poly, ncol=2, byrow=TRUE)
  Polygons(list(Polygon(xy)), ID=id)
}, 
split(square, row(square)), ID),
proj4string=crs(soillocs))


# Create SpatialPolygonDataFrame -- this step is required to output multiple polygons.
polys.df <- SpatialPolygonsDataFrame(polys, data.frame(id=ID, row.names=ID))

#add siteID to the polygon spatial dataframe whenever "samplID" matches
polys.df$siteID <- siteID


polys2 <- as.data.frame(polys.df)

# write the shapefiles 
writeOGR(polys.df, '.', '3mbuff', 'ESRI Shapefile', overwrite_layer = TRUE)



