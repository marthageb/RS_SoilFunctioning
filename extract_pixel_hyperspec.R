# Load required packages
library(raster)
library(rhdf5)

####CREATE A NEW SHAPE FILE FOR SAMPLE POINTS
points <- read.csv("C:/Users/marthag/Desktop/envi_spectralJan2020/xycoordsNEW.csv")

coordinates(points)= ~ xcoord + ycoord

pointsold <- readOGR("Z:/SRER/Martha/arcGIS/plots/ENVIpts.shp")

proj4string(points) <- proj4string(pointsold)

writeOGR(points, "Z:/SRER/Martha/arcGIS/plots", "ENVIptsNEW", driver="ESRI Shapefile", overwrite_layer = TRUE)


####EXTRACT PIXEL HYPERSPEC READINGS
#load ENVI point coords
points <- read.csv("C:/Users/marthag/Desktop/envi_spectralJan2020/xycoordsNEW.csv")
#points <- subset(points, type %in% c("grass", "mesq", "cres"))

setwd("Z:/SRER/Martha/hyperspec/flightlines/output/noNDVI_thick_sparse/")
refl.path <- "BRDF/Correction"

#define the hdf file and plot ID
    #site=="A1" & raster =="NEON_D14_SRER_DP1_20170829_180310_reflectance.h5")
    #site=="A2" & raster =="NEON_D14_SRER_DP1_20170825_165904_reflectance.h5")
    #site=="B1" & raster =="NEON_D14_SRER_DP1_20170824_202731_reflectance.h5")
    #site=="B2" & raster =="NEON_D14_SRER_DP1_20170824_195807_reflectance.h5")
    #site=="C1" & raster =="NEON_D14_SRER_DP1_20170829_161226_reflectance.h5")
    #site %in% c("C2", "E1", "G2") & raster =="NEON_D14_SRER_DP1_20170824_161655_reflectance.h5")
    #site=="D1" & raster =="NEON_D14_SRER_DP1_20170824_202008_reflectance.h5")
    #site %in% c("D2", "G1") & raster =="NEON_D14_SRER_DP1_20170824_170703_reflectance.h5")
    #site=="E2" & raster =="NEON_D14_SRER_DP1_20170824_162407_reflectance.h5")
    #site=="E3" & raster =="NEON_D14_SRER_DP1_20170825_165121_reflectance.h5")
    #site=="F1" & raster =="NEON_D14_SRER_DP1_20170824_163106_reflectance.h5")
    #site=="F2" & raster =="NEON_D14_SRER_DP1_20170824_174415_reflectance.h5")



f <- "NEON_D14_SRER_DP1_20170824_161655_reflectance.h5"
xy <- subset(points, site == "G2")


#define the function to extract reflectance readings at each pixel
extract_pixel <- function(z){
  colid <- as.numeric(z['col'])
  rowid <- as.numeric(z['row'])
  pixel <- h5read(f, refl.path, index = list(1:426, colid, rowid))
  return(pixel)
  h5closeAll()
}

#apply the "extract_pixel" function over each row in 'xy'
vals <- apply(xy,1,extract_pixel)

vals_t <- t(vals)
vals_t <- as.data.frame(vals_t)
#pull wavelength values from the .h5 file and make them the column names
wavelengths <- h5read(f,"/SRER/Reflectance/Metadata/Spectral_Data/Wavelength")
colnames(vals_t) <- wavelengths
final <- cbind.data.frame(xy[,c(2,5:8)], vals_t)

#combine all of the sites together
final_all <- final #for the first one 
#final_all <- rbind.data.frame(final_all, final)

##################################################################################
#################after getting readings extracted from all sites##################
##################################################################################

#replace -10000 values with 'NA'
final_all[final_all == -10000] <- NA

#divide by the scale value and get rid of NA observations
#remove rows with NA readings
ref <- final_all[,10:422]
noNA <- ref[complete.cases(ref),]
names <- rownames(noNA)
data <- subset(final_all, rownames(final_all) %in% names)
#divide by scale factor
data[,6:431] <- data[,6:431]/10000

#make plotID column
plotID <- substr(data$ENVI_ID, 1, 6)
data <- cbind.data.frame(plotID, data)


#load foliar chem data
chem <- read.csv("Z:/SRER/Martha/hyperspec/flightlines/FoliarN/plotlvlN.csv")
chem <- chem[,-c(1,5:7)]
chem$plotID <- substr(chem$plotID, 1,6)

newsite <- merge(chem, data, by="plotID", all.y=TRUE)
newsite$plotN[is.na(newsite$plotN)] <- 0
plants <- rbind.data.frame(plants, newsite)
#plants <- plants[!(plants$plotID %in% c("E2bare", "E2gras", "E2mesq", "E2cres")),] #remove certain observations

#write.csv(plants, "Z:/SRER/Martha/hyperspec/flightlines/FoliarN/NDVI0.2_TS_envi.csv")
