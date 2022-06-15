library(data.table)
library(dplyr)
#######calculate sample-averaged hyperspec reflectance and output file to "3mbuff" folder###########
setwd("Z:/SRER/Martha/hyperspec/flightlines/extract_hyperspec/NDVI_0.2_thick_sparse_7x7plots/")

#create a list all of the file names
files <- list.files(pattern = ".csv$")

#define a function to remove all observations that have no reflectance readings
ref_vals <- function(x){
  data <- read.csv(x)
  #replace -9999 values with 'NA'
  data[data == -10000] <- NA
  
  #remove rows with NA readings
  ref <- data[,15:439]
  noNA <- ref[complete.cases(ref),]
  names <- rownames(noNA)
  data1 <- subset(data, rownames(data) %in% names)
  
  #create a column of the plant/plot type
  type <- substr(data1$plot, 3,6)
  data1 <- cbind.data.frame(type,data1)
  return(data1)
}

#apply the "ref_val" function over the csv files
ref_tmp <- lapply(files, ref_vals)
#Bind rows together
reflect <- do.call(rbind, ref_tmp)

#select only the creseote, grass, and mesquite plots
plants <- subset(reflect, type %in% c("gras", "mesq", "cres"))
plants <- plants[,-c(1:10)]
plants <- plants[,c(429,430,1:428)]

site <- substr(plants$plot, 1,2)
plants <- cbind.data.frame(site,plants)

A1 <- subset(plants, site=="A1" & raster =="NEON_D14_SRER_DP1_20170829_180310_reflectance.h5")
A2 <- subset(plants, site=="A2" & raster =="NEON_D14_SRER_DP1_20170825_165904_reflectance.h5")
B1 <- subset(plants, site=="B1" & raster =="NEON_D14_SRER_DP1_20170824_202731_reflectance.h5")
B2 <- subset(plants, site=="B2" & raster =="NEON_D14_SRER_DP1_20170824_195807_reflectance.h5")
C1 <- subset(plants, site=="C1" & raster =="NEON_D14_SRER_DP1_20170829_161226_reflectance.h5")
C2 <- subset(plants, site %in% c("C2", "E1", "G2") & raster =="NEON_D14_SRER_DP1_20170824_161655_reflectance.h5")
D1 <- subset(plants, site=="D1" & raster =="NEON_D14_SRER_DP1_20170824_202008_reflectance.h5")
D2 <- subset(plants, site %in% c("D2", "G1") & raster =="NEON_D14_SRER_DP1_20170824_170703_reflectance.h5")
E2 <- subset(plants, site=="E2" & raster =="NEON_D14_SRER_DP1_20170824_162407_reflectance.h5")
E3 <- subset(plants, site=="E3" & raster =="NEON_D14_SRER_DP1_20170825_165121_reflectance.h5")
F1 <- subset(plants, site=="F1" & raster =="NEON_D14_SRER_DP1_20170824_163106_reflectance.h5")
F2 <- subset(plants, site=="F2" & raster =="NEON_D14_SRER_DP1_20170824_174415_reflectance.h5")

all <- rbind.data.frame(A1, A2, B1, B2, C1, C2, D1, D2, E2, E3, F1, F2)

write.csv(plants, 'Z:/SRER/Martha/hyperspec/flightlines/FoliarN/NDVI0.2_thick_sparse.csv')



###############################################################################
###############combine reflectance vals with foliar chem data vars#############
###############################################################################
#load the descriptive variables
chem <- read.csv("Z:/SRER/Martha/hyperspec/flightlines/FoliarN/plotlvlN.csv")
chem <- chem[,-c(1,5:7)]



data <- merge(chem, plants, by="plotID")
data <- data[,-c(4,7)]
#divide hyperspectral readings by 10,000 to get actual reflectance readings
refl <- data[,c(6:431)] /10000
desc <- data[,c(1:5)]

final <- cbind.data.frame(desc,refl)

write.csv(final, "Z:/SRER/Martha/hyperspec/flightlines/FoliarN/NDVI0.2_thick_sparse.csv")



xycoords <- read.csv("C:/Users/marthag/Desktop/envi_spectralJan2020/xycoords.csv")
xy <- xycoords[,c(2,7,8)]
colnames(xy) <- c("ENVI_ID", "x", "y")
new <- match_df(final,xy)
new$plot <- as.factor(new$plot)
summary(new$plot)

new <- subset(new, !plot =="D1mesq")
new <- subset(new, !plot =="E3cres")

d1 <- subset(xycoords, site=="D1" & type=="mesq")
xy <- d1[,c(2,7,8)]
colnames(xy) <- c("ENVI_ID", "x", "y")
newd1 <- match_df(final,xy)

new <- rbind.data.frame(new, newd1)

cres <- subset(final, type=='cres')


new <- rbind.data.frame(new,cres)
