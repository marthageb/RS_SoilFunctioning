library(plyr)

setwd("Z:/SRER/Martha/RFdata/extracted_hyperspec")

#load extracted reflectance data (output from extract_hyperspec.py)
ref <- ldply(.data = list.files("Z:/SRER/Martha/RFdata/extracted_hyperspec"),
             .fun = read.csv)
  
#make site column
ref$site <- substr(ref$id, 1,2)

#subset the readings for each site for the specific flightline
A1 <- subset(ref, site=="A1" & raster =="NEON_D14_SRER_DP1_20170829_180310_reflectance.h5")
A2 <- subset(ref, site=="A2" & raster =="NEON_D14_SRER_DP1_20170825_165904_reflectance.h5")
B1 <- subset(ref, site=="B1" & raster =="NEON_D14_SRER_DP1_20170824_202731_reflectance.h5")
B2 <- subset(ref, site=="B2" & raster =="NEON_D14_SRER_DP1_20170824_195807_reflectance.h5")
C1 <- subset(ref, site=="C1" & raster =="NEON_D14_SRER_DP1_20170829_161226_reflectance.h5")
C2 <- subset(ref, site %in% c("C2", "E1", "G2") & raster =="NEON_D14_SRER_DP1_20170824_161655_reflectance.h5")
D1 <- subset(ref, site=="D1" & raster =="NEON_D14_SRER_DP1_20170824_202008_reflectance.h5")
D2 <- subset(ref, site %in% c("D2", "G1") & raster =="NEON_D14_SRER_DP1_20170824_170703_reflectance.h5")
E2 <- subset(ref, site=="E2" & raster =="NEON_D14_SRER_DP1_20170824_162407_reflectance.h5")
E3 <- subset(ref, site=="E3" & raster =="NEON_D14_SRER_DP1_20170825_165121_reflectance.h5")
F1 <- subset(ref, site=="F1" & raster =="NEON_D14_SRER_DP1_20170824_163106_reflectance.h5")
F2 <- subset(ref, site=="F2" & raster =="NEON_D14_SRER_DP1_20170824_174415_reflectance.h5")

all <- rbind.data.frame(A1, A2, B1, B2, C1, C2, D1, D2, E2, E3, F1, F2)

#apply the scale factor
all[,5:430] <- all[,5:430]/10000

#reorder columns
all <- all[c(431:433,1:4,5:430)]


#aggregate -- calculate mean for each plot
plot <- aggregate(all, by=list(all$id), FUN=mean)
names(plot)[names(plot) == "Group.1"] <- "sampleID"
plot <- plot[,-c(4:8)]
#write.csv
write.csv(plot, "Z:/SRER/Martha/RFdata/3mbuff/ref.csv")
