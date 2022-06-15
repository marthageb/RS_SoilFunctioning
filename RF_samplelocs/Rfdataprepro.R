setwd("Z:/SRER/Martha/RFdata")
####################################################################################################
############single pixel processing first, 3x3m buffer processing starting at line 78###############
####################################################################################################
###SINGLE PIXEL###
#load data where raster values were extracts from soil sample locs (output from "raster value extraction.R")
rasdata <- read.csv("rfdata_allsamples.csv")
#load hyperspectral reflectance readings
hyper <- read.csv("./ref/ref_all.csv")
#load predicted foliar N values (using the PLSR model coefficients)
foliarN <- read.csv("PLSRpredN.csv")

#change sampleIDS to having the 1st letter capitilized (to match with sampleIDs on rasdata)
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

hyper$IDs <- sapply(hyper$sampleID, CapStr)
foliarN$IDs <- sapply(foliarN$sampleID, CapStr)

#change all negative foliar N values to 0
foliarN$Nnoneg <- foliarN$Npred
foliarN$Nnoneg[foliarN$Nnoneg<0] <- 0

#there isn't SSURGO soil data for "C" sites. So that they can be included in the RF I will use averaged texture values measured in the lab
text <- read.csv("soiltext.csv")
text$siteID <- as.factor(text$siteID)
text$land <- text$geomorph
aov1 <- aov(clay~siteID, text)
summary(aov1)
print(model.tables(aov1, "means", se=T))
###assign average sand, silt, clay, BD vals to C1 and C2 samples
  #C2
rasdata[c(1,16,17,41,61,80),17] <- 81.52
rasdata[c(1,16,17,41,61,80),18] <- 14.74
rasdata[c(1,16,17,41,61,80),19] <- 3.743
rasdata[c(1,16,17,41,61,80),20] <- 1.48
  #C1
rasdata[c(14,15,53,54,60,77,78,79),17] <- 69.15
rasdata[c(14,15,53,54,60,77,78,79),18] <- 26.02
rasdata[c(14,15,53,54,60,77,78,79),19] <- 4.824
rasdata[c(14,15,53,54,60,77,78,79),20] <- 1.52


#subset to include only wanted columns
rasdata1 <- rasdata[,c(2:4,7:20)]
hyper1 <- hyper[,c(429,3:428)]
names(hyper1)[names(hyper1) == "IDs"] <- "sampleID"
foliarN1 <- foliarN[,5:6]
names(foliarN1)[names(foliarN1) == "Nnoneg"] <- "foliarN"
names(foliarN1)[names(foliarN1) == "IDs"] <- "sampleID"
desc <- text[,c(4,2,5,6,8,17,11,12)] #descriptive variables from soil texture dataset
names(desc)[names(desc) == "IDs"] <- "sampleID"

all <- merge(rasdata1,desc, by="sampleID", all.y=FALSE)
all <- all[!duplicated(all$sampleID),]#delete the extra grass samples
all <- merge(all, foliarN1, by="sampleID")
all <- merge(all,hyper1, by="sampleID")

#######################
#load soil lab measurements and aggregate the grass soil samples by taking the average
soil <- read.csv("soildata.csv")
IDs <- text[,c(2,4)]
soil <- merge(IDs,soil, by="uniqueID")
soil1 <- soil[,c(2,9:31)]

agg <- aggregate(soil1, by=list(soil1$IDs), FUN=mean)
agg <- agg[,-2]
names(agg)[names(agg) == "Group.1"] <- "sampleID"

final <- merge(agg, all, by="sampleID")

write.csv(final, "rfdata.csv")

#####################################################
####################3x3m buffers#####################
#####################################################
#load data where raster values were extracts from soil sample locs (output from "raster value extraction.R")
rasdata <- read.csv("./3mbuff/rasterdata.csv")
#load hyperspectral reflectance readings
hyper <- read.csv("./3mbuff/ref.csv")

#get rid of extra columns and rename sampleID
rasdata <- rasdata[,-1]
hyper <- hyper[,-1]


all <- merge(rasdata,hyper[,-c(2:3)], by="sampleID", all.y=FALSE)

#######################
#load soil lab measurements and aggregate the grass soil samples by taking the average
soil <- read.csv("soildata.csv")

#create sampleID column to match the pother datasets
soil$mesq_type[soil$mesq_type == "bole"] <- "mesq"
soil$rep <- substr(soil$sampleID, 7,9)
soil$rep[soil$rep == "s"] <- ""
soil$rep[soil$rep == "1-1"] <- ""
soil$rep[soil$rep == "s1-"] <- ""
soil$rep[soil$rep == "s2-"] <- ""

soil$sampleID <- paste(soil$siteID, soil$mesq_type, soil$rep, sep="")

#remove unwanted columns
soil <- soil[,-c(1,4:7,31)]


#calcualte mean values for grass samples
agg <- aggregate(soil, by=list(soil$sampleID), FUN=mean)
agg <- agg[,-c(2:3)]
names(agg)[names(agg) == "Group.1"] <- "sampleID"

final <- merge(agg, all, by="sampleID")

write.csv(final, "./3mbuff/rfdata.csv")



#############################################################################################
###########calculate soil functioning index##################################################
#############################################################################################
#calculate log normalized variables
data$lnDOC <- log(data$DOC)
data$lnMBC <- log(data$MBC + 1)
data$lnDON <- log(data$DON + 1)
data$lnMBN <- log(data$MBN + 1)
data$lndna <- log(data$dna)
data$lnOM <- log(data$OM. + 1)
data$lnnh4 <- log(data$NH4i + 1)
data$lnno3 <- log(data$NO3i +1)

data$lnag <- log(data$AG+1)
data$lnbg <- log(data$BG+1)
data$lncb <- log(data$CB+1)
data$lnxyl <- log(data$XYL+1)
data$lnnag <- log(data$NAG+1)
data$lnlap <- log(data$LAP)
data$lnphos <- log(data$PHOS)

#soil nutrient and multifunctionality indices w zcore transformations
data$z_om <- scale(data$lnOM, center=TRUE, scale=TRUE)
data$z_doc <- scale(data$lnDOC, center=TRUE, scale=TRUE)
data$z_don <- scale(data$lnDON, center=TRUE, scale=TRUE)
data$z_nh4 <- scale(data$lnnh4, center=TRUE, scale=TRUE)
data$z_no3 <- scale(data$lnno3, center=TRUE, scale=TRUE)

data$z_ag <- scale(data$lnag, center=TRUE, scale=TRUE)
data$z_bg <- scale(data$lnbg, center=TRUE, scale=TRUE)
data$z_cb <- scale(data$lncb, center=TRUE, scale=TRUE)
data$z_xyl <- scale(data$lnxyl, center=TRUE, scale=TRUE)
data$z_nag <- scale(data$lnnag, center=TRUE, scale=TRUE)
data$z_lap <- scale(data$lnlap, center=TRUE, scale=TRUE)
data$z_phos <- scale(data$lnphos, center=TRUE, scale=TRUE)

data$z_dna <- scale(data$lndna, center=TRUE, scale=TRUE)
data$z_mbc <- scale(data$lnMBC, center=TRUE, scale=TRUE)
data$z_mbn <- scale(data$lnMBN, center=TRUE, scale=TRUE)

#indices
data$mfi <- (data$z_ag + data$z_bg + data$z_cb + data$z_xyl + data$z_nag + data$z_lap + data$z_phos + data$z_om + data$z_doc + data$z_don + data$z_nh4 + data$z_no3 + data$z_dna + data$z_mbn + data$z_mbc ) / 15
