library(randomForest)
library(rfPermute)
library(VSURF)
library(dplyr)

setwd("Z:/SRER/Martha/RFdata")
#single pixels
data <- read.csv("rfdata.csv")
#3x3mbuff
data <- read.csv("./3mbuff/rfdata.csv")

#delete duplicate columns
data <- data[,-c(1,24)]

#replace negative MBC values with NA
data$MBN[data$MBN < 0] <- NA
data$MBC[data$MBC < 0] <- NA

str(data)

#brightness normalization per Feilhauer et al 2010
ref <- data[,40:465]
ref <- as.matrix(ref)
#create the brightness normalized function
bnorm <- function (X)
{ X <- X / sqrt (rowSums (X ^ 2, na.rm=T))}
#apply the brightness normalized function
bnref <- bnorm(ref)
bnref <- as.data.frame(bnref)
#exclude the water absorption regions at 937 and 1127
x1 <- dplyr::select(bnref, `b_436.362396`:`b_922.151001`)
x2 <- dplyr::select(bnref, `b_962.216003`:`b_1107.451782`)
x3 <- dplyr::select(bnref, `b_1162.54126`:`b_1317.793213`)
x4 <- dplyr::select(bnref, `b_1488.06958`:`b_1768.524902`)
x5 <- dplyr::select(bnref, `b_2038.963867`:`b_2434.606201`)
x <- cbind(x1, x2, x3, x4, x5)



#combine brightness normalized values with descriptive variables
descvars <- data[,1:39]
data <- cbind(descvars, x)

str(data)

#change character variables into factors
data$landform <- as.factor(data$landform)
data$PM <- as.factor(data$PM)
data$soilage <- as.factor(data$soilage)


####CALCULATE COLLECTIVE SOIL FUNCTIONALITY INDEX####
summary(data$NO3i)
data$lnom <- log(data$OM.+1)
data$lndoc <- log(data$DOC)
data$lndon <- log(data$DON +1)
data$lnnh4 <- log(data$NH4i +1)
data$lnno3 <- log(data$NO3i +1)
data$lnmbc <- log(data$MBC +1)
data$lnmbn <- log(data$MBN +1)
data$lndna <- log(data$dna)
data$lnag <- log(data$AG +1)
data$lnbg <- log(data$BG +1)
data$lncb <- log(data$CB +1)
data$lnxyl <- log(data$XYL +1)
data$lnnag <- log(data$NAG +1)
data$lnlap <- log(data$LAP)
data$lnphos <- log(data$PHOS)

data$zom <- (data$lnom - mean(data$lnom, na.rm=TRUE)) / sd(data$lnom, na.rm=TRUE)
data$zdoc <- (data$lndoc - mean(data$lndoc, na.rm=TRUE)) / sd(data$lndoc, na.rm=TRUE)
data$zdon <- (data$lndon - mean(data$lndon, na.rm=TRUE)) / sd(data$lndon, na.rm=TRUE)
data$znh4 <- (data$lnnh4 - mean(data$lnnh4, na.rm=TRUE)) / sd(data$lnnh4, na.rm=TRUE)
data$zno3 <- (data$lnno3 - mean(data$lnno3, na.rm=TRUE)) / sd(data$lnno3, na.rm=TRUE)
data$zmbc <- (data$lnmbc - mean(data$lnmbc, na.rm=TRUE)) / sd(data$lnmbc, na.rm=TRUE)
data$zmbn <- (data$lnmbn - mean(data$lnmbn, na.rm=TRUE)) / sd(data$lnmbn, na.rm=TRUE)
data$zdna <- (data$lndna - mean(data$lndna, na.rm=TRUE)) / sd(data$lndna, na.rm=TRUE)
data$zag <- (data$lnag - mean(data$lnag, na.rm=TRUE)) / sd(data$lnag, na.rm=TRUE)
data$zbg <- (data$lnbg - mean(data$lnbg, na.rm=TRUE)) / sd(data$lnbg, na.rm=TRUE)
data$zcb <- (data$lncb - mean(data$lncb, na.rm=TRUE)) / sd(data$lncb, na.rm=TRUE)
data$zxyl <- (data$lnxyl - mean(data$lnxyl, na.rm=TRUE)) / sd(data$lnxyl, na.rm=TRUE)
data$znag <- (data$lnnag - mean(data$lnnag, na.rm=TRUE)) / sd(data$lnnag, na.rm=TRUE)
data$zlap <- (data$lnlap - mean(data$lnlap, na.rm=TRUE)) / sd(data$lnlap, na.rm=TRUE)
data$zphos <- (data$lnphos - mean(data$lnphos, na.rm=TRUE)) / sd(data$lnphos, na.rm=TRUE)

data$cfi <- (data$zom + data$zdoc + data$zdon + data$znh4 + data$zno3 + data$zmbc + data$zmbn + data$zdna +
               data$zag + data$zbg + data$zcb + data$zxyl + data$znag + data$zlap + data$zphos) / 15


#select the data needed for random forest
#rfdata <- data[,c(473,26:38,43:472)] #single pixels
rfdata <- data[,c(5,24:336)]#3mbuffs, all exp vars
rfdata <- rfdata[complete.cases(rfdata),]
respvar <- "pH"
#rfdata$zRV <- scale(rfdata$lnRV, center=TRUE, scale=TRUE)
#rfdata <- rfdata[,-1]

###############DO THE LOOP###################
set.seed(100)
#NUMBER OF ITERATIONS
n.loop<-100
#9:48
##running loop on 92 obs/443 varibales, 100 times takes about 14 minutes
#create matrices for storage
rsq.frame <- matrix(NA, nrow=n.loop, ncol=1)
mse.frame <- matrix(NA, nrow=n.loop, ncol=1)

pb <- txtProgressBar(min = 2, max = 100, style = 3)

#start loop
for(i in 1:n.loop){
  model1 <- randomForest(pH ~ ., data = rfdata, ntree=2000, importance = TRUE, parallel=TRUE)
  #add R2 to otuput matrix
  rsq.frame[i,] <- model1$rsq[2000]*100
  #add error val to output matrix
  mse.frame[i,] <- model1$mse[2000]
  #show progress bar
  setTxtProgressBar(pb, i)
}

rsq.frame <- data.frame(rsq.frame)
mse.frame <- data.frame(mse.frame)

write.csv(rsq.frame, paste("./3mbuff/RFresults/7April/r2_", respvar, ".csv", sep=""))
write.csv(mse.frame, paste("./3mbuff/RFresults/7April/mse_", respvar, ".csv", sep=""))

model1
mean(rsq.frame$rsq.frame)
mean(mse.frame$mse.frame)


#write.csv(rsq.frame, "./3mbuff/RFresults/try2/vsurf_red/r2_om.csv")
#write.csv(mse.frame, "./3mbuff/RFresults/try2/vsurf_red/mse_om.csv")


######Variable reduction/importance with the VSURF package#########
#rfdata <- data[,c(4,24:36,41:470)] #3mbuffs
#rfdata <- rfdata[complete.cases(rfdata),]

#model1 <- randomForest(XYL ~ ., data = rfdata, ntree=2000, importance = TRUE, parallel=TRUE)
#model1

varsel <- VSURF(pH ~ ., data=rfdata, parallel=TRUE, ncores=32)
summary(varsel)
number <- c(2:314)
thresvars <- number[varsel$varselect.thres]
interpvars <- number[varsel$varselect.interp]
predvars <- number[varsel$varselect.pred]

col.names <- names(rfdata)
col.names[c(thresvars)]
col.names[c(interpvars)]
col.names[c(predvars)]

plot(varsel)
plot(varsel, step="pred", var.names=TRUE)

rfdata2 <- rfdata[,c(1,predvars)]
model2 <- randomForest(pH ~ ., data = rfdata2, ntree=2000, importance = TRUE, parallel=TRUE)
model2



##############################################################################################
#######extracting pH from the raster prediction and adding it to RF predictor vars###########
setwd("Z:/SRER/Martha/RFdata")
#load raster and soil sample locs
r <- raster("./RFresults/correct_extent/ph.tif")
pts <- readOGR("soildata.shp")

#extract raster values for soil sample locs; sp=TRUE adds the extracted value to the other values of that spatial points df
rasValue <- raster::extract(r, pts, sp=TRUE)

#convert into a data frame
soildf <- as.data.frame(rasValue)

ph <- soildf[,c(1,26)]
names(ph)[names(ph) == "samplID"] <- "sampleID"

#load RFdata
rfdata <- read.csv("./3mbuff/rfdata.csv")

data <- merge(rfdata, ph, by="sampleID")

#delete duplicate columns
data <- data[,-c(2,6,24)]
#replace negative MBC values with NA
data$MBN[data$MBN < 0] <- NA
data$MBC[data$MBC < 0] <- NA

str(data)

#brightness normalization per Feilhauer et al 2010
ref <- data[,39:464]
ref <- as.matrix(ref)
#create the brightness normalized function
bnorm <- function (X)
{ X <- X / sqrt (rowSums (X ^ 2, na.rm=T))}
#apply the brightness normalized function
bnref <- bnorm(ref)
bnref <- as.data.frame(bnref)
#exclude the water absorption regions at 937 and 1127
x1 <- dplyr::select(bnref, `b_436.362396`:`b_922.151001`)
x2 <- dplyr::select(bnref, `b_962.216003`:`b_1107.451782`)
x3 <- dplyr::select(bnref, `b_1162.54126`:`b_1317.793213`)
x4 <- dplyr::select(bnref, `b_1488.06958`:`b_1768.524902`)
x5 <- dplyr::select(bnref, `b_2038.963867`:`b_2434.606201`)
x <- cbind(x1, x2, x3, x4, x5)

#combine brightness normalized values with descriptive variables
descvars <- data[,c(1,2,465,3:38)]
data <- cbind(descvars, x)

str(data)

#change character variables into factors
data$landform <- as.factor(data$landform)
data$PM <- as.factor(data$PM)
data$soilage <- as.factor(data$soilage)

#summary(data$OM.)
#data$lnRV <- log(data$OM.)

rfdata <- data[,c(4,24:39,3,40:336)]#3mbuffs, all exp vars
respvar = "SWC"
rfdata <- rfdata[complete.cases(rfdata),]

#rfdata$zRV <- scale(rfdata$lnRV, center=TRUE, scale=TRUE)
#rfdata <- rfdata[,-1]

###############DO THE LOOP###################
set.seed(100)
#NUMBER OF ITERATIONS
n.loop<-100
##running loop on 92 obs/443 varibales, 100 times takes about 14 minutes
#create matrices for storage
rsq.frame <- matrix(NA, nrow=n.loop, ncol=1)
mse.frame <- matrix(NA, nrow=n.loop, ncol=1)

pb <- txtProgressBar(min = 2, max = 100, style = 3)

#start loop
for(i in 1:n.loop){
  model1 <- randomForest(SWC ~ ., data = rfdata, ntree=2000, importance = TRUE, parallel=TRUE)
  #add R2 to otuput matrix
  rsq.frame[i,] <- model1$rsq[2000]*100
  #add error val to output matrix
  mse.frame[i,] <- model1$mse[2000]
  #show progress bar
  setTxtProgressBar(pb, i)
}

rsq.frame <- data.frame(rsq.frame)
mse.frame <- data.frame(mse.frame)

#vsurf_red
#allexp
write.csv(rsq.frame, paste("./3mbuff/RFresults/7April/with_pH/r2_", respvar, ".csv", sep=""))
write.csv(mse.frame, paste("./3mbuff/RFresults/7April/with_pH/mse_", respvar, ".csv", sep=""))

#model1
mean(rsq.frame$rsq.frame)
mean(mse.frame$mse.frame)

##VSURF VARIABLE REDUCTION##
varsel <- VSURF(SWC ~ ., data=rfdata, parallel=TRUE, ncores=32)
summary(varsel)
number <- c(2:315)
thresvars <- number[varsel$varselect.thres]
interpvars <- number[varsel$varselect.interp]
predvars <- number[varsel$varselect.pred]

col.names <- names(rfdata)
col.names[c(thresvars)]
col.names[c(interpvars)]
col.names[c(predvars)]

plot(varsel)
plot(varsel, step="pred", var.names=TRUE)

rfdata2 <- rfdata[,c(1,predvars)]

###RE-run the RF loop with reduced exp vars
rsq.frame <- matrix(NA, nrow=n.loop, ncol=1)
mse.frame <- matrix(NA, nrow=n.loop, ncol=1)

pb <- txtProgressBar(min = 2, max = 100, style = 3)

#start loop
for(i in 1:n.loop){
  model1 <- randomForest(SWC ~ ., data = rfdata2, ntree=2000, importance = TRUE, parallel=TRUE)
  #add R2 to otuput matrix
  rsq.frame[i,] <- model1$rsq[2000]*100
  #add error val to output matrix
  mse.frame[i,] <- model1$mse[2000]
  #show progress bar
  setTxtProgressBar(pb, i)
}

rsq.frame <- data.frame(rsq.frame)
mse.frame <- data.frame(mse.frame)

mean(rsq.frame$rsq.frame)
mean(mse.frame$mse.frame)



















###########Variable importance with RFpermute package#######
library(rfPermute)
set.seed(100)
model1.rfP <- rfPermute(BG ~ ., data = rfdata, ntree=2000, importance=TRUE)
# Plot the unscaled importance distributions and highlight significant predictors
plot(rp.importance(model1.rfP, scale = TRUE))
rp.importance(model1.rfP, scale = TRUE)


# To check important variables
importance(model1)        
varImpPlot(model1) 

