library(pls)
library(plsVarSel)
library(dplyr)
library(ape) #for Moran's I calculation


setwd("Z:/SRER/Martha/hyperspec/flightlines")
Nref <- read.csv("noNDVI_TS_withbare.csv")
Nref <- Nref[,-1]

#delete B2 and D1 grass
#bad <- c("B2gras", "D1gras")
#Nref <- Nref[!(Nref$plotID %in% bad),]


#change col names (all say 'b_' right now)
#names(Nref) <- gsub("b_", "", names(Nref))

#brightness normalization per Feilhauer et al 2010
ref <- Nref[,9:434]
ref <- as.matrix(ref)
  #create the brightness normalized function
  bnorm <- function (X)
  { X <- X / sqrt (rowSums (X ^ 2, na.rm=T))}
#apply the brightness normalized function
bnref <- bnorm(ref)

#combine brightness normalized values with descriptive variables
descvars <- Nref[,1:8]
NrefBN <- cbind(descvars, bnref)

#trim data to remove begining/end portions and water absorption bands from the brighness normalized reflectance vals
x1 <- dplyr::select(NrefBN, `X436.362396`:`X1317.793213`)
x2 <- dplyr::select(NrefBN, `X1478.053345`:`X1768.524902`)
x3 <- dplyr::select(NrefBN, `X2038.963867`:`X2434.606201`)
x <- cbind(x1, x2, x3)

#exclude the water absorption regions at 937 and 1127
x1 <- dplyr::select(NrefBN, `X436.362396`:`X922.151001`)
x2 <- dplyr::select(NrefBN, `X962.216003`:`X1107.451782`)
x3 <- dplyr::select(NrefBN, `X1162.54126`:`X1317.793213`)
x4 <- dplyr::select(NrefBN, `X1488.06958`:`X1768.524902`)
x5 <- dplyr::select(NrefBN, `X2038.963867`:`X2434.606201`)
x <- cbind(x1, x2, x3, x4, x5)

############################PLSR Leave-one-PIXEL-out#######################
#combine plot lvl leaf N content and selected wavelengths
##convert selected wavelengths into SPC (format needed for PLSR)
foliar <- data.frame(NrefBN[,c(1,2)], SPC=I(as.matrix(x)))

rownames(foliar) <- seq(1:nrow(foliar))
#SET SEED TO INITIALLY CONTROL "RANDOMNESS"
set.seed(18)
#NUMBER OF ITERATIONS
n.loop<-100

#create matrices for storage

#unsure if length of TestIDs will always be the same
min.err.frame<-matrix(NA, nrow=n.loop, ncol=1)
min.comp.frame<-matrix(NA, nrow=n.loop, ncol=1)
var.frame<-matrix(NA, nrow=n.loop, ncol=1)
pred.err.frame<-matrix(NA, nrow=n.loop, ncol=1)
pred.var.frame<-matrix(NA, nrow=n.loop, ncol=1)
regcoef.frame <- matrix(NA, nrow=length(x)+1, ncol=n.loop)
rn <- c("yint", colnames(x))
row.names(regcoef.frame) <- rn
residuals.df <- matrix(NA, nrow=294, ncol=n.loop)
VIP.frame <- matrix(NA, nrow=length(x), ncol=n.loop)

pb <- txtProgressBar(min = 1, max = n.loop, style = 3)
#start loop
for(i in 1:n.loop){
  #create a training and testing dataset
  testIDs <- sample(1:420, 126, replace=FALSE)
  datTrain <- subset(foliar, !rownames(foliar) %in% testIDs)
  datTest <- subset(foliar, rownames(foliar) %in% testIDs)
  #str(datTrain)
  #str(datTest)
  
  TestIDs <- rownames(datTest)
  
  
  #define the model parameters and run the plsr model
  #(response_var ~ predictor_var, data)
  #10-fold cross validation (validation="CV", segments= 10)
  plsFit <- plsr(plotN ~ SPC, data = datTrain, validation="CV", segments= 10)
  # summary(plsFit)
  #calculate the minimum RMSEP
  error <-  RMSEP(plsFit, estimate="adjCV")
  min_err <- min(error$val)
  #min_err
  min.err.frame[i,]<-min_err
  
  #determine the # of comps that corresponds with the minimum RMSEP
  #need to subtract 1 bacause indexing starts at intercept (not the 1st comp)
  min_comp <- which.min(error$val)-1
  #min_comp
  min.comp.frame[i,]<-min_comp
  
  if (min_comp > 0){
    #calculate the R2 for # of components defined above
    var <- R2(plsFit, estimate="train", comps=1:min_comp)
    var
    #r2 for model
    var.frame[i,] <- var$val
    #use the plsr model and the optimal # of comps IDed to predict plot N values that were left out of model development
    plsPred <- predict(plsFit, datTest, ncomp=min_comp)
    # evaluate how well the model performed
    pls.eval=data.frame(obs=datTest$plotN, pred=plsPred[,1,1])
    ##needs caret package... will give you RMSE and R2 for the test data
    #DOUBLE COLON!
    info <- caret::defaultSummary(pls.eval)
    pred_err <- as.numeric(info[1:1])
    pred_var <- as.numeric(info[2:2])
    pred.err.frame[i,] <- pred_err
    pred.var.frame[i,] <- pred_var
    
    #get reg coefficients and residuals from plsr model based on optimal # of comps IDed
    regcoef.frame[,i]<- coef(plsFit, ncomp=min_comp, intercept=TRUE)
    residuals.df[,i] <- plsFit$residuals[,,min_comp]
    
    #get the VIPscore for the minimum # of components defined
    VIP.frame[,i] <- VIP(plsFit, min_comp)
  } else {
    var.frame[i,] <- NA
    pred.err.frame[i,] <- NA
    pred.var.frame[i,] <- NA
    regcoef.frame[,i] <- NA
    residuals.df[,i] <- NA
    VIP.frame[,i] <- NA
  }
  setTxtProgressBar(pb, i)
}
  
#can convert to data.frame
min.err.frame <- data.frame(min.err.frame)
min.comp.frame <- data.frame(min.comp.frame)
var.frame <- data.frame(var.frame)
pred.err.frame <- data.frame(pred.err.frame)
pred.var.frame <- data.frame(pred.var.frame)
regcoef.frame <- data.frame(regcoef.frame)
regcoef.frame <- t(regcoef.frame)
VIP.frame <- data.frame(VIP.frame)

results <- cbind(min.err.frame, min.comp.frame, var.frame, pred.err.frame, pred.var.frame)
colnames(results) <- c("trainRMSE", "ncomp", "trainr2", "testRMSE", "testr2")

mean(results$trainr2)
mean(results$trainRMSE)
mean(results$testr2)
mean(results$testRMSE)


#write.csv(results, "./PLSR/PixelLOO/results.csv")
#write.csv(regcoef.frame, "./PLSR/PixelLOO/regcoeff.csv")
#write.csv(VIP.frame, "./PLSR/PixelLOO/VIP.csv")


# create df of results where [R2 validation > avg R2] and [RMSE validation < avg RMSE] (per Chadwick & Asner 2016 RSE)
goodmodel <- subset(results, testRMSE < median(testRMSE))
goodmodel <- subset(goodmodel, testr2 > median(results$testr2))

nrow(goodmodel)
mean(goodmodel$testr2)
mean(goodmodel$testRMSE)



############################PLSR Leave-one-PLOT-out#######################
#create a list of the plots
plots <- as.character(NrefBN$plotID)
plots <- unique(plots)


#combine plot lvl leaf N content and selected wavelengths
##convert selected wavelengths into SPC (format needed for PLSR)
foliar <- data.frame(NrefBN[,c(1,2)], SPC=I(as.matrix(x)))
foliar$IDs <- rownames(foliar)

#SET SEED TO INITIALLY CONTROL "RANDOMNESS"
set.seed(18)
#NUMBER OF ITERATIONS
n.loop<-100

#create matrices for storage

#unsure if length of TestIDs will always be the same
min.err.frame<-matrix(NA, nrow=n.loop, ncol=1)
min.comp.frame<-matrix(NA, nrow=n.loop, ncol=1)
var.frame<-matrix(NA, nrow=n.loop, ncol=1)
pred.err.frame<-matrix(NA, nrow=n.loop, ncol=1)
pred.var.frame<-matrix(NA, nrow=n.loop, ncol=1)
regcoef.frame <- matrix(NA, nrow=length(x)+1, ncol=n.loop)
row.names(regcoef.frame) <- c("intercept",colnames(x))
VIP.frame <- matrix(NA, nrow=length(x), ncol=n.loop)
row.names(VIP.frame) <- colnames(x)
residuals.df <- matrix(NA, nrow=290, ncol=n.loop)
trainID.df <- matrix(NA, nrow=290, ncol=n.loop)

pb <- txtProgressBar(min = 1, max = n.loop, style = 3)
#start loop
for(i in 1:n.loop){
  
  #create a training and testing dataset
    #Shuffle plots and select 30% to be the testing IDs: 
    testIDs <- sample(plots,13)
    #select the remaining IDs to be the training plots:
    trainIDs <- subset(plots, !(plots%in% testIDs))
    
    #Divide data set into training and test set (70% training and 30% testing)
    datTrain <- subset(foliar, plotID %in% trainIDs)
    datTest <- subset(foliar, plotID %in% testIDs)
    trainID.df[,i] <- datTrain$IDs
  #str(datTrain)
  #str(datTest)
  
  #define the model parameters and run the plsr model
  #(response_var ~ predictor_var, data)
  #10-fold cross validation (validation="CV", segments= 10)
  plsFit <- plsr(plotN ~ SPC, data = datTrain,validation="CV", segments= 10)
  # summary(plsFit)
  #calculate the minimum RMSEP
  error <-  RMSEP(plsFit, estimate="adjCV")
  min_err <- min(error$val)
  #min_err
  min.err.frame[i,]<-min_err
  
  #determine the # of comps that corresponds with the minimum RMSEP
  #need to subtract 1 bacause indexing starts at intercept (not the 1st comp)
  min_comp <- which.min(error$val)-1
  #min_comp
  min.comp.frame[i,]<-min_comp
  
  if (min_comp > 0){
    #calculate the R2 for # of components defined above
    var <- R2(plsFit, estimate="train", comps=1:min_comp)
    var
    #r2 for model
    var.frame[i,] <- var$val
    #use the plsr model and the optimal # of comps IDed to predict plot N values that were left out of model development
    plsPred <- predict(plsFit, datTest, ncomp=min_comp)
    # evaluate how well the model performed
    pls.eval=data.frame(obs=datTest$plotN, pred=plsPred[,1,1])
    ##needs caret package... will give you RMSE and R2 for the test data
    #DOUBLE COLON!
    info <- caret::defaultSummary(pls.eval)
    pred_err <- as.numeric(info[1:1])
    pred_var <- as.numeric(info[2:2])
    pred.err.frame[i,] <- pred_err
    pred.var.frame[i,] <- pred_var
    
    #get reg coefficients and residuals from plsr model based on optimal # of comps IDed
    regcoef.frame[,i]<- coef(plsFit, ncomp=min_comp, intercept=TRUE)
    residuals.df[,i] <- plsFit$residuals[,,min_comp]
      resd <- as.data.frame(plsFit$residuals[,,min_comp])
    #get the VIPscore for the minimum # of components defined
    VIP.frame[,i] <- VIP(plsFit, min_comp)
  } else {
    var.frame[i,] <- NA
    pred.err.frame[i,] <- NA
    pred.var.frame[i,] <- NA
    regcoef.frame[,i] <- NA
    residuals.df[,i] <- NA
    VIP.frame[,i] <- NA
  }
  setTxtProgressBar(pb, i)
  }


#can convert to data.frame
min.err.frame <- data.frame(min.err.frame)
min.comp.frame <- data.frame(min.comp.frame)
var.frame <- data.frame(var.frame)
pred.err.frame <- data.frame(pred.err.frame)
pred.var.frame <- data.frame(pred.var.frame)
regcoef.frame <- data.frame(regcoef.frame)
VIP.frame <- data.frame(VIP.frame)

regcoef.frame <- t(regcoef.frame)


results <- cbind(min.err.frame, min.comp.frame, var.frame, pred.err.frame, pred.var.frame, regcoef.frame)
colnames(results) <- c("trainRMSE", "ncomp", "trainr2", "testRMSE", "testr2", colnames(regcoef.frame))

results <- results[complete.cases(results), ]

mean(results$trainr2)
mean(results$trainRMSE)
mean(results$testr2)
mean(results$testRMSE)

#write.csv(results, "Z:/SRER/Martha/hyperspec/flightlines/FoliarN/PLSR/PlotLOO/results.csv")
#write.csv(regcoef.frame, "Z:/SRER/Martha/hyperspec/flightlines/FoliarN/PLSR/PlotLOO/regcoeff.csv")
#write.csv(VIP.frame, "Z:/SRER/Martha/hyperspec/flightlines/FoliarN/PLSR/PlotLOO/VIP.csv")


# create df of results where [R2 validation > avg R2] and [RMSE validation < avg RMSE] (per Chadwick & Asner 2016 RSE)
goodmodel <- subset(results, testRMSE < median(testRMSE))
goodmodel <- subset(goodmodel, testr2 > median(results$testr2))
nrow(goodmodel)

mean(goodmodel$testr2)
mean(goodmodel$testRMSE)


###########TESTING RESIDUALS###################
#convert the residuals and training IDS results to dataframes
residuals.df <- data.frame(residuals.df)
trainID.df <- data.frame(trainID.df)

#only keep the results that were included in the "good model"
goodrun <- rownames(goodmodel)
goodresiduals <- residuals.df %>% dplyr::select(all_of(goodrun))
goodtrainIDs <- trainID.df %>% dplyr::select(all_of(goodrun))
res_test <- cbind.data.frame(goodresiduals[,1],goodtrainIDs[,1])
colnames(res_test) <- c("residual", "IDs")
#res_test$ENVI_ID <- as.character(res_test$ENVI_ID)                          


xycoords <- Nref[,7:8]
xycoords$IDs <- as.numeric(rownames(xycoords))
goodcords <- subset(xycoords, rownames(xycoords) %in% res_test[,2])

sds <- merge(goodcords, res_test, by = "IDs")

#create a distance matrix
sds.dist <- as.matrix(dist(cbind(sds$x, sds$y)))
#take the inverse of the matrix values
sds.dist.inv <- 1/sds.dist
sds.dist.inv[is.infinite(sds.dist.inv)] <- 0
#replace the diagonal entries with zero:
diag(sds.dist.inv) <- 0
sds.dist.inv [1:5,1:5]
#calculate Moran's I
Moran.I(sds$residual, sds.dist.inv, scaled=FALSE) 


####do the above in a loop for all good model results
moranI.frame <- matrix(NA, nrow=length(goodresiduals), ncol=1)

pb <- txtProgressBar(min = 1, max = n.loop, style = 3)

for (i in 1:length(moranI.frame)){
  #dataframe with the residuals and IDs for each of the PLSR runs that were deemed "good"
  res_test <- cbind.data.frame(goodresiduals[,i],goodtrainIDs[,i])
  colnames(res_test) <- c("residual", "IDs")

  #get the coordinates for these pixels and merge coordinates with residual values
  goodcords <- subset(xycoords, rownames(xycoords) %in% res_test[,2])
  sds <- merge(goodcords, res_test, by = "IDs")
  #create a distance matrix
  sds.dist <- as.matrix(dist(cbind(sds$x, sds$y)))
  #take the inverse of the matrix values
  sds.dist.inv <- 1/sds.dist
  sds.dist.inv[is.infinite(sds.dist.inv)] <- 0
  #replace the diagonal entries with zero:
  diag(sds.dist.inv) <- 0
  
  #calculate Moran's I based on the inverse distance matrix and residual values
  morI <- Moran.I(sds$residual, sds.dist.inv)
  moranI.frame[i,] <- morI$p.value
  setTxtProgressBar(pb, i)
}

moranI.frame <- data.frame(moranI.frame)
moranI.frame$PLSRID <- colnames(goodresiduals)
nonsig <- filter_at(moranI.frame, vars(moranI.frame), any_vars(.>0.05))

###########subset results to include only models that are robust AND no spatial auto correlation

finalvals <- subset(goodmodel, rownames(goodmodel) %in% nonsig$PLSRID)
mean(finalvals$testr2)                     
mean(finalvals$testRMSE) 
#mean(finalvals$trainr2)
#mean(finalvals$trainRMSE)

finalcoef <- finalvals[,6:303]
coeffs <- as.matrix(colMeans(finalcoef))    
sdcoef <- as.matrix(apply(finalcoef, 2, sd))


vip.good <- t(VIP.frame)
vip.good <- subset(vip.good, rownames(vip.good) %in% rownames(finalcoef))
vips <- as.matrix(colMeans(vip.good))    
sdvip <- as.matrix(apply(vip.good, 2, sd))


results2 <- cbind.data.frame(coeffs[-1,], sdcoef[-1,], vips, sdvip)
colnames(results2) <- c("coef", "sd_coef", "vip", "sd_vip")


write.csv(finalcoef, "Z:/SRER/Martha/hyperspec/flightlines/FoliarN/PLSR/PlotLOO/PlotLOO_goodcoeff_noNDVI_withbare.csv")
write.csv(coeffs, "Z:/SRER/Martha/hyperspec/flightlines/FoliarN/PLSR/PlotLOO/PlotLOO_avggoodcoeff_noNDVI_withbare.csv")
write.csv(results2, "Z:/SRER/Martha/hyperspec/flightlines/FoliarN/PLSR/PlotLOO/results2_noNDVI_withbare.csv")






###get PLSR coeffs in correct format for 'analyze_reflect.py'
coef <- read.csv("Z:/SRER/Martha/hyperspec/flightlines/FoliarN/PLSR/PlotLOO/PlotLOO_avggoodcoeff_noNDVI_withbare.csv")

#get the y intercept value (copy and paste onto 'analyze_reflect')
int <- coef[1,2]

#remove the first row
coef <- coef[-1,]

#rename columns
colnames(coef)<- c("lambda", "PLSRcoef")

#remove 'X' from wavelengths
coef$lambda <- sub("^X", "", coef$lambda)
rownames(coef) <- coef$lambda

#add the wavelength regions that were excluded
#create vector of excluded bands
wavelengths <- sub("^X", "", colnames(ref))
lambda <- wavelengths[!wavelengths %in% coef$lambda]

#add "NA" columns
colNAs <-matrix(NaN, nrow=length(lambda), ncol=2)
colNAs <- as.data.frame(colNAs)
rownames(colNAs) <- lambda
colnames(colNAs) <- c("lambda", "PLSRcoef")
colNAs$lambda <- rownames(colNAs)

allvals <- rbind.data.frame(coef,colNAs)

#reorder data in wavelength order
lamnum <- rownames(allvals)
lamnum <- as.numeric(lamnum)
#assign the numeric wavelength as column and row names
rownames(allvals) <- lamnum
#order the df in decending wavelength order
order <- allvals[ order(as.numeric(row.names(allvals))), ]

#reassign NaN values 0
order$PLSRcoef[is.nan(order$PLSRcoef)] <- 0

write.csv(order, "Z:/SRER/Martha/hyperspec/flightlines/PLSRcoef.csv")

format(round(int, 14), nsmall = 14)
