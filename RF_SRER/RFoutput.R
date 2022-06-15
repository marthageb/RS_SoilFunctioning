library(tidyverse)
library(caret)
library(raster)
library(rgdal)


mast <- raster("./same_extent/NEON_AGB.tif")
soil <- raster("./RFresults/ph.tif")

#files <- list()
#files[[1]] <- mast
#files[[2]] <- soil
###test rasters
test <- list()

#for(i in 1:length(files)){
#  test[[i]] <- capture_warnings(compareRaster(mast,files[[i]], res=T, orig=T, stopiffalse=F, showwarning=T))
#}
test

soilr <- raster(vals=values(soil), ext=extent(mast), crs=crs(mast), nrows=dim(mast)[1], ncols=dim(mast)[2])

e <- extent(mast)
rc <- crop(soilr, e)
if(sum(as.matrix(extent(rc))!=as.matrix(e)) == 0){ # edited
  rc <- mask(rc, mast) # You can't mask with extent, only with a Raster layer, RStack or RBrick
}else{
  rc <- extend(rc,mast)
  rc<- mask(rc, mast)
}


files <- list()
files[[1]] <- mast
files[[2]] <- rc
rfdata<- stack(files)

writeRaster(rc, "./RFresults/correct_extent/ph.tif", overwrite=TRUE)


########extract RF predicted values at soil sample locs###########
setwd("Z:/SRER/Martha/RFdata")
#load raster and soil sample locs
r <- raster("./RFresults/correct_extent/ph.tif")
pts <- readOGR("soildata.shp")

#extract raster values for soil sample locs; sp=TRUE adds the extracted value to the other values of that spatial points df
rasValue <- raster::extract(r, pts, sp=TRUE)

#convert into a data frame
soildf <- as.data.frame(rasValue)

###perform the linear regressions and look at results###
xvar <- soildf$pH
yvar <- soildf$ph

#xvar <- xvar * 1000 #for SWC
#yvar <- yvar * 1000 #for SWC

#calculate R2
cor1 <- (cor(xvar, yvar, use="pairwise.complete.obs"))^2
cor1
# Build the model
regmodel <- lm(yvar ~ xvar, na.action = na.omit)
summary(regmodel)
#residual sum of squares
RSS <- c(crossprod(regmodel$residuals))
#mean squared error
MSE <- RSS / length(regmodel$residuals)
MSE
#root mean squared error
RMSE <- sqrt(MSE)
RMSE

###THEN GO TO "REGPLOTPARMS.R TO CREATE PLOTS###

#set the font
windowsFonts(A=windowsFont("Calibri Light"))

#practice plot to see axes ranges
plot(xvar,yvar)

plot(x=1, xlab="", ylab = "", xlim = c(0,2.5), ylim = c(0,2.5), axes=FALSE, type="n")
#add y-axis title
title(ylab="Predicted OM %", family="A", line=2.2)
#add x-axis title
title(xlab="Actual OM %", family="A", line=2.2)
#add y-axis labels
axis(side=2,at=ylabloc,labels=ylabs, family="A")
#add x-axis labels
axis(side=1,at=xlabloc,labels=xlabs, family="A")
#add the box around the plot area
box(which="plot")

#add points and lines
points(xvar, yvar, pch= 16, cex = 0.5)
abline(regmodel, lwd=1.5)
# abline(0,1, lwd=1.5, lty=2) #if wanting to include 1:1 line

#add reg. stats
Rsq <- expression(paste("R"^"2"*" = 0.921"))
#pval <- expression(paste("p < 2.22 x 10"^"-16"))
mselab <- expression(paste("mse = 0.0111"))
text(2.32,0.3, labels=Rsq, family="A")
#text(2,0.18, labels=pval, family="A")
text(2.25,0.1, labels=mselab, family="A")

mtext("a.", side=3, at=0, cex=1.25, family="A")



