library(FactoMineR)

setwd("Z:/SRER/Martha/hyperspec/flightlines/FoliarN")
data <- read.csv("NDVI0.2_thick_sparseALL.csv")

data[data == -1] <-NA


data <- plants
#select a single plot
plot <- subset(data, plotID == "G2mesq")
plot <- plot[complete.cases(plot[,c(11:434)]), ]
#select 10 rows from that plot at random
#plot <- plot[sample(nrow(plot), 10), ]

#select only the reflectance readings
tmp <- plot[,9:434]
#remove reflectance at 381 since it is 'NA' for alot of samples
tmp <- tmp[,-1]
#run the PCA
pca1 <- PCA(tmp, scale=F, graph=F)

summary(pca1)
#construct confidence ellipse values
ID <- rep("plot1", times=nrow(tmp))
pca_vals <- cbind.data.frame(ID, pca1$ind$coord)
pca_vals$ID <- as.factor(pca_vals$ID)
confellipse <- coord.ellipse(pca_vals, bary=F)

###For some reason the 'coord.ellipse' function isn't plotting correctly so skip to line 27 to plot manually
#ID <- rep("A2mesq", times=nrow(tmp))
#pca_vals <- cbind.data.frame(ID, pca1$ind$coord)
#pca_vals$ID <- as.factor(pca_vals$ID)
#plot.PCA(pca1, ellipse=confellipse$res, xlim=c(-10,10), ylim=c(-10,10))

#make a dataframe of the plot data locations
plotdata <- as.data.frame(pca1$ind$coord)
plot(plotdata$Dim.1, plotdata$Dim.2, xlim=c(-5,5), ylim=c(-5,5))
#label the points on the plot
text(Dim.2 ~Dim.1, labels=rownames(plotdata),data=plotdata, cex=0.9, font=2)
#add the confidence ellipse to the plot
polygon(confellipse$res$Dim.1, confellipse$res$Dim.2, density=0, col='red')

###ONCE YOU ENSURE ALL OF THE PIXELS FALL INSIDE THE CONF ELLIPSE THEN RUN THE NEXT 2 LINES
#append dataframe with the 'good samples'
good <- data[rownames(plot),]
#finalPCA <- good #for the 1st one
finalPCA <- rbind.data.frame(good,finalPCA)


#write.csv(finalPCA, "Z:/SRER/Martha/hyperspec/flightlines/noNDVI_TS_withbare.csv")




summary(data$plotID)
















library(plyr)
data <- read.csv("allplants.csv")
xycoords <- read.csv("C:/Users/marthag/Desktop/envi_spectralJan2020/xycoords.csv")
xy <- xycoords[,c(2,7,8)]
colnames(xy) <- c("ENVI_ID", "x", "y")
new <- match_df(data,xy)

new$plot <- as.factor(new$plot)
summary(new$plot)


dups <- data[duplicated(data[,c("x","y")]),]
no_dups <- data[!duplicated(data[,c("x","y")]),]
data <- no_dups
