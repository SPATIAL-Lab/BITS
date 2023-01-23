library(coda)
library(rjags)
library(R2jags)
library(mcmcplots)
library(bayestestR)
library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)
library(zoo)

plot.col<-viridis(7)

setwd("C:/Users/ydmag/Google Drive/U of U/Elephant movement/Sr-in-ivory")

source("code/1 Helper functions.R")

####supplementary figures######

#Fig S1 
#no associated code 

#Fig S2 
#comparing solution methods to laser ablation methods
misha.micromill <- read.csv("data/Misha micromill 400-600.csv")

misha.micromill<-na.omit(misha.micromill)
misha.micromill.dist <- misha.micromill$Position * 547

misha.micromill.sr <- misha.micromill$corr..87Sr.86Sr
misha.micromill.sr.err <- misha.micromill$comb..Err

plot(n.avg.misha.50.dist, n.avg.misha.50.sr,type = "l",col="#00b4ffff",lwd=1.5,
     xlim=c(23000,0),ylim=c(0.705,0.712),xlab="Distance from pulp cavity (microns)",ylab="87Sr/86Sr",
     main="Comparing data from the LA-ICP-MS and the solution methods")
PlotPE(n.avg.misha.50.dist, n.avg.misha.50.sr,n.sd.misha.50.sr,col="#00b4ffff")
points(n.avg.misha.50.dist, n.avg.misha.50.sr, col="#00b4ffff",pch=18,cex=1.1)
PlotPE(misha.micromill.dist,misha.micromill.sr,misha.micromill.sr.err,col="red")
legend(10000, 0.707, c("50-point avg. LA-ICP-MS","Micromill"),
       pch = c(16, 16), col=c("#00b4ffff","red"))

#######Fig S3 sensitivity test Sensitivity to excursion
#excursion makes rate parameter smaller

#Fig S3
#sensitivity tests section 3#
#cauchy vs cauchy with autocorrelation, and show distribution of the autocorrelation term  



#Fig S4: #cauchy vs normal with autocorrelation


#sensitivity to Re (this might be a thing!)

#Fig S5: forward model

