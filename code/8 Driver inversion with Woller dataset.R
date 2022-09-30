library(coda)
library(rjags)
library(R2jags)
library(mcmcplots)
library(bayestestR)
library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)

plot.col<-viridis(7)

setwd("C:/Users/ydmag/Google Drive/U of U/Elephant movement/Sr-in-ivory")

####invterting laser abalition data from Wooller et al., 2021####
Wooller <- read.csv("data/Wooller_Data_S3.csv")

#dist is in cm, convert to mm
wooller.micron <- Wooller$Dist_Seg01*10000
n.avg <- 100 #now many data points to average
n <- floor(length(Wooller$Sr_Seg01)/n.avg) #4287 data points, discard the last <100 data points
Wooller.sr <- Wooller$Sr_Seg01
#initiate vectors
n.avg.Wooller.sr <- rep(0, n)
n.sd.Wooller.sr <- rep(0, n) 

n.avg.Wooller.dist <- rep(0, n) 
n.med.Wooller.dist <- rep(0, n) 
  
for(i in 1:n){
  x <- ((i-1)*n.avg + 1):(i*n.avg)
  n.avg.Wooller.sr[i] <- mean(Wooller.sr[x])
  n.sd.Wooller.sr[i] <- sd(Wooller.sr[x])
  
  n.avg.Wooller.dist[i] <- mean(wooller.micron[x])#mean and median are the same
  n.med.Wooller.dist[i] <- median(wooller.micron[x])
}

plot(n.avg.Wooller.dist, n.avg.Wooller.sr,type="l")
lines(n.med.Wooller.dist, n.avg.Wooller.sr, col="red")

####subseting the entire dataset
sub <- 1000:1300
sub.n.avg.dist <- rev(n.avg.Wooller.dist[sub])
sub.n.avg.sr <- rev(n.avg.Wooller.sr[sub])
sub.n.sd.sr <- rev(n.sd.Wooller.sr[sub])

plot(sub.n.avg.dist, sub.n.avg.sr, type="l", xlim=c(max(sub.n.avg.dist),min(sub.n.avg.dist)))

##ivory rate

##body mass scaling

#extract mean values every 100 data points or every 7 days worth of ivory
