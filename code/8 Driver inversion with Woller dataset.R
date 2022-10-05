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
}

plot(n.avg.Wooller.dist, n.avg.Wooller.sr,type="l")
lines(n.med.Wooller.dist, n.avg.Wooller.sr, col="red")

####subseting the entire dataset
sub <- 1001:1300
sub.n.avg.dist <- rev(n.avg.Wooller.dist[sub])
sub.n.avg.sr <- rev(n.avg.Wooller.sr[sub])
sub.n.sd.sr <- rev(n.sd.Wooller.sr[sub])

plot(sub.n.avg.dist, sub.n.avg.sr, type="l", xlim=c(max(sub.n.avg.dist),min(sub.n.avg.dist)))

#foward model simulating micromill results (400 micron band)
mm.bwidth <- 400 #microns
index.wooller.dist<- ceiling(wooller.micron/mm.bwidth) #this is about the same as averaging per 100 data points!

#number of micromill simulations, discarding the last bit of data in the sequence
mm.sim.n <- max(index.wooller.dist) - 1 

mm.sim.avg.Wooller.sr <- rep(0,mm.sim.n)#initiate vectors
mm.sim.sd.Wooller.sr <- rep(0,mm.sim.n)
mm.sim.avg.Wooller.dist <- rep(0,mm.sim.n)

for(i in 1:mm.sim.n){
  temp.sr <- subset(Wooller.sr,index.wooller.dist==i)
  temp.dist <- subset(wooller.micron,index.wooller.dist==i)

  mm.sim.avg.Wooller.sr[i] <- mean(temp.sr)
  mm.sim.sd.Wooller.sr[i] <- sd(temp.sr)
  
  mm.sim.avg.Wooller.dist[i] <- median(temp.dist)#median is less sensitive to potential data gaps
}#takes ~30s 


plot(mm.sim.avg.Wooller.dist, mm.sim.avg.Wooller.sr,type="l")

####subseting the entire dataset
sub <- 1001:1300
sub.mm.sim.avg.dist <- rev(mm.sim.avg.Wooller.dist[sub])
sub.mm.sim.avg.sr <- rev(mm.sim.avg.Wooller.sr[sub])
sub.mm.sim.sd.sr <- rev(mm.sim.sd.Wooller.sr[sub])

plot(sub.mm.sim.avg.dist, sub.mm.sim.avg.sr, type="l", 
     xlim=c(max(sub.mm.sim.avg.dist),min(sub.mm.sim.avg.dist)))


######ivory extension rate Wooller et al 2021####
wooller.COSr<-read.csv("data/Wooller_isotope_data.csv")

wooller.rate <- rep(NA,27)

#first 2 years of life is neonate (Wooller et al 2021)
#last year is close to death, so these data are omitted
for(i in 3:27){
  test.sub.last <- subset(wooller.COSr$d,wooller.COSr$year==(i-1))
  test.sub <- subset(wooller.COSr$d,wooller.COSr$year==i)
  #d is in an increasing order, so get the first element of last year
  test.comb <- rbind(test.sub,test.sub.next[1])
  wooller.rate[i] <- 10000*(max(test.comb)- min(test.comb))/365 #distance in cm, 365 days in a year
  #results are in microns/day
}

wooller.rate<- na.omit(wooller.rate)
mean.wooller.rate <- mean(wooller.rate)
sd.wooller.rate <- sd(wooller.rate)
hist(wooller.rate[3:27])
plot(3:27,wooller.rate[3:27])

mean.wooller.rate*365*28 #this is close to the total length 1.7 meters

##body mass scaling



#extract mean values every 400 microns, if rate is 170, then 400 microns is around 3 days.
