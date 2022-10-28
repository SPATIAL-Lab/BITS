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

setwd("C:/Users/ydmag/Google Drive/U of U/Elephant movement/Sr-in-ivory")

##############################Try one pool model with 25 pt average and half the data#############
####also use 25 point average plus cut off at 14000#
#one pool scenario is supported#
misha.raw <- read.csv("data/Misha raw.csv")
misha.raw.Sr <- misha.raw$X87Sr.86Sr
misha.raw.dist <- misha.raw$dist

n.avg <- 25 #now many data points to average
n <- floor(length(misha.raw$X87Sr.86Sr)/n.avg) #99 data points, discard the last <100 data points

#initiate vectors
n.avg.misha.25.sr <- rep(0, n)
n.sd.misha.25.sr <- rep(0, n)

n.avg.misha.25.dist <- rep(0, n)

for(i in 1:n){
  x <- ((i-1)*n.avg + 1):(i*n.avg)
  n.avg.misha.25.sr[i] <- mean(misha.raw.Sr[x])
  n.sd.misha.25.sr[i] <- sd(misha.raw.Sr[x])
  
  n.avg.misha.25.dist[i] <- mean(misha.raw.dist[x])#mean and median are the same
}
plot(misha.raw$dist, misha.raw$X87Sr.86Sr, col=alpha("dark blue", 0.1))#plot all raw data points
plot(n.avg.misha.25.dist, n.avg.misha.25.sr, main="25 pt average",
     xlim=c(max(n.avg.misha.25.dist),min(n.avg.misha.25.dist)))
lines(n.avg.misha.25.dist, n.avg.misha.25.sr)

###########remove anomaly#######
ind.remv <- which(n.avg.misha.25.dist > 17000 & n.avg.misha.25.dist < 17500 & n.avg.misha.25.sr < 0.709)

index.25.anom.remv1 <- c(1:70)
index.25.anom.remv2 <- c(78:200)
index.25.anom.remv <- c(1:70,78:200)
n.avg.misha.25.dist.rmv <- n.avg.misha.25.dist[index.25.anom.remv]
n.avg.misha.25.sr.rmv <- n.avg.misha.25.sr[index.25.anom.remv]
n.sd.misha.25.sr.rmv <- n.sd.misha.25.sr[index.25.anom.remv]

plot(n.avg.misha.25.dist.rmv, n.avg.misha.25.sr.rmv, main="25 pt average",
     xlim=c(max(n.avg.misha.25.dist.rmv),min(n.avg.misha.25.dist.rmv)))
lines(n.avg.misha.25.dist[index.25.anom.remv1], n.avg.misha.25.sr[index.25.anom.remv1])
lines(n.avg.misha.25.dist[index.25.anom.remv2], n.avg.misha.25.sr[index.25.anom.remv2])
points(n.avg.misha.25.dist[ind.remv],n.avg.misha.25.sr[ind.remv],col= "red")
lines(n.avg.misha.25.dist[70:78],n.avg.misha.25.sr[70:78],col= "red")

###prep data###
R.sd.mea <- n.sd.misha.25.sr.rmv
dist.mea <- n.avg.misha.25.dist.rmv
R.mea <- n.avg.misha.25.sr.rmv
n.mea = length(n.avg.misha.25.sr.rmv)

R0 <- 0.707

#Re is the mean ratio of end value  
Re <- 0.711
s.intv <- 26.2

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.6
max.dist.mea <- max(n.avg.misha.25.dist) + 30 #adding a bit of distance leading into the series
#parameters to save
parameters <- c("a", "Ivo.rate", "Rs.m","Rin","dist","h.l",
                "Sr.pre", "Ps", "Fin", "Re.mean", "R0.mean", "switch")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 400, n.mea = n.mea, 
           R0 = R0, Re = Re, s.intv = s.intv, Ivo.rate.mean = Ivo.rate.mean,
           Ivo.rate.sd = Ivo.rate.sd, max.dist.mea = max.dist.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/400

#Run it
post.misha.fdnb1pr = do.call(jags.parallel,list(model.file = "code/Gabe JAGS 1 pool.R", 
                                                parameters.to.save = parameters, 
                                                data = dat, n.chains=5, n.iter = n.iter, 
                                                n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 3.5 hours

post.misha.fdnb1pr$BUGSoutput$summary
traplot(post.misha.fdnb1pr,parms = c("a"))

##plot calibration curves, run Helper functions first
plot(0,0, xlim = c(20000,13000), ylim = c(0.7066, 0.712), xlab = "distance", ylab ="Sr 87/86",
     main="One-pool model")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.fdnb1pr$BUGSoutput$sims.list$Rs.m,
               post.misha.fdnb1pr$BUGSoutput$sims.list$dist)
points(n.avg.misha.25.dist[index.25.anom.remv1], n.avg.misha.25.sr[index.25.anom.remv1],
       pch=18, col = "#00b4ffff")
points(n.avg.misha.25.dist[index.25.anom.remv2], n.avg.misha.25.sr[index.25.anom.remv2],
       pch=18, col = "#00b4ffff")

lines(n.avg.misha.25.dist[index.25.anom.remv1], n.avg.misha.25.sr[index.25.anom.remv1],
      lwd=1.5, col = "#00b4ffff")
lines(n.avg.misha.25.dist[index.25.anom.remv2], n.avg.misha.25.sr[index.25.anom.remv2],
      lwd=1.5, col = "#00b4ffff")