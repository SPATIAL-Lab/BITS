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

####supplementary figures######

#one pool model fit#
misha.raw <- read.csv("data/Misha raw.csv")

#calculating moving average (400 micron band)
rm.misha.25 <- rollmean(misha.raw, k = 25)

plot(rm.misha.25[,1],rm.misha.25[,2], type="l", main="rolling mean 25",
     xlim=c(max(rm.misha.25[,1]),min(rm.misha.25[,1])))

rm.misha.100 <- rollmean(misha.raw, k = 100)

plot(rm.misha.100[,1],rm.misha.100[,2], type="l", main="rolling mean 100",
     xlim=c(max(rm.misha.100[,1]),min(rm.misha.100[,1])))

##calculate 100 point average to reduce data amount
misha.raw.Sr <- misha.raw$X87Sr.86Sr
misha.raw.dist <- misha.raw$dist

n.avg <- 100 #now many data points to average
n <- floor(length(misha.raw$X87Sr.86Sr)/n.avg) #99 data points, discard the last <100 data points

#initiate vectors
n.avg.misha.100.sr <- rep(0, n)
n.sd.misha.100.sr <- rep(0, n)

n.avg.misha.100.dist <- rep(0, n)

for(i in 1:n){
  x <- ((i-1)*n.avg + 1):(i*n.avg)
  n.avg.misha.100.sr[i] <- mean(misha.raw.Sr[x])
  n.sd.misha.100.sr[i] <- sd(misha.raw.Sr[x])
  
  n.avg.misha.100.dist[i] <- mean(misha.raw.dist[x])#mean and median are the same
}

plot(misha.raw$dist, misha.raw$X87Sr.86Sr, col=alpha("dark blue", 0.2),
     xlim=c(20000,8000),ylim=c(0.705,0.715),pch=16,main="100 pt average")#plot all raw data points
lines(n.avg.misha.100.dist, n.avg.misha.100.sr,col="#00b4ffff",lwd=2)
points(n.avg.misha.100.dist, n.avg.misha.100.sr, col="#00b4ffff",pch=18)

##calculate 50 point average to reduce data amount
n.avg <- 50 #now many data points to average
n <- floor(length(misha.raw$X87Sr.86Sr)/n.avg) #99 data points, discard the last <100 data points

#initiate vectors
n.avg.misha.50.sr <- rep(0, n)
n.sd.misha.50.sr <- rep(0, n)

n.avg.misha.50.dist <- rep(0, n)

for(i in 1:n){
  x <- ((i-1)*n.avg + 1):(i*n.avg)
  n.avg.misha.50.sr[i] <- mean(misha.raw.Sr[x])
  n.sd.misha.50.sr[i] <- sd(misha.raw.Sr[x])
  
  n.avg.misha.50.dist[i] <- mean(misha.raw.dist[x])#mean and median are the same
}

#check data range
plot(n.avg.misha.50.dist, n.avg.misha.50.sr, main="50 pt average",
     xlim=c(max(n.avg.misha.50.dist),min(n.avg.misha.50.dist)))
lines(n.avg.misha.50.dist, n.avg.misha.50.sr)
hist(n.sd.misha.50.sr)

#########1 pool turnover model with excursion removed (50 pt average dataset)###################
ind.remv.50 <- which(n.avg.misha.50.dist > 17000 & n.avg.misha.50.dist < 17500 & n.avg.misha.50.sr < 0.709)

index.50.anom.remv1 <- c(1:35)
index.50.anom.remv2 <- c(39:length(n.avg.misha.50.dist))
index.50.anom.remv <- c(1:35,39:length(n.avg.misha.50.dist))
n.avg.misha.50.dist.rmv <- n.avg.misha.50.dist[index.50.anom.remv]
n.avg.misha.50.sr.rmv <- n.avg.misha.50.sr[index.50.anom.remv]
n.sd.misha.50.sr.rmv <- n.sd.misha.50.sr[index.50.anom.remv]

R.sd.mea <- n.sd.misha.50.sr.rmv
dist.mea <- n.avg.misha.50.dist.rmv
R.mea <- n.avg.misha.50.sr.rmv
n.mea = length(n.avg.misha.50.sr.rmv)

R0 <- 0.707

#Re is the mean ratio of end value  
Re <- 0.711
s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.25.dist) + 30
#parameters to save
parameters <- c("a", "Ivo.rate", "rs.m","Rin","dist","h.l",
                "Sr.pre", "Ps", "Fin", "Re.mean", "R0.mean", "switch")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 750, n.mea = n.mea, 
           R0 = R0, Re = Re, s.intv = s.intv, Ivo.rate.mean = Ivo.rate.mean,
           Ivo.rate.sd = Ivo.rate.sd, max.dist.mea = max.dist.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/400

#Run it
post.misha.1p50r = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS no bone woi.R", 
                                              parameters.to.save = parameters, 
                                              data = dat, n.chains=5, n.iter = n.iter, 
                                              n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 3.5 hours

post.misha.1p50r$BUGSoutput$summary

save(post.misha.1p50r, file = "out/post.misha.1p50r.RData")
load("out/post.misha.1p50r.RData")

traplot(post.misha.1p50r,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.1p50r,parms = c("a", "b","c"))

summary(post.misha.1p50r$BUGSoutput$sims.list$switch)
summary(post.misha.1p50r$BUGSoutput$sims.list$Ivo.rate)

#plotting modeled serum values mapped onto ivory and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.1p50r$BUGSoutput$sims.list$rs.m,
               post.misha.1p50r$BUGSoutput$sims.list$dist)
lines(n.avg.misha.50.dist,n.avg.misha.50.sr,lwd = 2, col = "red")

#########1 pool turnover model with 50 pt average dataset###################
###prep data###
R.sd.mea <- n.sd.misha.50.sr
dist.mea <- n.avg.misha.50.dist
R.mea <- n.avg.misha.50.sr
n.mea = length(n.avg.misha.50.sr)

R0 <- 0.707

#Re is the mean ratio of end value  
Re <- 0.711
s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.25.dist) + 30
#parameters to save
parameters <- c("a", "Ivo.rate", "rs.m","Rin","dist","h.l",
                "Sr.pre", "Ps", "Fin", "Re.mean", "R0.mean", "switch")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 750, n.mea = n.mea, 
           R0 = R0, Re = Re, s.intv = s.intv, Ivo.rate.mean = Ivo.rate.mean,
           Ivo.rate.sd = Ivo.rate.sd, max.dist.mea = max.dist.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/400

#Run it
post.misha.2p50 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS no bone woi.R", 
                                             parameters.to.save = parameters, 
                                             data = dat, n.chains=5, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 3.5 hours

post.misha.2p50$BUGSoutput$summary

save(post.misha.2p50, file = "out/post.misha.2p50.RData")
load("out/post.misha.2p50.RData")

traplot(post.misha.2p50,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.2p50,parms = c("a"))

summary(post.misha.2p50$BUGSoutput$sims.list$switch)
summary(post.misha.2p50$BUGSoutput$sims.list$Ivo.rate)

#plotting some parameters
plot(density(post.misha.2p50$BUGSoutput$sims.list$Ivo.rate))
map_estimate(post.misha.2p50$BUGSoutput$sims.list$Ps) #0.66 mmol
map_estimate(post.misha.2p50$BUGSoutput$sims.list$Fin) #0.02 mmol/day
map_estimate(post.misha.2p50$BUGSoutput$sims.list$h.l) #36.7 day
hdi(post.misha.2p50$BUGSoutput$sims.list$Fin, 0.89)[[2]]
hdi(post.misha.2p50$BUGSoutput$sims.list$Fin, 0.89)[[3]]
plot(density(post.misha.2p50$BUGSoutput$sims.list$Fin))
plot(density(post.misha.2p50$BUGSoutput$sims.list$h.l))

#plotting modeled serum values mapped onto ivory and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.2p50$BUGSoutput$sims.list$rs.m,
               post.misha.2p50$BUGSoutput$sims.list$dist)
lines(n.avg.misha.50.dist,n.avg.misha.50.sr,lwd = 2, col = "red")

#######Fig S1 Position of cracks and how they have little effect on measured Sr 87/86######

#######Fig S2 Log-normal goodness of fit for parameter a#######
#check q-q plot for fit
qqPlot(post.misha.fdnb1pr$BUGSoutput$sims.list$a[,1],distribution = "lnorm",
       param.list=list(mean=a.param.nb1pr$parameters[1],sd=a.param.nb1pr$parameters[2]),add.line=T)

#######Fig S3 results of intake model###########
#use a different JAGS file for evaluation!

#######Fig S4 sensitivity test 1 Sensitivity to excursion
#excursion makes rate parameter smaller
#########1 pool turnover model with full 50 pt average dataset###################
###prep data###
R.sd.mea <- n.sd.misha.50.sr
dist.mea <- n.avg.misha.50.dist
R.mea <- n.avg.misha.50.sr
n.mea = length(n.avg.misha.50.sr)

R0 <- 0.707

#Re is the mean ratio of end value  
Re <- 0.711
s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30
#parameters to save
parameters <- c("a", "Ivo.rate", "r1.m","Rin","dist","h.l",
                "Sr.pre", "Re.mean", "R0.mean", "switch")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 750, n.mea = n.mea, 
           R0 = R0, Re = Re, s.intv = s.intv, Ivo.rate.mean = Ivo.rate.mean,
           Ivo.rate.sd = Ivo.rate.sd, max.dist.mea = max.dist.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/400

#Run it
post.misha.2p50 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS no bone woi.R", 
                                             parameters.to.save = parameters, 
                                             data = dat, n.chains=5, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 3.5 hours

post.misha.2p50$BUGSoutput$summary

save(post.misha.2p50, file = "out/post.misha.2p50.RData")
load("out/post.misha.2p50.RData")

traplot(post.misha.2p50,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.2p50,parms = c("a"))

summary(post.misha.2p50$BUGSoutput$sims.list$switch)
summary(post.misha.2p50$BUGSoutput$sims.list$Ivo.rate)

#plotting some parameters
plot(density(post.misha.2p50$BUGSoutput$sims.list$Ivo.rate))
map_estimate(post.misha.2p50$BUGSoutput$sims.list$h.l)
plot(density(post.misha.2p50$BUGSoutput$sims.list$h.l))

#plotting modeled serum values mapped onto ivory and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.2p50$BUGSoutput$sims.list$R1.m,
               post.misha.2p50$BUGSoutput$sims.list$dist)
lines(n.avg.misha.50.dist,n.avg.misha.50.sr,lwd = 2, col = "red")

#######Fig S5 sensitivity test 2 Sensitivity to switch date
#how do you determine switch date?

######Fig S6: Sensitivity test 3 Sensitivity to precision terms?


