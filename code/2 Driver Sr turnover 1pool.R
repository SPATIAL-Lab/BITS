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
#######plotting 50 pt average#######
plot(misha.raw$dist, misha.raw$X87Sr.86Sr, col=alpha("black", 0.25),
     xlim=c(20000,8000),ylim=c(0.705,0.714),pch=16,main="50 pt average")#plot all raw data points
lines(n.avg.misha.50.dist, n.avg.misha.50.sr,col="#00b4ffff",lwd=1.5)
points(n.avg.misha.50.dist, n.avg.misha.50.sr, col="#00b4ffff",pch=18)

misha <- read.csv("data/Misha ivory.csv")

R.sd.mea <- misha$sd
dist.mea <- misha$dist
R.mea <- misha$mean
n.mea = length(dist.mea)

#adding the model component of food and water mixture as intake#
intake <- read.csv("data/intake.csv")

Hay <- subset(intake, intake$type=="H")
Pellet <- subset(intake, intake$type=="P")
Supplement <- subset(intake, intake$type=="S")
Water <- subset(intake, intake$type=="W")

#estimate mean and sd for hay Sr87/86
Sr.hay <- enorm(Hay$X87Sr.86Sr)
Sr.hay.mean <- Sr.hay$parameters[1]
Sr.hay.sd <- Sr.hay$parameters[2]

#estimate mean and sd for pellet Sr87/86
Sr.pel <- enorm(Pellet$X87Sr.86Sr)
Sr.pel.mean <- Sr.pel$parameters[1]
Sr.pel.sd <- Sr.pel$parameters[2]

#estimate mean and sd for alfalfa Sr87/86
Sr.sup <- enorm(Supplement$X87Sr.86Sr)
Sr.sup.mean <- Sr.sup$parameters[1]
Sr.sup.sd <- Sr.sup$parameters[2]

#estimate mean and sd for water Sr87/86
Sr.w <- enorm(Water$X87Sr.86Sr)
Sr.w.mean <- Sr.w$parameters[1]
Sr.w.sd <- Sr.w$parameters[2]

#estimate mean and sd for hay Sr concentration
#log-normal distribution
conc.hay <- elnorm(Hay$Sr_conc)
conc.hay.mean <- conc.hay$parameters[1]
conc.hay.sd <- conc.hay$parameters[2]

conc.sup <- elnorm(Supplement$Sr_conc)
conc.sup.mean <- conc.sup$parameters[1]
conc.sup.sd <- conc.sup$parameters[2]

#estimate mean and sd for pellet concentration
#log-normal distribution
conc.pel <- elnorm(Pellet$Sr_conc)
conc.pel.mean <- conc.pel$parameters[1]
conc.pel.sd <- conc.pel$parameters[2]

#estimate mean and sd for water concentration
#log-normal distribution
conc.w <- elnorm(Water$Sr_conc)
conc.w.mean <- conc.w$parameters[1]
conc.w.sd <- conc.w$parameters[2]

##########one pool turnover model########

#R0 is the mean ratio of initial value 
R0 <- 0.7071

m.feed <- 40 #mixing of x kg of hay, this is to constrain Sr variability of hay

Ivo.rate.mean <- 14.7
samp.interval <- mean(dist.mea[1:(n.mea - 1)] - dist.mea[2:n.mea])
n.days.bef.aft <- trunc(samp.interval/2/Ivo.rate.mean) #this parameter has to be supplied to the model

#parameters to save
parameters <- c("a", "Ivo.rate", "Rs.m","Rin","mod.index","mod.dist","f.h","f.pel","h.l",
                "Sr.pre", "Ps", "Fin", "Re.mean", "R0.mean", "switch", "w.contrib", "h.contrib")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 1050, n.mea = n.mea, 
           R0 = R0, Sr.hay.mean = Sr.hay.mean, Sr.hay.sd = Sr.hay.sd, 
           Sr.pel.mean = Sr.pel.mean, Sr.pel.sd = Sr.pel.sd, 
           Sr.sup.mean = Sr.sup.mean, Sr.sup.sd = Sr.sup.sd,
           Sr.w.mean = Sr.w.mean, Sr.w.sd = Sr.w.sd, m.feed = m.feed,
           conc.hay.mean = conc.hay.mean, conc.hay.sd = conc.hay.sd, 
           conc.pel.mean = conc.pel.mean, conc.pel.sd = conc.pel.sd,
           conc.sup.mean = conc.sup.mean, conc.sup.sd = conc.sup.sd,
           conc.w.mean = conc.w.mean, conc.w.sd = conc.w.sd,
           n.days.bef.aft = n.days.bef.aft)

#Start time
t1 = proc.time()

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 2e3
n.thin = floor(n.iter-n.burnin)/1000

#Run it
post.misha.nb5 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS no bone v5.R", 
                                            parameters.to.save = parameters, 
                                            data = dat, n.chains=3, n.iter = n.iter, 
                                            n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 4 hours

post.misha.nb5$BUGSoutput$summary

save(post.misha.nb5, file = "out/post.misha.nb5.RData")

load("out/post.misha.nb5.RData")

plot(density(post.misha.nb5$BUGSoutput$sims.list$a))

summary(post.misha.nb5$BUGSoutput$sims.list$switch)

plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rs.index.nb5 <- MCMC.ts.dist(post.misha.nb5$BUGSoutput$sims.list$Rs.m, 
                                     post.misha.nb5$BUGSoutput$sims.list$mod.dist,
                                     post.misha.nb5$BUGSoutput$sims.list$mod.index)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
MCMC.ts.RS.89.nb5 <- MCMC.ts(MCMC.ts.Rs.index.nb5)

lines(dist.mea, MCMC.ts.RS.89.nb5[[1]], lwd = 1, col = "cyan")
lines(dist.mea, MCMC.ts.RS.89.nb5[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, MCMC.ts.RS.89.nb5[[3]], lwd = 1, lty = 2, col = "cyan")

plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rin.index.nb5 <- MCMC.ts.dist(post.misha.nb5$BUGSoutput$sims.list$Rin, 
                                      post.misha.nb5$BUGSoutput$sims.list$mod.dist,
                                      post.misha.nb5$BUGSoutput$sims.list$mod.index)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
MCMC.ts.Rin.89.nb5 <- MCMC.ts(MCMC.ts.Rin.index.nb5)

lines(dist.mea, MCMC.ts.Rin.89.nb5[[1]], lwd = 1, col = "cyan")
lines(dist.mea, MCMC.ts.Rin.89.nb5[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, MCMC.ts.Rin.89.nb5[[3]], lwd = 1, lty = 2, col = "cyan")

plot(0,0, xlim = c(20000,8000), ylim = c(0, 1), xlab = "distance", ylab ="r")
MCMC.ts.fh.index.nb5 <- MCMC.ts.dist(post.misha.nb5$BUGSoutput$sims.list$f.h, 
                                     post.misha.nb5$BUGSoutput$sims.list$mod.dist,
                                     post.misha.nb5$BUGSoutput$sims.list$mod.index)

plot(0,0, xlim = c(20000,8000), ylim = c(0, 1), xlab = "distance", ylab ="r")
MCMC.ts.fpel.index.nb5 <- MCMC.ts.dist(post.misha.nb5$BUGSoutput$sims.list$f.pel, 
                                       post.misha.nb5$BUGSoutput$sims.list$mod.dist,
                                       post.misha.nb5$BUGSoutput$sims.list$mod.index)

plot(density(post.misha.nb5$BUGSoutput$sims.list$w.contrib))
plot(density(post.misha.nb5$BUGSoutput$sims.list$h.contrib))


###prep data###
R.sd.mea <- n.sd.misha.50.sr
dist.mea <- n.avg.misha.50.dist
R.mea <- n.avg.misha.50.sr
n.mea = length(n.avg.misha.50.sr)
##########one pool turnover model with 50 pt average full dataset########


#R0 is the mean ratio of initial value 
R0 <- 0.7071

m.feed <- 40 #mixing of x kg of hay, this is to constrain Sr variability of hay

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.6

s.intv <- mean(n.avg.misha.50.dist[1:(n.mea - 1)] - n.avg.misha.50.dist[2:n.mea])#52.4 microns

max.dist.mea <- max(n.avg.misha.50.dist) + 300 #adding some 20 days worth of distance prior to the start

#parameters to save
parameters <- c("a", "Ivo.rate", "Rs.m","Rin","mod.index","mod.dist","f.h","f.pel","h.l",
                "Sr.pre", "Ps", "Fin", "Re.mean", "R0.mean", "switch", "w.contrib", "h.contrib")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 1050, n.mea = n.mea, 
           R0 = R0, Sr.hay.mean = Sr.hay.mean, Sr.hay.sd = Sr.hay.sd, max.dist.mea = max.dist.mea,
           Sr.pel.mean = Sr.pel.mean, Sr.pel.sd = Sr.pel.sd, 
           Sr.sup.mean = Sr.sup.mean, Sr.sup.sd = Sr.sup.sd,
           Sr.w.mean = Sr.w.mean, Sr.w.sd = Sr.w.sd, m.feed = m.feed,
           conc.hay.mean = conc.hay.mean, conc.hay.sd = conc.hay.sd, 
           conc.pel.mean = conc.pel.mean, conc.pel.sd = conc.pel.sd,
           conc.sup.mean = conc.sup.mean, conc.sup.sd = conc.sup.sd,
           conc.w.mean = conc.w.mean, conc.w.sd = conc.w.sd,
           s.intv = s.intv, Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 2e3
n.thin = floor(n.iter-n.burnin)/400

#Run it
post.misha.fdnb = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS no bone v5.R", 
                                            parameters.to.save = parameters, 
                                            data = dat, n.chains=3, n.iter = n.iter, 
                                            n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 4 hours

post.misha.fdnb$BUGSoutput$summary

save(post.misha.fdnb, file = "out/post.misha.fdnb.RData")

load("out/post.misha.fdnb.RData")

plot(density(post.misha.fdnb$BUGSoutput$sims.list$a))

summary(post.misha.fdnb$BUGSoutput$sims.list$switch)

#check calibration
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.dist.plot(post.misha.fdnb$BUGSoutput$sims.list$Rs.m,
               post.misha.fdnb$BUGSoutput$sims.list$dist.m)
lines(n.avg.misha.50.dist,n.avg.misha.50.sr,lwd = 2, col = "red")


plot(density(post.misha.nb5$BUGSoutput$sims.list$w.contrib))
plot(density(post.misha.nb5$BUGSoutput$sims.list$h.contrib))

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


###prep data###
R.sd.mea <- n.sd.misha.25.sr[1:200]
dist.mea <- n.avg.misha.25.dist[1:200]
R.mea <- n.avg.misha.25.sr[1:200]
n.mea = length(n.avg.misha.25.sr[1:200])

R0 <- 0.707

#Re is the mean ratio of end value  
Re <- 0.7112
s.intv <- mean(dist.mea[1:(n.mea - 1)] - dist.mea[2:n.mea])

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.6
max.dist.mea <- max(n.avg.misha.25.dist) + 300
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
post.misha.fdnb1p = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS no bone woi.R", 
                                              parameters.to.save = parameters, 
                                              data = dat, n.chains=5, n.iter = n.iter, 
                                              n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 3.5 hours

post.misha.fdnb1p$BUGSoutput$summary

save(post.misha.fdnb1p, file = "out/post.misha.fdnb1p.RData")
load("out/post.misha.fdnb1p.RData")

traplot(post.misha.fdnb1p,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.fdnb1p,parms = c("a", "b","c"))

summary(post.misha.fdnb1p$BUGSoutput$sims.list$switch)
summary(post.misha.fdnb1p$BUGSoutput$sims.list$Ivo.rate)

#make a contour map
# contour.flux.pool <- kde2d(post.misha.5$BUGSoutput$sims.list$flux.ratio[,1], 
#                            post.misha.5$BUGSoutput$sims.list$pool.ratio[,1], n = 64)
# image(contour.flux.pool, col=viridis(64), xlab="Flux ratio: Fs/Fb", ylab="Pool ratio: Ps/Pb",xlim =c(0,200))
# contour(contour.flux.pool,lwd = 1.5, add = TRUE, labcex = 1)

#plotting some parameters
plot(density(post.misha.fdnb1p$BUGSoutput$sims.list$Ivo.rate))
map_estimate(post.misha.fdnb1p$BUGSoutput$sims.list$Ps) #0.66 mmol
map_estimate(post.misha.fdnb1p$BUGSoutput$sims.list$Fin) #0.02 mmol/day
hdi(post.misha.fdnb1p$BUGSoutput$sims.list$Fin, 0.89)[[2]]
hdi(post.misha.fdnb1p$BUGSoutput$sims.list$Fin, 0.89)[[3]]
plot(density(post.misha.fdnb1p$BUGSoutput$sims.list$Fin))
plot(density(post.misha.fdnb1p$BUGSoutput$sims.list$h.l)) #30 days

#plotting modeled serum values mapped onto ivory and checking the fit of the data
plot(0,0, xlim = c(20000,13000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.fdnb1p$BUGSoutput$sims.list$Rs.m,
               post.misha.fdnb1p$BUGSoutput$sims.list$dist)
lines(n.avg.misha.25.dist[1:200],n.avg.misha.25.sr[1:200],lwd = 2, col = "red")

#prior and posterior for R0 and Re
plot(density(post.misha.fdnb1p$BUGSoutput$sims.list$R0.mean),col="red",xlim=c(0.7065,0.7075))#posterior
lines(seq(0.706,0.708,0.00001),dnorm(seq(0.706,0.708,0.00001),mean = R0, sd= 1/sqrt(100/2e-6)))#prior

plot(density(post.misha.fdnb1p$BUGSoutput$sims.list$Re.mean),col="red",xlim=c(0.7095,0.7125))#posterior
lines(seq(0.7095,0.7125,0.00001),dnorm(seq(0.7095,0.7125,0.00001),mean = Re, sd= 1/sqrt(100/2e-5)))#prior

post.misha.fdnb1p.Rs.m.89<- MCMC.CI.bound(post.misha.fdnb1p$BUGSoutput$sims.list$Rs.m, 0.89)
#extract median distance from MCMC(t) results
med.dist.fdnb1p<- MCMC.dist.median(post.misha.fdnb1p$BUGSoutput$sims.list$dist)
lines(med.dist.fdnb1p, post.misha.fdnbh.Rs.m.89[[1]], lwd = 1, col = "cyan")
lines(med.dist.fdnb1p, post.misha.fdnbh.Rs.m.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(med.dist.fdnb1p, post.misha.fdnbh.Rs.m.89[[3]], lwd = 1, lty = 2, col = "cyan")

####check the Rin values####
plot(0,0, xlim = c(0,t), ylim = c(0.706, 0.713), xlab = "days", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
post.misha.fdnb1p.Rin.89<- MCMC.CI.bound(post.misha.fdnb1p$BUGSoutput$sims.list$Rin, 0.89)
lines(1:t, post.misha.fdnb1p.Rin.89[[1]], lwd = 2, col = "ff00ffff")
lines(1:t, post.misha.fdnb1p.Rin.89[[2]], lwd = 1, lty = 2, col = "ff00ffff")
lines(1:t, post.misha.fdnb1p.Rin.89[[3]], lwd = 1, lty = 2, col = "ff00ffff")

#reconstructed Rs.m history
plot(0,0, xlim = c(0,t), ylim = c(0.706, 0.713), xlab = "days", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.tl.plot(post.misha.fdnb1p$BUGSoutput$sims.list$Rs.m,750)
lines(1:750, post.misha.fdnb1p.Rs.m.89[[1]], lwd = 1, col = "cyan")
lines(1:750, post.misha.fdnb1p.Rs.m.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(1:750, post.misha.fdnb1p.Rs.m.89[[3]], lwd = 1, lty = 2, col = "cyan")

###########one pool model with anomaly removed#######
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
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.25.dist) + 30
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
post.misha.fdnb1pr = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS no bone woi.R", 
                                              parameters.to.save = parameters, 
                                              data = dat, n.chains=5, n.iter = n.iter, 
                                              n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 3.5 hours

post.misha.fdnb1pr$BUGSoutput$summary

save(post.misha.fdnb1pr, file = "out/post.misha.fdnb1pr.RData")
load("out/post.misha.fdnb1pr.RData")

traplot(post.misha.fdnb1pr,parms = c("a"))

summary(post.misha.fdnb1pr$BUGSoutput$sims.list$switch)
summary(post.misha.fdnb1pr$BUGSoutput$sims.list$Ivo.rate)

#plotting some parameters
plot(density(post.misha.fdnb1pr$BUGSoutput$sims.list$Ivo.rate))
map_estimate(post.misha.fdnb1pr$BUGSoutput$sims.list$Ps) #0.65 mmol
map_estimate(post.misha.fdnb1pr$BUGSoutput$sims.list$Fin) #0.02 mmol/day
plot(density(post.misha.fdnb1pr$BUGSoutput$sims.list$Fin))

a.param.nb1pr <- elnorm(post.misha.fdnb1pr$BUGSoutput$sims.list$a[,1])
#check data fit
plot(seq(0.015,0.045,0.0001), dlnorm(seq(0.015,0.045,0.0001),a.mean,a.sd), col = "blue", lwd = 2, type="l",
     xlim = c(0.015,0.045), xlab = "a", ylab= "density")
lines(density(post.misha.fdnb1pr$BUGSoutput$sims.list$a), col = "red", lwd = 2)





plot(density(post.misha.fdnb1p$BUGSoutput$sims.list$Fin),xlim= c(0.005, 0.035),ylim=c(0,180),lwd=2)
lines(density(post.misha.fdnb1pr$BUGSoutput$sims.list$Fin), col = "red",lwd=2)
lines(density(post.misha.fdnb1p50$BUGSoutput$sims.list$Fin), lwd=2, lty=2)
lines(density(post.misha.fdnbh$BUGSoutput$sims.list$Fin), lwd=2, lty=2, col = "blue")
lines(density(post.misha.fdnbhr$BUGSoutput$sims.list$Fin), lwd=2, lty=3, col = "blue")
lines(density(post.misha.fdnb2$BUGSoutput$sims.list$Fin), lwd=2, lty=3,col = "red")

legend(0.02, 180, c("25 pt avg with bone remodeling","25 pt avg w/o bone remodeling",
                     "50 pt avg with bone remodeling 1 pool", "50 pt avg with bone remodeling 2 pool"),
       lwd= c(2,2,2,2), col=c("black","red","black","red"), lty=c(1,1,2,3))


#plotting modeled serum values mapped onto ivory and checking the fit of the data
plot(0,0, xlim = c(20000,13000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.fdnb1pr$BUGSoutput$sims.list$Rs.m,
               post.misha.fdnb1pr$BUGSoutput$sims.list$dist)
lines(n.avg.misha.25.dist[index.25.anom.remv1], n.avg.misha.25.sr[index.25.anom.remv1],
      lwd = 2, col = "red")
lines(n.avg.misha.25.dist[index.25.anom.remv2], n.avg.misha.25.sr[index.25.anom.remv2],
      lwd = 2, col = "red")


post.misha.fdnbh.Rs.m.89<- MCMC.CI.bound(post.misha.fdnb1pr$BUGSoutput$sims.list$Rs.m, 0.89)
#extract median distance from MCMC(t) results
med.dist.fdnb1pr<- MCMC.dist.median(post.misha.fdnb1pr$BUGSoutput$sims.list$dist)
lines(med.dist.fdnb1pr, post.misha.fdnb1pr.Rs.m.89[[1]], lwd = 1, col = "cyan")
lines(med.dist.fdnb1pr, post.misha.fdnb1pr.Rs.m.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(med.dist.fdnb1pr, post.misha.fdnb1pr.Rs.m.89[[3]], lwd = 1, lty = 2, col = "cyan")

####check the Rin values####
plot(0,0, xlim = c(0,t), ylim = c(0.706, 0.713), xlab = "days", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
post.misha.fdnb1pr.Rin.89<- MCMC.CI.bound(post.misha.fdnb1pr$BUGSoutput$sims.list$Rin, 0.89)
lines(1:t, post.misha.fdnb1pr.Rin.89[[1]], lwd = 2, col = "ff00ffff")
lines(1:t, post.misha.fdnb1pr.Rin.89[[2]], lwd = 1, lty = 2, col = "ff00ffff")
lines(1:t, post.misha.fdnb1pr.Rin.89[[3]], lwd = 1, lty = 2, col = "ff00ffff")

#reconstructed Rs.m history
plot(0,0, xlim = c(0,t), ylim = c(0.706, 0.713), xlab = "days", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.tl.plot(post.misha.fdnb1pr$BUGSoutput$sims.list$Rs.m,750)
lines(1:750, post.misha.fdnb1pr.Rs.m.89[[1]], lwd = 1, col = "cyan")
lines(1:750, post.misha.fdnb1pr.Rs.m.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(1:750, post.misha.fdnb1pr.Rs.m.89[[3]], lwd = 1, lty = 2, col = "cyan")

##################50 pt average full dataset##############
misha.raw <- read.csv("data/Misha raw.csv")
misha.raw.Sr <- misha.raw$X87Sr.86Sr
misha.raw.dist <- misha.raw$dist

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

plot(n.avg.misha.50.dist, n.avg.misha.50.sr, main="50 pt average",
     xlim=c(max(n.avg.misha.50.dist),min(n.avg.misha.50.dist)))
lines(n.avg.misha.50.dist, n.avg.misha.50.sr)
hist(n.sd.misha.50.sr)
###prep data###
R.sd.mea <- n.sd.misha.50.sr
dist.mea <- n.avg.misha.50.dist
R.mea <- n.avg.misha.50.sr
n.mea = length(n.avg.misha.50.sr)

R0 <- 0.707

#Re is the mean ratio of end value  
Re <- 0.7112
s.intv <- mean(dist.mea[1:(n.mea - 1)] - dist.mea[2:n.mea])

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.25.dist) + 30
#parameters to save
parameters <- c("a", "Ivo.rate", "Rs.m","Rin","dist","h.l",
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

MCMC.dist.plot(post.misha.2p50$BUGSoutput$sims.list$Rs.m,
               post.misha.2p50$BUGSoutput$sims.list$dist)
lines(n.avg.misha.50.dist,n.avg.misha.50.sr,lwd = 2, col = "red")

#########1 pool with excursion removed
R.sd.mea <- n.sd.misha.50.sr
dist.mea <- n.avg.misha.50.dist
R.mea <- n.avg.misha.50.sr
n.mea = length(n.avg.misha.50.sr)

R0 <- 0.707

#Re is the mean ratio of end value  
Re <- 0.7112
s.intv <- mean(dist.mea[1:(n.mea - 1)] - dist.mea[2:n.mea])

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.25.dist) + 30
#parameters to save
parameters <- c("a", "Ivo.rate", "Rs.m","Rin","dist","h.l",
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
post.misha.fdnb1p50 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS no bone woi.R", 
                                                 parameters.to.save = parameters, 
                                                 data = dat, n.chains=5, n.iter = n.iter, 
                                                 n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 3.5 hours

post.misha.fdnb1p50$BUGSoutput$summary

save(post.misha.fdnb1p50, file = "out/post.misha.fdnb1p50.RData")
load("out/post.misha.fdnb1p50.RData")

traplot(post.misha.fdnb1p50,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.fdnb1p50,parms = c("a", "b","c"))

summary(post.misha.fdnb1p50$BUGSoutput$sims.list$switch)
summary(post.misha.fdnb1p50$BUGSoutput$sims.list$Ivo.rate)

#plotting some parameters
plot(density(post.misha.fdnb1p50$BUGSoutput$sims.list$Ivo.rate))
map_estimate(post.misha.fdnb1p50$BUGSoutput$sims.list$Ps) #0.66 mmol
map_estimate(post.misha.fdnb1p50$BUGSoutput$sims.list$Fin) #0.02 mmol/day
map_estimate(post.misha.fdnb1p50$BUGSoutput$sims.list$h.l) #36.7 day
hdi(post.misha.fdnb1p50$BUGSoutput$sims.list$Fin, 0.89)[[2]]
hdi(post.misha.fdnb1p50$BUGSoutput$sims.list$Fin, 0.89)[[3]]
plot(density(post.misha.fdnb1p50$BUGSoutput$sims.list$Fin))
plot(density(post.misha.fdnb1p50$BUGSoutput$sims.list$h.l))

#plotting modeled serum values mapped onto ivory and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.fdnb1p50$BUGSoutput$sims.list$Rs.m,
               post.misha.fdnb1p50$BUGSoutput$sims.list$dist)
lines(n.avg.misha.50.dist,n.avg.misha.50.sr,lwd = 2, col = "red")