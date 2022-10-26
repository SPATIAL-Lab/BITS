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

########inversion model without calibration, parameters taken from the two pool model#######
#inversion based on micromill results

misha <- read.csv("data/Misha ivory.csv")

R.sd.cal <- misha$sd
dist.cal <- misha$dist
R.cal <- misha$mean
n.cal = length(dist.cal)

micromill <- read.csv("data/Misha micromill.csv")

R.mic <- rev(micromill$X87Sr.86Sr)
R.sd.mic <- rev(micromill$Sd)
dist.mic <- rev(micromill$dist)
n.mic <- length(dist.mic)

#assign input values before and after the switch, but allows some variation
R0 <- 0.70706

#Re is the mean ratio of end value  
Re <- 0.7112

max.dist.cal <- 19200

max.dist.mea <- 12400

s.intv <- 400

cal.intv <- 100

#micromill rate
Ivo.rate.mean <- 19.9 #microns per day
Ivo.rate.sd <- 5.4

#calibration rate
Ivo.rate.cal.mean <- 14.7 #microns per day
Ivo.rate.cal.sd <- 0.6

#parameters to save
parameters <- c("Ivo.rate", "Rs.cal", "Rb.cal", "dist.cal.m","dist","Rin.cal","a","b","c",
                "Rs.m", "Rb.m","Rin.m", "a.m", "b.m", "c.m","Rin.m.cps.ac","Ivo.rate.cal")

##Data to pass to the model
#compared to the turnover model that is essentially the .cal part here 
#the inversion takes the measured value of potentially a different ivory series
dat = list(R.cal = R.cal, dist.cal = dist.cal, R.sd.cal = R.sd.cal, t.cal = 750, n.cal = n.cal, 
           R0 = R0, Re = Re, cal.intv = cal.intv, s.intv = s.intv,
           max.dist.cal = max.dist.cal, max.dist.mea = max.dist.mea,
           Ivo.rate.cal.mean = Ivo.rate.cal.mean, Ivo.rate.cal.sd = Ivo.rate.cal.sd,
           Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
           R.mea = R.mic, dist.mea = dist.mic, R.sd.mea = R.sd.mic, t = 700, n.mea = n.mic)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 2e3
n.thin = floor(n.iter-n.burnin)/400

#Run it
post.mic.inv.woi = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS wo intake model.R", 
                                                   parameters.to.save = parameters, 
                                                   data = dat, n.chains=5, n.iter = n.iter, 
                                                   n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 25 hours

save(post.mic.inv.woi, file = "out/post.mic.inv.woi.RData")

post.mic.inv.woi$BUGSoutput$summary

load("out/post.mic.inv.woi.RData")

plot(density(post.mic.inv.woi$BUGSoutput$sims.list$a))

#plot calibration curve
par(mfrow=c(1,1))
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.mic.inv.woi$BUGSoutput$sims.list$Rs.cal,
               post.mic.inv.woi$BUGSoutput$sims.list$dist.cal.m)
lines(misha$dist,misha$mean,lwd = 2, col = "red")

post.mic.inversion.woi.Rs.cal.89<- MCMC.CI.bound(post.mic.inv.woi$BUGSoutput$sims.list$Rs.cal, 0.89)

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,500), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to days using rate Ivo.rate.mic
lines((max.dist.mea - dist.mic)/Ivo.rate.mic + 1, R.mic, lwd = 2, col = "blue")#results from micromill

MCMC.ts.Rin.m.89.woi <- MCMC.ts(post.mic.inv.woi$BUGSoutput$sims.list$Rin.m)
lines(1:500,MCMC.ts.Rin.m.89.woi[[1]],lwd = 2, col = "firebrick4")
lines(1:500,MCMC.ts.Rin.m.89.woi[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:500,MCMC.ts.Rin.m.89.woi[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))

MCMC.ts.Rs.m.89.woi <- MCMC.ts(post.mic.inversion.woi$BUGSoutput$sims.list$Rs.m)
MCMC.ts.Rb.m.89.woi <- MCMC.ts(post.mic.inversion.woi$BUGSoutput$sims.list$Rb.m)

plot(0,0, xlim = c(1,t), ylim = c(0.705, 0.716), xlab = "days", ylab ="ivory rate")
lines((max.dist.mea - dist.mic)/Ivo.rate.mic + 1, R.mic, lwd = 2, col = "blue")#results from micromill
lines(1:t,MCMC.ts.Rs.m.89.woi[[1]],lwd = 2, col = "firebrick4")
lines(1:t,MCMC.ts.Rs.m.89.woi[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:t,MCMC.ts.Rs.m.89.woi[[3]], lwd = 1, lty = 2, col = "firebrick4")
points((max.dist.mea - dist.mic)/Ivo.rate.mic + 1, R.mic)

MCMC.ts.Ivo.rate.89.woi <- MCMC.ts(post.mic.inversion.woi$BUGSoutput$sims.list$Ivo.rate)

plot(0,0, xlim = c(1,t), ylim = c(10, 25), xlab = "days", ylab ="ivory rate")
lines(1:t,MCMC.ts.Ivo.rate.89.woi[[1]],lwd = 2, col = "firebrick4")

par(mfrow=c(1,3))

plot(abc.prior.params$x,abc.prior.params$y[1,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "a", ylab= "density")
lines(density(log(post.mic.inversion.woi$BUGSoutput$sims.list$a.m)), col = "red", lwd = 2)

plot(abc.prior.params$x,abc.prior.params$y[2,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "b", ylab= "density")
lines(density(log(post.mic.inversion.woi$BUGSoutput$sims.list$b.m)), col = "red", lwd = 2)

plot(abc.prior.params$x,abc.prior.params$y[3,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "c", ylab= "density")
lines(density(log(post.mic.inversion.woi$BUGSoutput$sims.list$c.m)), col = "red", lwd = 2)

########inversion model with calibration#######
#inversion based on laser ablation results

misha <- read.csv("data/Misha ivory.csv")

R.sd.cal <- misha$sd
dist.cal <- misha$dist
R.cal <- misha$mean
n.cal = length(dist.cal)

#assign input values before and after the switch, but allows some variation
R0 <- 0.70706

#Re is the mean ratio of end value  
Re <- 0.7112

cal.intv <- 100

max.dist.cal <- 19200

#calibration rate
Ivo.rate.cal.mean <- 14.7 #microns per day
Ivo.rate.cal.sd <- 0.6

#parameters to save
parameters <- c("Ivo.rate", "Rs.cal", "Rb.cal", "dist.cal.m","dist","Rin.cal","a","b","c",
                "Rs.m", "Rb.m","Rin.m", "a.m", "b.m", "c.m","Rin.m.cps.ac","Ivo.rate.cal")

##Data to pass to the model
#compared to the turnover model that is essentially the .cal part here 
#the inversion takes the measured value of potentially a different ivory series
dat = list(R.cal = R.cal, dist.cal = dist.cal, R.sd.cal = R.sd.cal, t.cal = 750, n.cal = n.cal, 
           R0 = R0, Re = Re, cal.intv = cal.intv, s.intv = cal.intv,
           max.dist.cal = max.dist.cal, max.dist.mea = max.dist.cal,
           Ivo.rate.cal.mean = Ivo.rate.cal.mean, Ivo.rate.cal.sd = Ivo.rate.cal.sd,
           Ivo.rate.mean = Ivo.rate.cal.mean, Ivo.rate.sd = Ivo.rate.cal.sd,
           R.mea = R.cal, dist.mea = dist.cal, R.sd.mea = R.sd.cal, t = 750, n.mea = n.cal)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 4e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/600

#Run it
post.misha.inversion.woi = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS wo intake model.R", 
                                                      parameters.to.save = parameters, 
                                                      data = dat, n.chains=5, n.iter = n.iter, 
                                                      n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 19 hours

save(post.misha.inversion.woi, file = "out/post.misha.inversion.woi.RData")

post.misha.inversion.woi$BUGSoutput$summary

load("out/post.misha.inversion.woi.RData")
par(mfrow=c(1,3))
plot(abc.prior.params$x,abc.prior.params$y[1,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "a", ylab= "density")
lines(density(log(post.misha.inversion.woi$BUGSoutput$sims.list$a.m)), col = "red", lwd = 2)

plot(abc.prior.params$x,abc.prior.params$y[2,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "b", ylab= "density")
lines(density(log(post.misha.inversion.woi$BUGSoutput$sims.list$b.m)), col = "red", lwd = 2)

plot(abc.prior.params$x,abc.prior.params$y[3,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "c", ylab= "density")
lines(density(log(post.misha.inversion.woi$BUGSoutput$sims.list$c.m)), col = "red", lwd = 2)

plot(density(post.misha.inversion.woi$BUGSoutput$sims.list$a))
plot(density(post.misha.inversion.woi$BUGSoutput$sims.list$b))
plot(density(post.misha.inversion.woi$BUGSoutput$sims.list$c))

#plot calibration curve
par(mfrow=c(1,1))
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.inversion.woi$BUGSoutput$sims.list$Rs.cal,
               post.misha.inversion.woi$BUGSoutput$sims.list$dist.cal.m)
lines(misha$dist,misha$mean,lwd = 2, col = "red")

post.misha.inversion.woi.Rs.cal.89<- MCMC.CI.bound(post.misha.inversion.woi$BUGSoutput$sims.list$Rs.cal, 0.89)

med.dist.inversion.woi<- MCMC.dist.median(post.misha.inversion.woi$BUGSoutput$sims.list$dist)
lines(med.dist.inversion.woi, post.misha.inversion.woi.Rs.cal.89[[1]], lwd = 1, col = "cyan")
lines(med.dist.inversion.woi, post.misha.inversion.woi.Rs.cal.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(med.dist.inversion.woi, post.misha.inversion.woi.Rs.cal.89[[3]], lwd = 1, lty = 2, col = "cyan")

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,750), ylim = c(0.706, 0.715), xlab = "days", ylab ="Sr 87/86")
points((max(misha$dist)-misha$dist)/14.7,misha$mean, pch= 18)
#converting misha distance to days using rate Ivo.rate
MCMC.ts.Rin.cal.89.woi <- MCMC.ts(post.misha.inversion.woi$BUGSoutput$sims.list$Rin.cal)
lines(1:750,MCMC.ts.Rin.cal.89.woi[[1]],lwd = 2, col = "blue")
lines(1:750,MCMC.ts.Rin.cal.89.woi[[2]], lwd = 1, lty = 2, col = "blue")
lines(1:750,MCMC.ts.Rin.cal.89.woi[[3]], lwd = 1, lty = 2, col = "blue")

MCMC.ts.Rin.m.89.woi <- MCMC.ts(post.misha.inversion.woi$BUGSoutput$sims.list$Rin.m)
lines(1:750,MCMC.ts.Rin.m.89.woi[[1]],lwd = 2, col = "firebrick4")
lines(1:750,MCMC.ts.Rin.m.89.woi[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:750,MCMC.ts.Rin.m.89.woi[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.715, c("Model input","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))


#########################inversion 1 pool misha LA##########
########inversion model with calibration#######
#inversion based on laser ablation results
R.sd.cal <- n.sd.misha.25.sr.rmv
dist.cal <- n.avg.misha.25.dist.rmv
R.cal <- n.avg.misha.25.sr.rmv
n.cal = length(n.avg.misha.25.sr.rmv)

R0 <- 0.707

#Re is the mean ratio of end value  
Re <- 0.711
cal.intv <- mean(dist.cal[1:(n.cal - 1)] - dist.cal[2:n.cal])

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.6

max.dist.cal <- max(n.avg.misha.25.dist) + 300

R.sd.mea <- n.sd.misha.50.sr
dist.mea <- n.avg.misha.50.dist
R.mea <- n.avg.misha.50.sr
n.mea = length(n.avg.misha.50.sr)

s.intv <- mean(dist.mea[1:(n.mea - 1)] - dist.mea[2:n.mea])

max.dist.mea <- max(n.avg.misha.50.dist) + 300

#parameters to save
parameters <- c("Ivo.rate", "Rs.cal", "dist.cal.m","dist","Rin.cal","a",
                "Rs.m","Rin.m", "a.m","Ivo.rate.cal")

##Data to pass to the model
#compared to the turnover model that is essentially the .cal part here 
#the inversion takes the measured value of potentially a different ivory series
dat = list(R.cal = R.cal, dist.cal = dist.cal, R.sd.cal = R.sd.cal, t.cal = 400, n.cal = n.cal, 
           R0 = R0, Re = Re, cal.intv = cal.intv, s.intv = s.intv,
           max.dist.cal = max.dist.cal, max.dist.mea = max.dist.mea,
           Ivo.rate.cal.mean = Ivo.rate.mean, Ivo.rate.cal.sd = Ivo.rate.sd,
           Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
           R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 780, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/400

#Run it
post.misha.inv1p.woi = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS 1p woi.R", 
                                                      parameters.to.save = parameters, 
                                                      data = dat, n.chains=5, n.iter = n.iter, 
                                                      n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 19 hours

save(post.misha.inv1p.woi, file = "out/post.misha.inv1p.woi.RData")

post.misha.inv1p.woi$BUGSoutput$summary

load("out/post.misha.inv1p.woi.RData")

plot(density(post.misha.inversion.woi$BUGSoutput$sims.list$a))

#plot calibration curve
par(mfrow=c(1,1))
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.inv1p.woi$BUGSoutput$sims.list$Rs.cal,
               post.misha.inv1p.woi$BUGSoutput$sims.list$dist.cal.m)
lines(misha$dist,misha$mean,lwd = 2, col = "red")

post.misha.inv1p.woi.Rs.cal.89<- MCMC.CI.bound(post.misha.inv1p.woi$BUGSoutput$sims.list$Rs.cal, 0.89)

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,750), ylim = c(0.706, 0.715), xlab = "days", ylab ="Sr 87/86")
points((max(misha$dist)-misha$dist)/14.7,misha$mean, pch= 18)
#converting misha distance to days using rate Ivo.rate
MCMC.ts.Rin.cal.89.woi <- MCMC.ts(post.misha.inversion.woi$BUGSoutput$sims.list$Rin.cal)
lines(1:750,MCMC.ts.Rin.cal.89.woi[[1]],lwd = 2, col = "blue")
lines(1:750,MCMC.ts.Rin.cal.89.woi[[2]], lwd = 1, lty = 2, col = "blue")
lines(1:750,MCMC.ts.Rin.cal.89.woi[[3]], lwd = 1, lty = 2, col = "blue")

MCMC.ts.Rin.m.89.woi <- MCMC.ts(post.misha.inversion.woi$BUGSoutput$sims.list$Rin.m)
lines(1:750,MCMC.ts.Rin.m.89.woi[[1]],lwd = 2, col = "firebrick4")
lines(1:750,MCMC.ts.Rin.m.89.woi[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:750,MCMC.ts.Rin.m.89.woi[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.715, c("Model input","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))

###################Inversion with prior distribution of parameter a########

#########################inversion 1 pool misha LA##########
########inversion model with calibration#######
#inversion based on laser ablation results

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.6

R.sd.mea <- n.sd.misha.50.sr
dist.mea <- n.avg.misha.50.dist
R.mea <- n.avg.misha.50.sr
n.mea = length(n.avg.misha.50.sr)

s.intv <- mean(dist.mea[1:(n.mea - 1)] - dist.mea[2:n.mea])

max.dist.mea <- max(n.avg.misha.50.dist) + 30

a.mean <- a.param.nb1pr$parameters[1]
a.sd <- a.param.nb1pr$parameters[2]
#parameters to save
parameters <- c("Ivo.rate", "Rs.cal", "dist.cal.m","dist","Rin.cal","a",
                "Rs.m","Rin.m", "a.m")

##Data to pass to the model
#compared to the turnover model that is essentially the .cal part here 
#the inversion takes the measured value of potentially a different ivory series
dat = list( s.intv = s.intv, max.dist.mea = max.dist.mea, a.mean=a.mean, a.sd=a.sd,
           Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
           R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 750, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/400

#Run it
post.misha.inv1p.param = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS 1p param.R", 
                                                  parameters.to.save = parameters, 
                                                  data = dat, n.chains=5, n.iter = n.iter, 
                                                  n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 12 hours

save(post.misha.inv1p.param, file = "out/post.misha.inv1p.param.RData")

post.misha.inv1p.param$BUGSoutput$summary

load("out/post.misha.inv1p.param.RData")

plot(density(post.misha.inv1p.param$BUGSoutput$sims.list$Ivo.rate))

#check prior vs posterior parameters
plot(seq(0.015,0.045,0.0001), dlnorm(seq(0.015,0.045,0.0001),a.mean,a.sd), col = "blue", lwd = 2, type="l",
     xlim = c(0.015,0.045), xlab = "a", ylab= "density")
lines(density(post.misha.fdnb1pr$BUGSoutput$sims.list$a), col = "red", lwd = 2)
lines(density(post.misha.inv1p.param$BUGSoutput$sims.list$a.m), col = "red", lwd = 2,lty=2)
#slight deviation from prior

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,750), ylim = c(0.706, 0.715), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(max.dist.mea)-n.avg.misha.50.dist)/14.7,n.avg.misha.50.sr, pch= 18, col="#00b3ffff")
#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.fdnb1pr.Rin.89<- MCMC.CI.bound(post.misha.fdnb1pr$BUGSoutput$sims.list$Rin, 0.89)
lines(1:400,post.misha.fdnb1pr.Rin.89[[1]],lwd = 2, col = "black")
lines(1:400,post.misha.fdnb1pr.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:400,post.misha.fdnb1pr.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

#estimated input series
MCMC.ts.Rin.m.inv1p.param.89<- MCMC.CI.bound(post.misha.inv1p.param$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:750,MCMC.ts.Rin.m.inv1p.param.89[[1]],lwd = 2, col = "magenta")
lines(1:750,MCMC.ts.Rin.m.inv1p.param.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:750,MCMC.ts.Rin.m.inv1p.param.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Model input","Reconstructed input"),lwd = c(2, 2), col=c("black","magenta"))

###################Inversion with prior distribution of parameter a########

#########################inversion 1 pool misha LA##########
########inversion model with calibration#######
#inversion based on laser ablation results 25 pt average and first half of the data

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.6

R.sd.mea <- n.sd.misha.25.sr[1:200]
dist.mea <- n.avg.misha.25.dist[1:200]
R.mea <- n.avg.misha.25.sr[1:200]
n.mea = length(n.avg.misha.25.sr[1:200])

s.intv <- 26.2

max.dist.mea <- max(n.avg.misha.25.dist) + 30

a.mean <- a.param.nb1pr$parameters[1]
a.sd <- a.param.nb1pr$parameters[2]
#parameters to save
parameters <- c("Ivo.rate", "Rs.cal", "dist.cal.m","dist","Rin.cal","a",
                "Rs.m","Rin.m", "a.m")

##Data to pass to the model
#compared to the turnover model that is essentially the .cal part here 
#the inversion takes the measured value of potentially a different ivory series
dat = list( s.intv = s.intv, max.dist.mea = max.dist.mea, a.mean=a.mean, a.sd=a.sd,
            Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 400, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/400

#Run it
post.misha.inv1ph.param = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS 1p param.R", 
                                                    parameters.to.save = parameters, 
                                                    data = dat, n.chains=5, n.iter = n.iter, 
                                                    n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 12 hours

save(post.misha.inv1ph.param, file = "out/post.misha.inv1ph.param.RData")

post.misha.inv1ph.param$BUGSoutput$summary

load("out/post.misha.inv1ph.param.RData")

plot(density(post.misha.inv1ph.param$BUGSoutput$sims.list$Ivo.rate))

#check prior vs posterior parameters
plot(seq(0.015,0.045,0.0001), dlnorm(seq(0.015,0.045,0.0001),a.mean,a.sd), col = "blue", lwd = 2, type="l",
     xlim = c(0.015,0.045), xlab = "a", ylab= "density")
lines(density(post.misha.fdnb1pr$BUGSoutput$sims.list$a), col = "red", lwd = 2)
lines(density(post.misha.inv1ph.param$BUGSoutput$sims.list$a.m), col = "red", lwd = 2,lty=2)
#slight deviation from prior

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,400), ylim = c(0.706, 0.715), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(max.dist.mea)-n.avg.misha.25.dist[1:200])/14.7,n.avg.misha.25.sr[1:200], pch= 18, col="#00b3ffff")
#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.fdnb1pr.Rin.89<- MCMC.CI.bound(post.misha.fdnb1pr$BUGSoutput$sims.list$Rin, 0.89)
lines(1:400,post.misha.fdnb1pr.Rin.89[[1]],lwd = 2, col = "black")
lines(1:400,post.misha.fdnb1pr.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:400,post.misha.fdnb1pr.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

#estimated input series
MCMC.ts.Rin.m.inv1ph.param.89<- MCMC.CI.bound(post.misha.inv1ph.param$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:400,MCMC.ts.Rin.m.inv1ph.param.89[[1]],lwd = 2, col = "magenta")
lines(1:400,MCMC.ts.Rin.m.inv1ph.param.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:400,MCMC.ts.Rin.m.inv1ph.param.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Model input","Reconstructed input"),lwd = c(2, 2), col=c("black","magenta"))

###################Inversion with prior distribution of parameter a########

#########################inversion 1 pool misha LA##########
########inversion model with calibration#######
#inversion based on laser ablation results 25 pt average and first half of the data with excursion removed

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.6

R.sd.mea <- n.sd.misha.25.sr.rmv
dist.mea <- n.avg.misha.25.dist.rmv
R.mea <- n.avg.misha.25.sr.rmv
n.mea = length(n.avg.misha.25.sr.rmv)

s.intv <- 26.2

max.dist.mea <- max(n.avg.misha.25.dist) + 30

a.mean <- a.param.nb1pr$parameters[1]
a.sd <- a.param.nb1pr$parameters[2]
#parameters to save
parameters <- c("Ivo.rate", "Rs.cal", "dist.cal.m","dist","Rin.cal","a",
                "Rs.m","Rin.m", "a.m")

##Data to pass to the model
#compared to the turnover model that is essentially the .cal part here 
#the inversion takes the measured value of potentially a different ivory series
dat = list( s.intv = s.intv, max.dist.mea = max.dist.mea, a.mean=a.mean, a.sd=a.sd,
            Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 400, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/400

#Run it
post.misha.inv1phr.param = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS 1p param.R", 
                                                     parameters.to.save = parameters, 
                                                     data = dat, n.chains=5, n.iter = n.iter, 
                                                     n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 12 hours

save(post.misha.inv1phr.param, file = "out/post.misha.inv1phr.param.RData")

post.misha.inv1phr.param$BUGSoutput$summary

load("out/post.misha.inv1phr.param.RData")

plot(density(post.misha.inv1phr.param$BUGSoutput$sims.list$Ivo.rate))

#check prior vs posterior parameters
plot(seq(0.015,0.045,0.0001), dlnorm(seq(0.015,0.045,0.0001),a.mean,a.sd), col = "blue", lwd = 2, type="l",
     xlim = c(0.015,0.045), xlab = "a", ylab= "density")
lines(density(post.misha.fdnb1pr$BUGSoutput$sims.list$a), col = "red", lwd = 2)
lines(density(post.misha.inv1phr.param$BUGSoutput$sims.list$a.m), col = "red", lwd = 2,lty=2)
#slight deviation from prior

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,400), ylim = c(0.706, 0.715), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(max.dist.mea)-n.avg.misha.25.dist.rmv)/14.7,n.avg.misha.25.sr.rmv, pch= 18, col="#00b3ffff")
#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.fdnb1pr.Rin.89<- MCMC.CI.bound(post.misha.fdnb1pr$BUGSoutput$sims.list$Rin, 0.89)
lines(1:400,post.misha.fdnb1pr.Rin.89[[1]],lwd = 2, col = "black")
lines(1:400,post.misha.fdnb1pr.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:400,post.misha.fdnb1pr.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

#estimated input series
MCMC.ts.Rin.m.inv1phr.param.89<- MCMC.CI.bound(post.misha.inv1phr.param$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:400,MCMC.ts.Rin.m.inv1phr.param.89[[1]],lwd = 2, col = "magenta")
lines(1:400,MCMC.ts.Rin.m.inv1phr.param.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:400,MCMC.ts.Rin.m.inv1phr.param.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Model input","Reconstructed input"),lwd = c(2, 2), col=c("black","magenta"))