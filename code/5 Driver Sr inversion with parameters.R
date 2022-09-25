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

########inversion model with calibration for either the one/two pool model#######
#inversion based on micromill results
micromill <- read.csv("data/Misha micromill.csv")

R.mic <- rev(micromill$X87Sr.86Sr)
R.sd.mic <- rev(micromill$Sd)
dist.mic <- rev(micromill$dist)
n.mic <- length(dist.mic)

Ivo.rate.mic <- 20

#sample interval
s.intv <- 400

#parameters to save
parameters <- c("Body.mass", "Body.mass.m","Rin.m.cps",
                "Rs.m", "Rb.m","Rin.m","dist", "a.m", "b.m", "c.m","Rin.m.cps.ac")

##Data to pass to the model
#omitting the turnover model here simplifies the model run time
#the inversion takes the measured value of potentially a different ivory series
dat = list(s.intv = s.intv, 
           params.mu = turnover.params.mu, params.vcov = turnover.params.vcov,
           R.mea = R.mic, dist.mea = dist.mic, R.sd.mea = R.sd.mic, t = 450, n.mea = n.mic)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 4e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/600

#Run it
post.misha.inv.p = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS with parameters.R", 
                                                   parameters.to.save = parameters, 
                                                   data = dat, n.chains=5, n.iter = n.iter, 
                                                   n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 8 hours

save(post.misha.inv.p, file = "out/post.misha.inv.p.RData")

post.misha.inv.p$BUGSoutput$summary

load("out/post.misha.inv.p.RData")

plot(density(post.misha.inv.p$BUGSoutput$sims.list$Rin.m.cps.ac))
plot(density(post.misha.inv.p$BUGSoutput$sims.list$Rin.m.cps))
plot(density(post.misha.inv.p$BUGSoutput$sims.list$Body.mass.m))

#plotting reconstructed Rin.m history
plot(0,0, xlim = c(1,450), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines((max(dist.mic+200) - dist.mic)/Ivo.rate.mic + 1, R.mic, lwd = 2, col = "blue") #approximate results from micromill

MCMC.inv.p.Rin.m.50 <- MCMC.CI.bound(post.misha.inv.p$BUGSoutput$sims.list$Rin.m, 0.5)
lines(1:450,MCMC.inv.p.Rin.m.50[[1]],lwd = 2, col = "firebrick4")
lines(1:450,MCMC.inv.p.Rin.m.50[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:450,MCMC.inv.p.Rin.m.50[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))

####brownian motion per step model####
#inversion of original misha series##
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
parameters <- c("Body.mass", "Body.mass.m","Rin.m.cps","Ivo.rate",
                "Rs.m", "Rb.m","Rin.m","dist", "a.m", "b.m", "c.m","Rin.m.cps.ac")

##Data to pass to the model
#omitting the turnover model here simplifies the model run time
#the inversion takes the measured value of potentially a different ivory series
dat = list(s.intv = cal.intv, 
           params.mu = turnover.params.mu, params.vcov = turnover.params.vcov,
           Ivo.rate.mean = Ivo.rate.cal.mean, Ivo.rate.sd = Ivo.rate.cal.sd,
           max.dist.mea = max.dist.cal,
           R.mea = R.cal, dist.mea = dist.cal, R.sd.mea = R.sd.cal, t = 750, n.mea = n.cal)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 4e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/600

#Run it
post.misha.inv.bm.l = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS BM with parameters.R", 
                                              parameters.to.save = parameters, 
                                              data = dat, n.chains=5, n.iter = n.iter, 
                                              n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 8 hours

save(post.misha.inv.bm.l, file = "out/post.misha.inv.bm.l.RData")

post.misha.inv.bm.l$BUGSoutput$summary

load("out/post.misha.inv.bm.l.RData")

#plotting modeled serum values mapped onto ivory and checking the fit of the data


plot(density(post.misha.inv.bm.l$BUGSoutput$sims.list$Rin.m.cps))
plot(density(post.misha.inv.bm.l$BUGSoutput$sims.list$Body.mass.m))

#Rs.m
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.inv.bm.l$BUGSoutput$sims.list$Rs.m,
               post.misha.inv.bm.l$BUGSoutput$sims.list$dist)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
post.misha.inversion.woi.Rs.cal.89<- MCMC.CI.bound(post.misha.inversion.woi$BUGSoutput$sims.list$Rs.m, 0.89)

#plotting reconstructed Rin.m history
plot(0,0, xlim = c(1,t), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines((max(misha$dist)-misha$dist)/Ivo.rate.cal.mean,misha$mean,lwd = 2, col = "blue") #approximate results from micromill

MCMC.inv.bm.l.Rin.m.89 <- MCMC.CI.bound(post.misha.inv.bm.l$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:t,MCMC.inv.bm.l.Rin.m.89[[1]],lwd = 2, col = "firebrick4")
lines(1:t,MCMC.inv.bm.l.Rin.m.89[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:t,MCMC.inv.bm.l.Rin.m.89[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))


plot(density(post.misha.inv.bm.l$BUGSoutput$sims.list$Ivo.rate))

#######brownian motion per step model with cauchy error####
#inversion of original misha series##
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
parameters <- c("Body.mass", "Body.mass.m","Rin.m.cps","Ivo.rate",
                "Rs.m", "Rb.m","Rin.m","dist", "a.m", "b.m", "c.m","Rin.m.cps.ac")

##Data to pass to the model
#omitting the turnover model here simplifies the model run time
#the inversion takes the measured value of potentially a different ivory series
dat = list(s.intv = cal.intv, 
           params.mu = turnover.params.mu, params.vcov = turnover.params.vcov,
           Ivo.rate.mean = Ivo.rate.cal.mean, Ivo.rate.sd = Ivo.rate.cal.sd,
           max.dist.mea = max.dist.cal,
           R.mea = R.cal, dist.mea = dist.cal, R.sd.mea = R.sd.cal, t = 750, n.mea = n.cal)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 4e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/600

#Run it
post.misha.inv.bmca.l = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS BM with parameters.R", 
                                                 parameters.to.save = parameters, 
                                                 data = dat, n.chains=5, n.iter = n.iter, 
                                                 n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 8 hours

save(post.misha.inv.bmca.l, file = "out/post.misha.inv.bm.l.RData")

post.misha.inv.bmca.l$BUGSoutput$summary

load("out/post.misha.inv.bm.l.RData")

#plotting modeled serum values mapped onto ivory and checking the fit of the data


plot(density(post.misha.inv.bmca.l$BUGSoutput$sims.list$Rin.m.cps))
plot(density(post.misha.inv.bmca.l$BUGSoutput$sims.list$Body.mass.m))

#Rs.m
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.inv.bmca.l$BUGSoutput$sims.list$Rs.m,
               post.misha.inv.bmca.l$BUGSoutput$sims.list$dist)
lines(misha$dist,misha$mean,lwd = 2, col = "red")

#plotting reconstructed Rin.m history
plot(0,0, xlim = c(1,t), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines((max(misha$dist)-misha$dist)/Ivo.rate.cal.mean,misha$mean,lwd = 2, col = "blue") #approximate results from micromill

MCMC.inv.bmca.l.Rin.m.89 <- MCMC.CI.bound(post.misha.inv.bmca.l$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:t,MCMC.inv.bmca.l.Rin.m.89[[1]],lwd = 2, col = "firebrick4")
lines(1:t,MCMC.inv.bmca.l.Rin.m.89[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:t,MCMC.inv.bmca.l.Rin.m.89[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))

par(mfrow=c(1,3))
plot(abc.prior.params$x,abc.prior.params$y[1,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "a", ylab= "density")
lines(density(log(post.misha.inv.bmca.l$BUGSoutput$sims.list$a.m)), col = "red", lwd = 2)

plot(abc.prior.params$x,abc.prior.params$y[2,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "b", ylab= "density")
lines(density(log(post.misha.inv.bmca.l$BUGSoutput$sims.list$b.m)), col = "red", lwd = 2)

plot(abc.prior.params$x,abc.prior.params$y[3,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "c", ylab= "density")
lines(density(log(post.misha.inv.bmca.l$BUGSoutput$sims.list$c.m)), col = "red", lwd = 2)

#########################################iversion with a different correlation structure####
#inversion of original misha series##
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
parameters <- c("Body.mass", "Body.mass.m","Rin.m.cps","Ivo.rate",
                "Rs.m", "Rb.m","Rin.m","dist", "a.m", "b.m", "c.m","Rin.m.cps.ac")

##Data to pass to the model
#omitting the turnover model here simplifies the model run time
#the inversion takes the measured value of potentially a different ivory series
dat = list(s.intv = cal.intv, 
           params.mu = turnover.params.mu, params.vcov = turnover.params.vcov,
           Ivo.rate.mean = Ivo.rate.cal.mean, Ivo.rate.sd = Ivo.rate.cal.sd,
           max.dist.mea = max.dist.cal,
           R.mea = R.cal, dist.mea = dist.cal, R.sd.mea = R.sd.cal, t = 750, n.mea = n.cal)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 4e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/600

#Run it
post.misha.inv.l = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS with parameters.R", 
                                                 parameters.to.save = parameters, 
                                                 data = dat, n.chains=5, n.iter = n.iter, 
                                                 n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 8 hours

save(post.misha.inv.l, file = "out/post.misha.inv.bm.l.RData")

post.misha.inv.l$BUGSoutput$summary

load("out/post.misha.inv.bm.l.RData")

#plotting modeled serum values mapped onto ivory and checking the fit of the data

plot(density(post.misha.inv.l$BUGSoutput$sims.list$Rin.m.cps.ac))
plot(density(post.misha.inv.l$BUGSoutput$sims.list$Rin.m.cps))
plot(density(post.misha.inv.l$BUGSoutput$sims.list$Body.mass.m))

#Rs.m
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.inv.l$BUGSoutput$sims.list$Rs.m,
               post.misha.inv.l$BUGSoutput$sims.list$dist)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
post.misha.inv.l.Rs.m.89<- MCMC.CI.bound(post.misha.inv.l$BUGSoutput$sims.list$Rs.m, 0.89)

#plotting reconstructed Rin.m history
plot(0,0, xlim = c(1,750), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines((max(misha$dist)-misha$dist)/Ivo.rate.cal.mean,misha$mean,lwd = 2, col = "blue") #approximate results from micromill

MCMC.inv.l.Rin.m.89 <- MCMC.CI.bound(post.misha.inv.l$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:750,MCMC.inv.l.Rin.m.89[[1]],lwd = 2, col = "firebrick4")
lines(1:750,MCMC.inv.l.Rin.m.89[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:750,MCMC.inv.l.Rin.m.89[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))
#brownian motion is actually okay...

####brownian motion per step model####
micromill <- read.csv("data/Misha micromill.csv")

R.mic <- rev(micromill$X87Sr.86Sr)
R.sd.mic <- rev(micromill$StdErr)
dist.mic <- rev(micromill$dist)
n.mic <- length(dist.mic)

Ivo.rate.mic <- 20

#sample interval
s.intv <- 400

max.dist.mea <- 8000

Ivo.rate.mean <- 19.9

Ivo.rate.sd <- 5.4

#parameters to save
parameters <- c("Body.mass", "Body.mass.m","Rin.m.cps",
                "Rs.m", "Rb.m","Rin.m","dist", "a.m", "b.m", "c.m","Rin.m.cps.ac")

##Data to pass to the model
#omitting the turnover model here simplifies the model run time
#the inversion takes the measured value of potentially a different ivory series
dat = list(s.intv = s.intv, 
           params.mu = turnover.params.mu, params.vcov = turnover.params.vcov,
           Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
           max.dist.mea = max.dist.mea,
           R.mea = R.mic, dist.mea = dist.mic, R.sd.mea = R.sd.mic, t = 450, n.mea = n.mic)



#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 4e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/600

#Run it
post.misha.inv.bm = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS BM with parameters.R", 
                                               parameters.to.save = parameters, 
                                               data = dat, n.chains=5, n.iter = n.iter, 
                                               n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 8 hours

save(post.misha.inv.bm, file = "out/post.misha.inv.bm.RData")

post.misha.inv.bm$BUGSoutput$summary

load("out/post.misha.inv.bm.RData")

#plotting modeled serum values mapped onto ivory and checking the fit of the data

plot(density(post.misha.inv.bm$BUGSoutput$sims.list$Rin.m.cps.ac))
plot(density(post.misha.inv.bm$BUGSoutput$sims.list$Rin.m.cps))
plot(density(post.misha.inv.bm$BUGSoutput$sims.list$Body.mass.m))

#plotting reconstructed Rin.m history
plot(0,0, xlim = c(1,450), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines((max(dist.mic+200) - dist.mic)/Ivo.rate.mic + 1, R.mic, lwd = 2, col = "blue") #approximate results from micromill

MCMC.inv.bm.Rin.m.89 <- MCMC.CI.bound(post.misha.inv.bm$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:450,MCMC.inv.bm.Rin.m.89[[1]],lwd = 2, col = "firebrick4")
lines(1:450,MCMC.inv.bm.Rin.m.89[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:450,MCMC.inv.bm.Rin.m.89[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))

####brownian motion with cauchy error per step model####
micromill <- read.csv("data/Misha micromill.csv")

R.mic <- rev(micromill$X87Sr.86Sr)
R.sd.mic <- rev(micromill$Sd)
dist.mic <- rev(micromill$dist)
n.mic <- length(dist.mic)

#sample interval
s.intv <- 400

max.dist.mea <- 8500 #add a bit of uncertainty before the evaluation

Ivo.rate.mean <- 19.9

Ivo.rate.sd <- 5.4

#parameters to save
parameters <- c("Body.mass","Rin.m.cps",
                "Rs.m", "Rb.m","Rin.m","dist", "a.m", "b.m", "c.m")

##Data to pass to the model
#omitting the turnover model here simplifies the model run time
#the inversion takes the measured value of potentially a different ivory series
dat = list(s.intv = s.intv, 
           params.mu = turnover.params.mu, params.vcov = turnover.params.vcov,
           Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
           max.dist.mea = max.dist.mea,
           R.mea = R.mic, dist.mea = dist.mic, R.sd.mea = R.sd.mic, t = 550, n.mea = n.mic)



#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 4e4
n.burnin = 4e3
n.thin = floor(n.iter-n.burnin)/500

#Run it
post.mic.inv.bmca = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS BM with parameters.R", 
                                                 parameters.to.save = parameters, 
                                                 data = dat, n.chains=4, n.iter = n.iter, 
                                                 n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 1.5 hours

save(post.mic.inv.bmca, file = "out/post.mic.inv.bm.RData")

post.mic.inv.bmca$BUGSoutput$summary

load("out/post.mic.inv.bmca.RData")

#plotting modeled serum values mapped onto ivory and checking the fit of the data

plot(density(post.mic.inv.bmca$BUGSoutput$sims.list$Rin.m.cps))
plot(density(post.mic.inv.bmca$BUGSoutput$sims.list$Body.mass.m))

#plotting reconstructed Rin.m history
par(mfrow=c(1,1))
plot(0,0, xlim = c(1,550), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines((max.dist.mea - dist.mic)/Ivo.rate.mic + 1, R.mic, lwd = 2, col = "blue") #approximate results from micromill

MCMC.inv.bmca.Rin.m.89 <- MCMC.CI.bound(post.mic.inv.bmca$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:550,MCMC.inv.bmca.Rin.m.89[[1]],lwd = 2, col = "firebrick4")
lines(1:550,MCMC.inv.bmca.Rin.m.89[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:550,MCMC.inv.bmca.Rin.m.89[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))

par(mfrow=c(1,1))
plot(0,0, xlim = c(1,550), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines((max.dist.mea - dist.mic)/Ivo.rate.mic + 1, R.mic, lwd = 2, col = "blue") #approximate results from micromill

MCMC.inv.bmca.Rs.m.89 <- MCMC.CI.bound(post.mic.inv.bmca$BUGSoutput$sims.list$Rs.m, 0.89)
lines(1:550,MCMC.inv.bmca.Rs.m.89[[1]],lwd = 2, col = "firebrick4")
lines(1:550,MCMC.inv.bmca.Rs.m.89[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:550,MCMC.inv.bmca.Rs.m.89[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed Rs"),lwd = c(2, 2), col=c("blue","firebrick4"))

par(mfrow=c(1,1))
plot(0,0, xlim = c(1,550), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines((max.dist.mea - dist.mic)/Ivo.rate.mic + 1, R.mic, lwd = 2, col = "blue") #approximate results from micromill

MCMC.inv.bmca.Rb.m.89 <- MCMC.CI.bound(post.mic.inv.bmca$BUGSoutput$sims.list$Rb.m, 0.89)
lines(1:550,MCMC.inv.bmca.Rb.m.89[[1]],lwd = 2, col = "firebrick4")
lines(1:550,MCMC.inv.bmca.Rb.m.89[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:550,MCMC.inv.bmca.Rb.m.89[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed Rb"),lwd = c(2, 2), col=c("blue","firebrick4"))

#check prior vs posterior params a, b, and c
abc.prior.params <- pri.multi.norm.den(-10,0,turnover.params.mu,turnover.params.vcov)
par(mfrow=c(1,3))

plot(abc.prior.params$x,abc.prior.params$y[1,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "a", ylab= "density")
lines(density(log(post.mic.inv.bmca$BUGSoutput$sims.list$a.m)), col = "red", lwd = 2)

plot(abc.prior.params$x,abc.prior.params$y[2,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "b", ylab= "density")
lines(density(log(post.mic.inv.bmca$BUGSoutput$sims.list$b.m)), col = "red", lwd = 2)

plot(abc.prior.params$x,abc.prior.params$y[3,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "c", ylab= "density")
lines(density(log(post.mic.inv.bmca$BUGSoutput$sims.list$c.m)), col = "red", lwd = 2)
