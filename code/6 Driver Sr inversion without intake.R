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
R.sd.mic <- rev(micromill$StdErr)
dist.mic <- rev(micromill$dist)
n.mic <- length(dist.mic)

#assign input values before and after the switch, but allows some variation
R0 <- 0.70706

#Re is the mean ratio of end value  
Re <- 0.7112

max.dist.cal <- 19200

max.dist.mea <- 8000

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
           R.mea = R.mic, dist.mea = dist.mic, R.sd.mea = R.sd.mic, t = 450, n.mea = n.mic)

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
proc.time() - t1 #~ 10 hours

save(post.misha.inversion.woi, file = "out/post.misha.inversion.woi.RData")

post.misha.inversion.woi$BUGSoutput$summary

load("out/post.misha.inversion.woi.RData")

plot(density(post.misha.inversion.woi$BUGSoutput$sims.list$a))

#plot calibration curve
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.inversion.woi$BUGSoutput$sims.list$Rs.cal,
               post.misha.inversion.woi$BUGSoutput$sims.list$dist.cal.m)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
post.misha.inversion.woi.Rs.cal.89<- MCMC.CI.bound(post.misha.inversion.woi$BUGSoutput$sims.list$Rs.cal, 0.89)

med.dist.inversion.woi<- MCMC.dist.median(post.misha.inversion.woi$BUGSoutput$sims.list$dist.cal.m)
lines(med.dist.inversion.woi, post.misha.inversion.woi.Rs.cal.89[[1]], lwd = 1, col = "cyan")
lines(med.dist.inversion.woi, post.misha.inversion.woi.Rs.cal.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(med.dist.inversion.woi, post.misha.inversion.woi.Rs.cal.89[[3]], lwd = 1, lty = 2, col = "cyan")

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,t), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to days using rate Ivo.rate.mic
lines((max(dist.mic+200) - dist.mic)/Ivo.rate.mic + 1, R.mic, lwd = 2, col = "blue")#results from micromill

MCMC.ts.Rin.m.89.woi <- MCMC.ts(post.misha.inversion.woi$BUGSoutput$sims.list$Rin.m)
lines(1:t,MCMC.ts.Rin.m.89.woi[[1]],lwd = 2, col = "firebrick4")
lines(1:t,MCMC.ts.Rin.m.89.woi[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:t,MCMC.ts.Rin.m.89.woi[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))


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
proc.time() - t1 #~ 10 hours

save(post.misha.inversion.woi, file = "out/post.misha.inversion.woi.RData")

post.misha.inversion.woi$BUGSoutput$summary

load("out/post.misha.inversion.woi.RData")

plot(density(post.misha.inversion.woi$BUGSoutput$sims.list$a))

#plot calibration curve
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
plot(0,0, xlim = c(1,t), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to days using rate Ivo.rate.mic
lines((max(misha$dist)-misha$dist)/Ivo.rate.cal.mean,misha$mean,lwd = 2, col = "blue")

MCMC.ts.Rin.m.89.woi <- MCMC.ts(post.misha.inversion.woi$BUGSoutput$sims.list$Rin.m)
lines(1:t,MCMC.ts.Rin.m.89.woi[[1]],lwd = 2, col = "firebrick4")
lines(1:t,MCMC.ts.Rin.m.89.woi[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:t,MCMC.ts.Rin.m.89.woi[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))
