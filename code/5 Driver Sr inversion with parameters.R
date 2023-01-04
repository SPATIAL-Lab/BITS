library(coda)
library(rjags)
library(R2jags)
library(mcmcplots)
library(bayestestR)
library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)

source("code/1 Helper functions.R")

plot.col<-viridis(7)

setwd("C:/Users/ydmag/Google Drive/U of U/Elephant movement/Sr-in-ivory")

############################# inversion of Misha using erm params#######
R.sd.mea <- n.sd.misha.50.sr.remv
dist.mea <- n.avg.misha.50.dist.remv
R.mea <- n.avg.misha.50.sr.remv
n.mea = length(n.avg.misha.50.sr.remv)

s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30

#posterior samples of parameters from Misha calibration

a.post <- post.misha.pc2p.erm$BUGSoutput$sims.list$a[,1]
b.post <- post.misha.pc2p.erm$BUGSoutput$sims.list$b[,1]
c.post <- post.misha.pc2p.erm$BUGSoutput$sims.list$c[,1]
post.leng <- length(a.post)

parameters <- c("Ivo.rate", "dist","Rin.m.pre",
                "R1.m","Rin.m", "R2.m","a.m","b.m", "c.m")

dat = list( s.intv = s.intv, max.dist.mea = max.dist.mea, post.leng=post.leng, 
            a.post=a.post, b.post=b.post, c.post=c.post,
            Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 680, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e3
n.burnin = 2e3
n.thin = 1#floor(n.iter-n.burnin)/400

#Run it
post.misha.inv2perm.param = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param.R", 
                                                     parameters.to.save = parameters, 
                                                     data = dat, n.chains=5, n.iter = n.iter, 
                                                     n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 11 hours

save(post.misha.inv2perm.param, file = "out/post.misha.inv2perm.param.RData")

post.misha.inv2perm.param$BUGSoutput$summary

#check posterior density of parameter a:
plot(density(post.misha.pc2p.erm$BUGSoutput$sims.list$a),xlim=c(0.01,0.05),ylim=c(0,120))
lines(density(post.misha.inv2perm.param$BUGSoutput$sims.list$a),col="red")

#######sensitivity test########
#####different time series error with auto correlation#####
############################# inversion of Misha using erm params#######
R.sd.mea <- n.sd.misha.50.sr.remv
dist.mea <- n.avg.misha.50.dist.remv
R.mea <- n.avg.misha.50.sr.remv
n.mea = length(n.avg.misha.50.sr.remv)

s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30

#posterior samples of parameters from Misha calibration

a.post <- post.misha.pc2p.erm$BUGSoutput$sims.list$a[,1]
b.post <- post.misha.pc2p.erm$BUGSoutput$sims.list$b[,1]
c.post <- post.misha.pc2p.erm$BUGSoutput$sims.list$c[,1]
post.leng <- length(a.post)

parameters <- c("Ivo.rate", "dist","Rin.m.pre","Rin.m.cps.ac",
                "R1.m","Rin.m", "R2.m","a.m","b.m", "c.m")

dat = list( s.intv = s.intv, max.dist.mea = max.dist.mea, post.leng=post.leng, 
            a.post=a.post, b.post=b.post, c.post=c.post,
            Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 680, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e3
n.burnin = 2e3
n.thin = 1#floor(n.iter-n.burnin)/400

#Run it
post.misha.inv2perm.tsrwca = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param tsrwca.R", 
                                                       parameters.to.save = parameters, 
                                                       data = dat, n.chains=5, n.iter = n.iter, 
                                                       n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 26 hours

save(post.misha.inv2perm.tsrwca, file = "out/post.misha.inv2perm.tsrwca.RData")

post.misha.inv2perm.tsrwca$BUGSoutput$summary

plot(density(post.misha.pc2p.erm$BUGSoutput$sims.list$a))
lines(density(post.misha.inv2perm.tsrwca$BUGSoutput$sims.list$a),col="red") #decrease error to 3e-8

plot(density(post.misha.inv2perm.tsrwca$BUGSoutput$sims.list$Rin.m.cps.ac)) #almost no autocorrelation


#plotting reconstructed Rin.m history
par(mfrow=c(1,1))
plot(0,0, xlim = c(1,650), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(n.avg.misha.50.dist) + 30-n.avg.misha.50.dist.remv)/14.7,n.avg.misha.50.sr.remv, pch= 18, col="#00b3ffff")
#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.pc2p.erm.Rin.89<- MCMC.CI.bound(post.misha.pc2p.erm$BUGSoutput$sims.list$Rin, 0.89)
lines(1:700,post.misha.pc2p.erm.Rin.89[[1]],lwd = 2, col = "black")
lines(1:700,post.misha.pc2p.erm.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:700,post.misha.pc2p.erm.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

#estimated input series
MCMC.misha.inv2perm.tsrwca.Rin.m.89 <- MCMC.CI.bound(post.misha.inv2perm.tsrwca$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:680,MCMC.misha.inv2perm.tsrwca.Rin.m.89[[1]],lwd = 2, col = "magenta")
lines(1:680,MCMC.misha.inv2perm.tsrwca.Rin.m.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:680,MCMC.misha.inv2perm.tsrwca.Rin.m.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Model input","Reconstructed input"),lwd = c(2, 2), col=c("black","magenta"))
