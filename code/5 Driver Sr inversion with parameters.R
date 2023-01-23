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

############################# inversion of Misha using params#######
R.sd.mea <- n.sd.misha.50.sr.rmv
dist.mea <- n.avg.misha.50.dist.rmv
R.mea <- n.avg.misha.50.sr.rmv
n.mea = length(n.avg.misha.50.sr.rmv)

s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30

#posterior samples of parameters from Misha calibration

a.post <- post.misha.pc2p3$BUGSoutput$sims.list$a[,1]
b.post <- post.misha.pc2p3$BUGSoutput$sims.list$b[,1]
c.post <- post.misha.pc2p3$BUGSoutput$sims.list$c[,1]
post.leng <- length(a.post)

parameters <- c("Ivo.rate", "dist","Rin.m.pre",
                "R1.m","Rin.m", "R2.m","a.m","b.m", "c.m")

dat = list( s.intv = s.intv, max.dist.mea = max.dist.mea, post.leng=post.leng, 
            a.post=a.post, b.post=b.post, c.post=c.post,
            Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 750, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e3
n.burnin = 2e3
n.thin = 1

#Run it
post.misha.inv2p3.param = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param.R", 
                                                     parameters.to.save = parameters, 
                                                     data = dat, n.chains=5, n.iter = n.iter, 
                                                     n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 51 hours

save(post.misha.inv2p3.param, file = "out/post.misha.inv2p3.param.RData")

post.misha.inv2p3.param$BUGSoutput$summary

#check posterior density of parameter a:
plot(density(post.misha.pc2p3$BUGSoutput$sims.list$a),xlim=c(0,0.04),ylim=c(0,120))
lines(density(post.misha.inv2p3.param$BUGSoutput$sims.list$a),col="red")

plot(0,0, xlim = c(1,700), ylim = c(0.704, 0.714), xlab = "days", ylab ="Sr 87/86",
     main="Fidelity test: model input series vs. estimated input series")
#converting misha distance to days using rate Ivo.rate

#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.pc2p3.Rin.89<- MCMC.CI.bound(post.misha.pc2p3$BUGSoutput$sims.list$Rin, 0.89)
lines(1:750,post.misha.pc2p3.Rin.89[[1]],lwd = 2, col = "black")
lines(1:750,post.misha.pc2p3.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:750,post.misha.pc2p3.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

points((max(n.avg.misha.50.dist) + 30-n.avg.misha.50.dist.rmv)/14.7,n.avg.misha.50.sr.rmv, pch= 18, col="#00b3ffff")

#estimated input series
#estimated input series
MCMC.misha.inv2p3.param.Rin.m.89 <- MCMC.CI.bound(post.misha.inv2p3.param$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:750,MCMC.misha.inv2p3.param.Rin.m.89[[1]],lwd = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2p3.param.Rin.m.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2p3.param.Rin.m.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(400, 0.708, c("Model input","Estimated input"),lwd = c(2, 2), col=c("black","magenta"))

############################# inversion of Misha using different error per step params#######
R.sd.mea <- n.sd.misha.50.sr.rmv
dist.mea <- n.avg.misha.50.dist.rmv
R.mea <- n.avg.misha.50.sr.rmv
n.mea = length(n.avg.misha.50.sr.rmv)

s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30

#posterior samples of parameters from Misha calibration

a.post <- post.misha.pc2p3$BUGSoutput$sims.list$a[,1]
b.post <- post.misha.pc2p3$BUGSoutput$sims.list$b[,1]
c.post <- post.misha.pc2p3$BUGSoutput$sims.list$c[,1]
post.leng <- length(a.post)

parameters <- c("Ivo.rate", "dist","Rin.m.pre",
                "R1.m","Rin.m", "R2.m","a.m","b.m", "c.m")

dat = list( s.intv = s.intv, max.dist.mea = max.dist.mea, post.leng=post.leng, 
            a.post=a.post, b.post=b.post, c.post=c.post,
            Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 750, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e3
n.burnin = 2e3
n.thin = 1

#Run it
post.misha.inv2p3.param.s = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param s.R", 
                                                     parameters.to.save = parameters, 
                                                     data = dat, n.chains=5, n.iter = n.iter, 
                                                     n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 51 hours

save(post.misha.inv2p3.param.s, file = "out/post.misha.inv2p3.param.s.RData")

post.misha.inv2p3.param.s$BUGSoutput$summary

#check posterior density of parameter a:
plot(density(post.misha.pc2p3$BUGSoutput$sims.list$a),xlim=c(0,0.04),ylim=c(0,120))
lines(density(post.misha.inv2p3.param.s$BUGSoutput$sims.list$a),col="red")

plot(0,0, xlim = c(1,700), ylim = c(0.704, 0.714), xlab = "days", ylab ="Sr 87/86",
     main="Fidelity test: model input series vs. estimated input series")
#converting misha distance to days using rate Ivo.rate

#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.pc2p3.Rin.89<- MCMC.CI.bound(post.misha.pc2p3$BUGSoutput$sims.list$Rin, 0.89)
lines(1:750,post.misha.pc2p3.Rin.89[[1]],lwd = 2, col = "black")
lines(1:750,post.misha.pc2p3.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:750,post.misha.pc2p3.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

points((max(n.avg.misha.50.dist) + 30-n.avg.misha.50.dist.rmv)/14.7,n.avg.misha.50.sr.rmv, pch= 18, col="#00b3ffff")

#estimated input series
#estimated input series
MCMC.misha.inv2p3.param.s.Rin.m.89 <- MCMC.CI.bound(post.misha.inv2p3.param.s$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:750,MCMC.misha.inv2p3.param.s.Rin.m.89[[1]],lwd = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2p3.param.s.Rin.m.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2p3.param.s.Rin.m.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(400, 0.708, c("Model input","Estimated input"),lwd = c(2, 2), col=c("black","magenta"))