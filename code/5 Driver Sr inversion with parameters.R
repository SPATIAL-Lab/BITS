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

############################# inversion of Misha#######
R.sd.mea <- n.sd.misha.50.sr.rmv
dist.mea <- n.avg.misha.50.dist.rmv
R.mea <- n.avg.misha.50.sr.rmv
n.mea = length(n.avg.misha.50.sr.rmv)

s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30

#posterior samples of parameters from Misha calibration

a.post <- post.misha.pc2p$BUGSoutput$sims.list$a[,1]
b.post <- post.misha.pc2p$BUGSoutput$sims.list$b[,1]
c.post <- post.misha.pc2p$BUGSoutput$sims.list$c[,1]
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
n.thin = floor(n.iter-n.burnin)/400

#Run it
post.misha.inv2pr.param = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param.R", 
                                                      parameters.to.save = parameters, 
                                                      data = dat, n.chains=5, n.iter = n.iter, 
                                                      n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 20 hours

save(post.misha.inv2pr.param, file = "out/post.misha.inv2pr.param.RData")

post.misha.inv2pr.param$BUGSoutput$summary


#plotting reconstructed Rin.m history
par(mfrow=c(1,1))
plot(0,0, xlim = c(1,750), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(n.avg.misha.50.dist) + 30-n.avg.misha.50.dist.rmv)/14.7,n.avg.misha.50.sr.rmv, pch= 18, col="#00b3ffff")
#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.pc2p.Rin.89<- MCMC.CI.bound(post.misha.pc2p$BUGSoutput$sims.list$Rin, 0.89)
lines(1:750,post.misha.pc2p.Rin.89[[1]],lwd = 2, col = "black")
lines(1:750,post.misha.pc2p.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:750,post.misha.pc2p.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

#estimated input series
MCMC.misha.inv2pr.param.Rin.m.89 <- MCMC.CI.bound(post.misha.inv2pr.param$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:750,MCMC.misha.inv2pr.param.Rin.m.89[[1]],lwd = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2pr.param.Rin.m.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2pr.param.Rin.m.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Model input","Reconstructed input"),lwd = c(2, 2), col=c("black","magenta"))

#####inversion of Misha full record 50 pt###########
R.sd.mea <- n.sd.misha.50.sr
dist.mea <- n.avg.misha.50.dist
R.mea <- n.avg.misha.50.sr
n.mea = length(n.avg.misha.50.sr)

s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30

#posterior samples of parameters from Misha calibration

a.post <- post.misha.pc2p$BUGSoutput$sims.list$a[,1]
b.post <- post.misha.pc2p$BUGSoutput$sims.list$b[,1]
c.post <- post.misha.pc2p$BUGSoutput$sims.list$c[,1]
post.leng <- length(a.post)

parameters <- c("Ivo.rate", "dist",
                "R1.m","Rin.m", "R2.m","a.m","b.m", "c.m")

dat = list( s.intv = s.intv, max.dist.mea = max.dist.mea, post.leng=post.leng, 
            a.post=a.post, b.post=b.post, c.post=c.post,
            Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 750, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/400

#Run it
post.misha.inv2pf.param = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param.R", 
                                                     parameters.to.save = parameters, 
                                                     data = dat, n.chains=5, n.iter = n.iter, 
                                                     n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 4.3 hours

save(post.misha.inv2pf.param, file = "out/post.misha.inv2pf.param.RData")

post.misha.inv2pf.param$BUGSoutput$summary


#plotting reconstructed Rin.m history
par(mfrow=c(1,1))
plot(0,0, xlim = c(1,750), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(n.avg.misha.50.dist) + 30-n.avg.misha.50.dist)/14.7,n.avg.misha.50.sr, pch= 18, col="#00b3ffff")
#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.pc2p.Rin.89<- MCMC.CI.bound(post.misha.pc2p$BUGSoutput$sims.list$Rin, 0.89)
lines(1:750,post.misha.pc2p.Rin.89[[1]],lwd = 2, col = "black")
lines(1:750,post.misha.pc2p.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:750,post.misha.pc2p.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

#estimated input series
MCMC.misha.inv2pf.param.Rin.m.89 <- MCMC.CI.bound(post.misha.inv2pf.param$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:750,MCMC.misha.inv2pf.param.Rin.m.89[[1]],lwd = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2pf.param.Rin.m.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2pf.param.Rin.m.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Model input","Reconstructed input"),lwd = c(2, 2), col=c("black","magenta"))

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
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 700, n.mea = n.mea)

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


#plotting reconstructed Rin.m history
par(mfrow=c(1,1))
plot(0,0, xlim = c(1,700), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(n.avg.misha.50.dist) + 30-n.avg.misha.50.dist.remv)/14.7,n.avg.misha.50.sr.remv, pch= 18, col="#00b3ffff")
#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.pc2p.erm.Rin.89<- MCMC.CI.bound(post.misha.pc2p.erm$BUGSoutput$sims.list$Rin, 0.89)
lines(1:700,post.misha.pc2p.erm.Rin.89[[1]],lwd = 2, col = "black")
lines(1:700,post.misha.pc2p.erm.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:700,post.misha.pc2p.erm.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

#estimated input series
MCMC.misha.inv2perm.param.Rin.m.89 <- MCMC.CI.bound(post.misha.inv2perm.param$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:700,MCMC.misha.inv2perm.param.Rin.m.89[[1]],lwd = 2, col = "magenta")
lines(1:700,MCMC.misha.inv2perm.param.Rin.m.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:700,MCMC.misha.inv2perm.param.Rin.m.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Model input","Reconstructed input"),lwd = c(2, 2), col=c("black","magenta"))

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
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 700, n.mea = n.mea)

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


#plotting reconstructed Rin.m history
par(mfrow=c(1,1))
plot(0,0, xlim = c(1,700), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(n.avg.misha.50.dist) + 30-n.avg.misha.50.dist.remv)/14.7,n.avg.misha.50.sr.remv, pch= 18, col="#00b3ffff")
#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.pc2p.erm.Rin.89<- MCMC.CI.bound(post.misha.pc2p.erm$BUGSoutput$sims.list$Rin, 0.89)
lines(1:700,post.misha.pc2p.erm.Rin.89[[1]],lwd = 2, col = "black")
lines(1:700,post.misha.pc2p.erm.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:700,post.misha.pc2p.erm.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

#estimated input series
MCMC.misha.inv2perm.tsrwca.Rin.m.89 <- MCMC.CI.bound(post.misha.inv2perm.tsrwca$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:700,MCMC.misha.inv2perm.tsrwca.Rin.m.89[[1]],lwd = 2, col = "magenta")
lines(1:700,MCMC.misha.inv2perm.tsrwca.Rin.m.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:700,MCMC.misha.inv2perm.tsrwca.Rin.m.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Model input","Reconstructed input"),lwd = c(2, 2), col=c("black","magenta"))
