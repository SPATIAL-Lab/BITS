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

#####################################################################################
################################      Section 1      ################################
#####################################################################################

###1 sensitivity to Re value################
#2 pool turnover with 50 pt average (excursion rmv) sensitivity test Re =0.7118##
###prep data###
R.sd.mea <- n.sd.misha.50.sr.rmv
dist.mea <- n.avg.misha.50.dist.rmv
R.mea <- n.avg.misha.50.sr.rmv
n.mea = length(n.avg.misha.50.sr.rmv)

R0 <- 0.706

#Re is the mean ratio of end value  
Re <- 0.7118
s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30
#parameters to save
parameters <- c("a", "b","c", "Ivo.rate", "R1.m","R2.m","Rin","dist.index",
                "Sr.pre", "Re.mean", "R0.mean", "switch","dist",
                "flux.ratio", "pool.ratio","Body.mass")
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
post.misha.pc2p8 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS wo intake model3.R", 
                                              parameters.to.save = parameters, 
                                              data = dat, n.chains=5, n.iter = n.iter, 
                                              n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 8 hours

save(post.misha.pc2p8, file = "out/post.misha.pc2p8.RData")
load("out/post.misha.pc2p8.RData")

traplot(post.misha.pc2p8,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.pc2p8,parms = c("a", "b","c"))

#convergence params are not so good!
subset(post.misha.pc2p8$BUGSoutput$summary,
       rownames(post.misha.pc2p4$BUGSoutput$summary)=="a")
subset(post.misha.pc2p8$BUGSoutput$summary,
       rownames(post.misha.pc2p4$BUGSoutput$summary)=="b")
subset(post.misha.pc2p8$BUGSoutput$summary,
       rownames(post.misha.pc2p4$BUGSoutput$summary)=="c")

#plotting modeled serum values mapped onto ivory and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.pc2p8$BUGSoutput$sims.list$R1.m,
               post.misha.pc2p8$BUGSoutput$sims.list$dist)
lines(n.avg.misha.50.dist.rmv,n.avg.misha.50.sr.rmv,lwd = 2, col = "#00b4ffff")


#compare calibration results with different priors#
plot(density(post.misha.pc2p3$BUGSoutput$sims.list$a, from = 0), xlim = c(0,0.06),ylim= c(0,180),
     lwd = 2, col = plot.col.6[3],main="Two-pool model", xlab="parameter estimate")
lines(density(post.misha.pc2p3$BUGSoutput$sims.list$b, from = 0),lwd = 2, col = plot.col.6[4])
lines(density(post.misha.pc2p3$BUGSoutput$sims.list$c, from = 0),lwd = 2, col = plot.col.6[5])
legend(0.02, 160, c("a, Two-pool","b, Two-pool", "c, Two-pool", "a, One-pool"),
       lwd = rep(2, 4), col=c(plot.col.6[3:5],plot.col.6[2]))

plot(density(post.misha.pc2p8$BUGSoutput$sims.list$a, from = 0), xlim = c(0,0.06),ylim= c(0,180),
     lwd = 2, col = plot.col.6[3],main="Two-pool model", xlab="parameter estimate")
lines(density(post.misha.pc2p8$BUGSoutput$sims.list$b, from = 0),lwd = 2, col = plot.col.6[4])
lines(density(post.misha.pc2p8$BUGSoutput$sims.list$c, from = 0),lwd = 2, col = plot.col.6[5])
legend(0.02, 160, c("a, Two-pool","b, Two-pool", "c, Two-pool", "a, One-pool"),
       lwd = rep(2, 4), col=c(plot.col.6[3:5],plot.col.6[2]))

#####sensitivity to the day of switch#######
#Try switch that is beyond 83 +- 5?
###prep data###
R.sd.mea <- n.sd.misha.50.sr.remv
dist.mea <- n.avg.misha.50.dist.remv
R.mea <- n.avg.misha.50.sr.remv
n.mea = length(n.avg.misha.50.sr.remv)

R0 <- 0.706

#Re is the mean ratio of end value  
Re <- 0.711
s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30

#Estimated date of switch
#by default, 83 is used, here we try 86, 80
#these values are expected to affect
switch <- 86 
#parameters to save
parameters <- c("a", "b","c", "Ivo.rate", "R1.m","R2.m","Rin","dist.index",
                "Sr.pre", "Re.mean", "R0.mean", "switch","dist",
                "flux.ratio", "pool.ratio","Body.mass")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 700, n.mea = n.mea, 
           R0 = R0, Re = Re, s.intv = s.intv, Ivo.rate.mean = Ivo.rate.mean,
           Ivo.rate.sd = Ivo.rate.sd, max.dist.mea = max.dist.mea, switch = switch)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e3
n.burnin = 4e2
n.thin = 1

#Run it
post.misha.woisw86 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS woi sw.R", 
                                              parameters.to.save = parameters, 
                                              data = dat, n.chains=5, n.iter = n.iter, 
                                              n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 8 hours

save(post.misha.woisw86, file = "out/post.misha.woisw86.RData")
load("out/post.misha.woisw86.RData")

traplot(post.misha.woisw86,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.woisw86,parms = c("a", "b","c"))

#convergence params are not so good!
subset(post.misha.woisw86$BUGSoutput$summary,
       rownames(post.misha.woisw86$BUGSoutput$summary)=="a")
subset(post.misha.woisw86$BUGSoutput$summary,
       rownames(post.misha.woisw86$BUGSoutput$summary)=="b")
subset(post.misha.woisw86$BUGSoutput$summary,
       rownames(post.misha.woisw86$BUGSoutput$summary)=="c")

plot(density(post.misha.woisw86$BUGSoutput$sims.list$Ivo.rate))
abline(v=Ivo.rate.mean)

#plotting modeled serum values mapped onto ivory and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.woisw86$BUGSoutput$sims.list$R1.m,
               post.misha.woisw86$BUGSoutput$sims.list$dist)
lines(n.avg.misha.50.dist.rmv,n.avg.misha.50.sr.rmv,lwd = 2, col = "#00b4ffff")

##############################################comparing parameter a #
plot(density(post.misha.pc2p3$BUGSoutput$sims.list$a, from = 0), xlim = c(0,0.06),ylim= c(0,180),
     lwd = 2, col = plot.col.6[3],main="Two-pool model", xlab="parameter estimate")
legend(0.02, 160, c("a, Two-pool","b, Two-pool", "c, Two-pool", "a, One-pool"),
       lwd = rep(2, 4), col=c(plot.col.6[3:5],plot.col.6[2]))

plot(density(post.misha.woisw86$BUGSoutput$sims.list$a, from = 0), xlim = c(0,0.06),ylim= c(0,180),
     lwd = 2, col = plot.col.6[3],main="Two-pool model", xlab="parameter estimate")
legend(0.02, 160, c("a, Two-pool","b, Two-pool", "c, Two-pool", "a, One-pool"),
       lwd = rep(2, 4), col=c(plot.col.6[3:5],plot.col.6[2]))

#####################################################################################
################################      Section 2      ################################
#####################################################################################

#sensitivity to parameter constraints
#-1 set no constrain on any parameter
#######try version 5, no constraint######
R.sd.mea <- n.sd.misha.50.sr.rmv
dist.mea <- n.avg.misha.50.dist.rmv
R.mea <- n.avg.misha.50.sr.rmv
n.mea = length(n.avg.misha.50.sr.rmv)

R0 <- 0.706

#Re is the mean ratio of end value  
Re <- 0.711
s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30
#parameters to save
parameters <- c("a", "b","c", "Ivo.rate", "R1.m","R2.m","Rin","dist.index",
                "Sr.pre", "Re.mean", "R0.mean", "switch","dist",
                "flux.ratio", "pool.ratio")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 750, n.mea = n.mea, 
           R0 = R0, Re = Re, s.intv = s.intv, Ivo.rate.mean = Ivo.rate.mean,
           Ivo.rate.sd = Ivo.rate.sd, max.dist.mea = max.dist.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e3
n.burnin = 2e3
n.thin = 1

#Run it
post.misha.pc2p5 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS wo intake model5.R", 
                                              parameters.to.save = parameters, 
                                              data = dat, n.chains=5, n.iter = n.iter, 
                                              n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 26 hours

save(post.misha.pc2p5, file = "out/post.misha.pc2p5.RData")
load("out/post.misha.pc2p5.RData")

#Convergence is better with param constrain a > b
traplot(post.misha.pc2p5,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.pc2p5,parms = c("a", "b","c"))

plot(density(post.misha.pc2p5$BUGSoutput$sims.list$Ivo.rate))
abline(v=14.7)

#convergence params are ok, but pool constraints on parameter c
subset(post.misha.pc2p5$BUGSoutput$summary,
       rownames(post.misha.pc2p$BUGSoutput$summary)=="a")
subset(post.misha.pc2p5$BUGSoutput$summary,
       rownames(post.misha.pc2p$BUGSoutput$summary)=="b")
subset(post.misha.pc2p5$BUGSoutput$summary,
       rownames(post.misha.pc2p$BUGSoutput$summary)=="c")

MAP.a <- map_estimate(post.misha.pc2p5$BUGSoutput$sims.list$a)
MAP.a[1]
HDI.a <- hdi(post.misha.pc2p5$BUGSoutput$sims.list$a,0.89)
HDI.a$CI_low
HDI.a$CI_high

plot(0,0, xlim = c(20000,8000), ylim = c(0.705, 0.711), xlab = "distance", ylab ="Sr 87/86",
     main="Calibration")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
#4000 lines are too many
#further thinning to 2000 lines
ind.pc2p<- sample(dim(post.misha.pc2p5$BUGSoutput$sims.list$R1.m)[1],500,replace = F)
MCMC.dist.plot(post.misha.pc2p5$BUGSoutput$sims.list$R1.m[ind.pc2p,],
               post.misha.pc2p5$BUGSoutput$sims.list$dist[ind.pc2p,])
points(n.avg.misha.50.dist[index.50.anom.remv1], n.avg.misha.50.sr[index.50.anom.remv1],
       pch=18, col = "#00b4ffff")
points(n.avg.misha.50.dist[index.50.anom.remv2], n.avg.misha.50.sr[index.50.anom.remv2],
       pch=18, col = "#00b4ffff")

lines(n.avg.misha.50.dist[index.50.anom.remv1], n.avg.misha.50.sr[index.50.anom.remv1],
      lwd=1.5, col = "#00b4ffff")
lines(n.avg.misha.50.dist[index.50.anom.remv2], n.avg.misha.50.sr[index.50.anom.remv2],
      lwd=1.5, col = "#00b4ffff")

par(mfrow=c(3,1))
par(mar = c(4.1, 4.1, 2.1, 4.1))
plot(density(post.misha.pc2p5$BUGSoutput$sims.list$a, from = 0), xlim = c(0,0.05),ylim= c(0,100),
     lwd = 2, col = plot.col.6[3],main="a & b", xlab="Parameter estimates")
lines(density(post.misha.pc2p5$BUGSoutput$sims.list$b, from = 0),
      lwd = 2, col = plot.col.6[4])
legend(0, 100, c("a","b"),
       lwd = rep(2, 2), col=c(plot.col.6[3:4]))

plot(density(post.misha.pc2p5$BUGSoutput$sims.list$c, from = 0, to = 0.015), xlim = c(0,0.015),
     lwd = 2, col = plot.col.6[5],main="c", xlab="Parameter estimate")

###-3 Constraints on parameters b, and c###
R.sd.mea <- n.sd.misha.50.sr.rmv
dist.mea <- n.avg.misha.50.dist.rmv
R.mea <- n.avg.misha.50.sr.rmv
n.mea = length(n.avg.misha.50.sr.rmv)

R0 <- 0.706

#Re is the mean ratio of end value  
Re <- 0.711
s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30
#parameters to save
parameters <- c("a", "b","c", "Ivo.rate", "R1.m","R2.m","Rin","dist.index",
                "Sr.pre", "Re.mean", "R0.mean", "switch","dist",
                "flux.ratio", "pool.ratio")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 750, n.mea = n.mea, 
           R0 = R0, Re = Re, s.intv = s.intv, Ivo.rate.mean = Ivo.rate.mean,
           Ivo.rate.sd = Ivo.rate.sd, max.dist.mea = max.dist.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e3
n.burnin = 2e3
n.thin = 1

#Run it
post.misha.pc2p = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS wo intake model4.R", 
                                             parameters.to.save = parameters, 
                                             data = dat, n.chains=5, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 26 hours

save(post.misha.pc2p, file = "out/post.misha.pc2p.RData")
load("out/post.misha.pc2p.RData")

traplot(post.misha.pc2p,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.pc2p,parms = c("a", "b","c"))

plot(density(post.misha.pc2p$BUGSoutput$sims.list$Ivo.rate))
abline(v=14.7)

#convergence params are not good with constraints on parameter c
subset(post.misha.pc2p$BUGSoutput$summary,
       rownames(post.misha.pc2p$BUGSoutput$summary)=="a")
subset(post.misha.pc2p$BUGSoutput$summary,
       rownames(post.misha.pc2p$BUGSoutput$summary)=="b")
subset(post.misha.pc2p$BUGSoutput$summary,
       rownames(post.misha.pc2p$BUGSoutput$summary)=="c")

MAP.a <- map_estimate(post.misha.pc2p$BUGSoutput$sims.list$a)
MAP.a[1]
HDI.a <- hdi(post.misha.pc2p$BUGSoutput$sims.list$a,0.89)
HDI.a$CI_low
HDI.a$CI_high

plot(0,0, xlim = c(20000,8000), ylim = c(0.7056, 0.712), xlab = "distance", ylab ="Sr 87/86",
     main="Calibration")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
#4000 lines are too many
#further thinning to 2000 lines
ind.pc2p<- sample(dim(post.misha.pc2p$BUGSoutput$sims.list$R1.m)[1],500,replace = F)
MCMC.dist.plot(post.misha.pc2p$BUGSoutput$sims.list$R1.m[ind.pc2p,],
               post.misha.pc2p$BUGSoutput$sims.list$dist[ind.pc2p,])
points(n.avg.misha.50.dist[index.50.anom.remv1], n.avg.misha.50.sr[index.50.anom.remv1],
       pch=18, col = "#00b4ffff")
points(n.avg.misha.50.dist[index.50.anom.remv2], n.avg.misha.50.sr[index.50.anom.remv2],
       pch=18, col = "#00b4ffff")

lines(n.avg.misha.50.dist[index.50.anom.remv1], n.avg.misha.50.sr[index.50.anom.remv1],
      lwd=1.5, col = "#00b4ffff")
lines(n.avg.misha.50.dist[index.50.anom.remv2], n.avg.misha.50.sr[index.50.anom.remv2],
      lwd=1.5, col = "#00b4ffff")

par(mfrow=c(3,1))
par(mar = c(4.1, 4.1, 2.1, 4.1))
plot(density(post.misha.pc2p$BUGSoutput$sims.list$a, from = 0), xlim = c(0,0.05),ylim= c(0,100),
     lwd = 2, col = plot.col.6[3],main="a & b", xlab="Parameter estimates")
lines(density(post.misha.pc2p$BUGSoutput$sims.list$b, from = 0),
      lwd = 2, col = plot.col.6[4])
legend(0, 100, c("a","b"),
       lwd = rep(2, 2), col=c(plot.col.6[3:4]))

plot(density(post.misha.pc2p$BUGSoutput$sims.list$c, from = 0, to = 0.015), xlim = c(0,0.015),
     lwd = 2, col = plot.col.6[5],main="c", xlab="Parameter estimate")

#####################################################################################
################################      Section 3      ################################
#####################################################################################
##############sensitivity to time series structure######
############################# inversion of Misha using different error per step params#######
##!! sensitive to precision term in Misha's inversion !!##
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
#turnover parameter is sensitive to the ts precision 

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

###################################################################################################
########################inversion with precision 2e-7########
####to test sensitivity to the precision term#####
Ivo.rate.mean <- mean.wooller.rate #microns per day
Ivo.rate.sd <- sd.wooller.rate

R.sd.mea <- sub.mm.sim.sd.sr
dist.mea <- sub.mm.sim.avg.dist
R.mea <- sub.mm.sim.avg.sr
n.mea = length(sub.mm.sim.avg.sr)

s.intv <- mm.bwidth

max.dist.mea <- max(sub.mm.sim.avg.dist)+ 800 #add some distance before the simulation

#posterior samples of parameters from Misha calibration

a.post <- post.misha.pc2p3$BUGSoutput$sims.list$a[,1]
b.post <- post.misha.pc2p3$BUGSoutput$sims.list$b[,1]
c.post <- post.misha.pc2p3$BUGSoutput$sims.list$c[,1]
post.leng <- length(a.post)

parameters <- c("Ivo.rate","dist", "R1.m","Rin.m", "R2.m","a","b","c","exp.ab",
                "Rin.m.pre","a.m","b.m","c.m","Body.mass.m", "Body.mass", "Rin.m.cps.ac")

dat = list( s.intv = s.intv, max.dist.mea = max.dist.mea, post.leng=post.leng, 
            a.post=a.post, b.post=b.post, c.post=c.post,
            Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 480, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e3
n.burnin = 2e3
n.thin = 1#floor(n.iter-n.burnin)/500

#Run it
post.misha.invmamm.s = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param mamms.R", 
                                                  parameters.to.save = parameters, 
                                                  data = dat, n.chains=5, n.iter = n.iter, 
                                                  n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 13.5 hours

save(post.misha.invmamm.s, file = "out/post.misha.invmamm.s.RData")

post.misha.invmamm.s$BUGSoutput$summary

load("out/post.misha.invmamm.s.RData")

plot(density(post.misha.invmamm.s$BUGSoutput$sims.list$Ivo.rate))

#check prior vs posterior parameters
plot(density(post.misha.pc2p3$BUGSoutput$sims.list$a[,1]), col = "black", lwd = 2, type="l",
     xlim = c(0.01,0.1), xlab = "a", ylab= "density")
#lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$a), col = "blue", lwd = 2)
lines(density(post.misha.invmamm.s$BUGSoutput$sims.list$a.m), col = "red", lwd = 2)
#slight deviation from prior


subset(post.misha.invmamm.s$BUGSoutput$summary,
       rownames(post.misha.invmamm.s$BUGSoutput$summary)=="a.m")#
subset(post.misha.invmamm.s$BUGSoutput$summary,
       rownames(post.misha.invmamm.s$BUGSoutput$summary)=="Rin.m.cps.ac")#poor convergence on autocorrelation

plot(density(post.misha.invmamm.s$BUGSoutput$sims.list$exp.ab))
plot(density(post.misha.invmamm.s$BUGSoutput$sims.list$Rin.m.cps.ac))

#do the posterior of a.m and c.m 

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,450), ylim = c(0.705, 0.715), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(sub.mm.sim.avg.dist)+ 800-sub.mm.sim.avg.dist)/mean.wooller.rate,
       sub.mm.sim.avg.sr, pch= 18, col="#00b4ffff")
lines((max(sub.mm.sim.avg.dist)+ 800-sub.mm.sim.avg.dist)/mean.wooller.rate,
      sub.mm.sim.avg.sr, lwd= 1.5, col="#00b4ffff")
#estimated input series
MCMC.ts.Rin.m.invmamm.s.89<- MCMC.CI.bound(post.misha.invmamm.s$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:480,MCMC.ts.Rin.m.invmamm.s.89[[1]],lwd = 2, col = "magenta")
lines(1:480,MCMC.ts.Rin.m.invmamm.s.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:480,MCMC.ts.Rin.m.invmamm.s.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Measured ivory","Reconstructed input"),lwd = c(1.5, 2), col=c("#00b3ffff","magenta"))

plot(density(post.misha.pc2p3$BUGSoutput$sims.list$a, from = 0), xlim = c(0.01,0.05),ylim= c(0,120),
     lwd = 2, col = "red",main="a", xlab="parameter estimate")
lines(density(post.misha.invmamm.s$BUGSoutput$sims.list$a.m, from = 0),lwd = 2,col="blue")

plot(density(post.misha.pc2p3$BUGSoutput$sims.list$c, from = 0), xlim = c(0.001,0.01),ylim= c(0,500),
     lwd = 2, col = "red",main="c", xlab="parameter estimate")
lines(density(post.misha.invmamm.s$BUGSoutput$sims.list$c.m, from = 0),lwd = 2,col="blue")

legend(0.04,120, c("Calibration","Mammoth"),lwd = c(2, 2), col=c("red","blue"))


#####################################################################################
################################      Section 4      ################################
#####################################################################################
##### inversion of Misha with autocorrelation and normal error per step: tsrw#######
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

parameters <- c("Ivo.rate", "dist","Rin.m.pre","Rin.m.cps.ac",
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
post.misha.inv2p.tsrw = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param tsrw.R", 
                                                     parameters.to.save = parameters, 
                                                     data = dat, n.chains=5, n.iter = n.iter, 
                                                     n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 26 hours

save(post.misha.inv2p.tsrw, file = "out/post.misha.inv2p.tsrw.RData")

post.misha.inv2p.tsrw$BUGSoutput$summary

#comparing parameter a
plot(density(post.misha.pc2p3$BUGSoutput$sims.list$a))
lines(density(post.misha.inv2p.tsrw$BUGSoutput$sims.list$a),col="red") 

plot(density(post.misha.inv2p.tsrw$BUGSoutput$sims.list$Rin.m.cps.ac)) #almost no autocorrelation


#plotting reconstructed Rin.m history
par(mfrow=c(1,1))
plot(0,0, xlim = c(1,700), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(n.avg.misha.50.dist) + 30-n.avg.misha.50.dist.rmv)/14.7,n.avg.misha.50.sr.rmv, pch= 18, col="#00b3ffff")
#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.pc2p3.Rin.89<- MCMC.CI.bound(post.misha.pc2p3$BUGSoutput$sims.list$Rin, 0.89)
lines(1:750,post.misha.pc2p3.Rin.89[[1]],lwd = 2, col = "black")
lines(1:750,post.misha.pc2p3.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:750,post.misha.pc2p3.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

#estimated input series
MCMC.misha.inv2p.tsrwca.Rin.m.89 <- MCMC.CI.bound(post.misha.inv2p.tsrwca$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:750,MCMC.misha.inv2p.tsrw.Rin.m.89[[1]],lwd = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2p.tsrw.Rin.m.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2p.tsrw.Rin.m.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Model input","Reconstructed input"),lwd = c(2, 2), col=c("black","magenta"))
##the normal error structure produce much more constrained random walk histories
##Which deviates from natural settings with possible leaps in 87Sr/86Sr

###################################################################################################
#####different time series error with auto correlation#####
##### inversion of Misha with autocorrelation and cauchy error per step: tsrwca#######
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

parameters <- c("Ivo.rate", "dist","Rin.m.pre","Rin.m.cps.ac",
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
post.misha.inv2p.tsrwca = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param tsrwca.R", 
                                                        parameters.to.save = parameters, 
                                                        data = dat, n.chains=5, n.iter = n.iter, 
                                                        n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 26 hours

save(post.misha.inv2p.tsrwca, file = "out/post.misha.inv2p.tsrwca.RData")

post.misha.inv2p.tsrwca$BUGSoutput$summary

#comparing parameter a
plot(density(post.misha.pc2p3$BUGSoutput$sims.list$a))
lines(density(post.misha.inv2p.tsrwca$BUGSoutput$sims.list$a),col="red") 

plot(density(post.misha.inv2p.tsrwca$BUGSoutput$sims.list$Rin.m.cps.ac)) #almost no autocorrelation


#plotting reconstructed Rin.m history
par(mfrow=c(1,1))
plot(0,0, xlim = c(1,700), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(n.avg.misha.50.dist) + 30-n.avg.misha.50.dist.rmv)/14.7,n.avg.misha.50.sr.rmv, pch= 18, col="#00b3ffff")
#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.pc2p3.Rin.89<- MCMC.CI.bound(post.misha.pc2p3$BUGSoutput$sims.list$Rin, 0.89)
lines(1:750,post.misha.pc2p3.Rin.89[[1]],lwd = 2, col = "black")
lines(1:750,post.misha.pc2p3.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:750,post.misha.pc2p3.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

#estimated input series
MCMC.misha.inv2p.tsrwca.Rin.m.89 <- MCMC.CI.bound(post.misha.inv2p.tsrwca$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:750,MCMC.misha.inv2p.tsrwca.Rin.m.89[[1]],lwd = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2p.tsrwca.Rin.m.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2p.tsrwca.Rin.m.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Model input","Reconstructed input"),lwd = c(2, 2), col=c("black","magenta"))
##the uncertainty of the estimated input is slightly wider that without autocorrelation (fig 4a)
##but the general patterns are very similar

###################################################################################################
##############case study: Mammoth##############
######Mammoth inversion with normal per step error: tsrw ########

Ivo.rate.mean <- mean.wooller.rate #microns per day
Ivo.rate.sd <- sd.wooller.rate

R.sd.mea <- sub.mm.sim.sd.sr
dist.mea <- sub.mm.sim.avg.dist
R.mea <- sub.mm.sim.avg.sr
n.mea = length(sub.mm.sim.avg.sr)

s.intv <- mm.bwidth

max.dist.mea <- max(sub.mm.sim.avg.dist)+ 800 #add some distance before the simulation

#posterior samples of parameters from Misha calibration

a.post <- post.misha.pc2p3$BUGSoutput$sims.list$a[,1]
b.post <- post.misha.pc2p3$BUGSoutput$sims.list$b[,1]
c.post <- post.misha.pc2p3$BUGSoutput$sims.list$c[,1]
post.leng <- length(a.post)

parameters <- c("Ivo.rate","dist", "R1.m","Rin.m", "R2.m","a","b","c","exp.ab",
                "Rin.m.pre","a.m","b.m","c.m","Body.mass.m", "Body.mass")

dat = list( s.intv = s.intv, max.dist.mea = max.dist.mea, post.leng=post.leng, 
            a.post=a.post, b.post=b.post, c.post=c.post,
            Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 500, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e3
n.burnin = 2e3
n.thin = 1

#Run it
post.invmamm.tsrw = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param mamm tsrw.R", 
                                                      parameters.to.save = parameters, 
                                                      data = dat, n.chains=5, n.iter = n.iter, 
                                                      n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 10 hours

save(post.invmamm.tsrw, file = "out/post.invmamm.tsrw.RData")

post.invmamm.tsrw$BUGSoutput$summary

load("out/post.misha.invmamm.param.RData")

plot(density(post.invmamm.tsrw$BUGSoutput$sims.list$Ivo.rate))

#check prior vs posterior parameters
plot(density(post.misha.pc2p3$BUGSoutput$sims.list$a[,1]), col = "black", lwd = 2, type="l",
     xlim = c(0.01,0.1), xlab = "a", ylab= "density")
lines(density(post.invmamm.tsrw$BUGSoutput$sims.list$a.m), col = "red", lwd = 2)
#slight deviation from prior

plot(density(post.misha.pc2p3$BUGSoutput$sims.list$c[,1]), col = "black", lwd = 2, type="l",
     xlim = c(0,0.1), xlab = "c", ylab= "density")
lines(density(post.invmamm.tsrw$BUGSoutput$sims.list$c.m), col = "red", lwd = 2)
#slight deviation from prior

subset(post.invmamm.tsrw$BUGSoutput$summary,
       rownames(post.invmamm.tsrw$BUGSoutput$summary)=="a.m")
subset(post.invmamm.tsrw$BUGSoutput$summary,
       rownames(post.invmamm.tsrw$BUGSoutput$summary)=="c.m")

plot(density(post.invmamm.tsrw$BUGSoutput$sims.list$exp.ab))

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,450), ylim = c(0.706, 0.715), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(sub.mm.sim.avg.dist)+ 800-sub.mm.sim.avg.dist)/mean.wooller.rate,
       sub.mm.sim.avg.sr, pch= 18, col="#00b4ffff")
lines((max(sub.mm.sim.avg.dist)+ 800-sub.mm.sim.avg.dist)/mean.wooller.rate,
      sub.mm.sim.avg.sr, lwd= 1.5, col="#00b4ffff")
#estimated input series
MCMC.ts.Rin.m.invmamm.tsrw.89<- MCMC.CI.bound(post.invmamm.tsrw$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:500,MCMC.ts.Rin.m.invmamm.tsrw.89[[1]],lwd = 2, col = "magenta")
lines(1:500,MCMC.ts.Rin.m.invmamm.tsrw.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:500,MCMC.ts.Rin.m.invmamm.tsrw.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Measured ivory","Reconstructed input"),lwd = c(1.5, 2), col=c("#00b3ffff","magenta"))