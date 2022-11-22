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

###1 sensitivity to Re value################
#2 pool turnover with 50 pt average (excursion rmv) sensitivity test Re =0.7118##
###prep data###
R.sd.mea <- n.sd.misha.50.sr.rmv
dist.mea <- n.avg.misha.50.dist.rmv
R.mea <- n.avg.misha.50.sr.rmv
n.mea = length(n.avg.misha.50.sr.rmv)

R0 <- 0.707

#Re is the mean ratio of end value  
Re <- 0.7118
s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30
#parameters to save
parameters <- c("a", "b","c", "Ivo.rate", "Rs.m","Rb.m","Rin","dist.index",
                "Sr.pre", "Ps", "Fin","Pb", "Fb", "Re.mean", "R0.mean", "switch","dist",
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

MCMC.dist.plot(post.misha.pc2p8$BUGSoutput$sims.list$Rs.m,
               post.misha.pc2p8$BUGSoutput$sims.list$dist)
lines(n.avg.misha.50.dist.rmv,n.avg.misha.50.sr.rmv,lwd = 2, col = "#00b4ffff")


#compare calibration results with different priors#
plot(density(post.misha.pc2p$BUGSoutput$sims.list$a, from = 0), xlim = c(0,0.06),ylim= c(0,180),
     lwd = 2, col = plot.col.6[3],main="Two-pool model", xlab="parameter estimate")
lines(density(post.misha.pc2p$BUGSoutput$sims.list$b, from = 0),lwd = 2, col = plot.col.6[4])
lines(density(post.misha.pc2p$BUGSoutput$sims.list$c, from = 0),lwd = 2, col = plot.col.6[5])
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
R.sd.mea <- n.sd.misha.50.sr.rmv
dist.mea <- n.avg.misha.50.dist.rmv
R.mea <- n.avg.misha.50.sr.rmv
n.mea = length(n.avg.misha.50.sr.rmv)

R0 <- 0.707

#Re is the mean ratio of end value  
Re <- 0.7118
s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30

#Estimated date of switch
#by default, 83 is used, here we try 85, 81
#these values are expected to affect
switch <- 83 
#parameters to save
parameters <- c("a", "b","c", "Ivo.rate", "Rs.m","Rb.m","Rin","dist.index",
                "Sr.pre", "Ps", "Fin","Pb", "Fb", "Re.mean", "R0.mean", "switch","dist",
                "flux.ratio", "pool.ratio","Body.mass")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 750, n.mea = n.mea, 
           R0 = R0, Re = Re, s.intv = s.intv, Ivo.rate.mean = Ivo.rate.mean,
           Ivo.rate.sd = Ivo.rate.sd, max.dist.mea = max.dist.mea, switch = switch)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/400

#Run it
post.misha.woisw = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS woi sw.R", 
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

MCMC.dist.plot(post.misha.pc2p8$BUGSoutput$sims.list$Rs.m,
               post.misha.pc2p8$BUGSoutput$sims.list$dist)
lines(n.avg.misha.50.dist.rmv,n.avg.misha.50.sr.rmv,lwd = 2, col = "#00b4ffff")


#compare calibration results with different priors#
plot(density(post.misha.pc2p$BUGSoutput$sims.list$a, from = 0), xlim = c(0,0.06),ylim= c(0,180),
     lwd = 2, col = plot.col.6[3],main="Two-pool model", xlab="parameter estimate")
lines(density(post.misha.pc2p$BUGSoutput$sims.list$b, from = 0),lwd = 2, col = plot.col.6[4])
lines(density(post.misha.pc2p$BUGSoutput$sims.list$c, from = 0),lwd = 2, col = plot.col.6[5])
legend(0.02, 160, c("a, Two-pool","b, Two-pool", "c, Two-pool", "a, One-pool"),
       lwd = rep(2, 4), col=c(plot.col.6[3:5],plot.col.6[2]))

plot(density(post.misha.pc2p8$BUGSoutput$sims.list$a, from = 0), xlim = c(0,0.06),ylim= c(0,180),
     lwd = 2, col = plot.col.6[3],main="Two-pool model", xlab="parameter estimate")
lines(density(post.misha.pc2p8$BUGSoutput$sims.list$b, from = 0),lwd = 2, col = plot.col.6[4])
lines(density(post.misha.pc2p8$BUGSoutput$sims.list$c, from = 0),lwd = 2, col = plot.col.6[5])
legend(0.02, 160, c("a, Two-pool","b, Two-pool", "c, Two-pool", "a, One-pool"),
       lwd = rep(2, 4), col=c(plot.col.6[3:5],plot.col.6[2]))

##############sensitivity to time series structure######
##############use a random walk per time step structure#######
#eq: y(t) = B1*X(t-1) + e(t)
#but this time series structure triggers a run-away process,
#that increases the change per step by a lot
########################inversion by sampling posterior of 2-pool model########

Ivo.rate.mean <- mean.wooller.rate #microns per day
Ivo.rate.sd <- sd.wooller.rate

R.sd.mea <- sub.mm.sim.sd.sr
dist.mea <- sub.mm.sim.avg.dist
R.mea <- sub.mm.sim.avg.sr
n.mea = length(sub.mm.sim.avg.sr)

s.intv <- mm.bwidth

max.dist.mea <- max(sub.mm.sim.avg.dist)+ 800 #add some distance before the simulation

#posterior samples of parameters from Misha calibration

a.post <- post.misha.pc2p$BUGSoutput$sims.list$a[,1]
b.post <- post.misha.pc2p$BUGSoutput$sims.list$b[,1]
c.post <- post.misha.pc2p$BUGSoutput$sims.list$c[,1]
post.leng <- length(a.post)

parameters <- c("Ivo.rate","dist", "Rs.m","Rin.m", "Rb.m","a","b","c","exp.ab",
                "Rin.m.pre","a.m","b.m","c.m","Body.mass.m", "Body.mass")

dat = list( s.intv = s.intv, max.dist.mea = max.dist.mea, post.leng=post.leng, 
            a.post=a.post, b.post=b.post, c.post=c.post,
            Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 500, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 3e3
n.burnin = 1e3
n.thin = floor(n.iter-n.burnin)/500

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
plot(density(post.misha.pc2p$BUGSoutput$sims.list$a[,1]), col = "black", lwd = 2, type="l",
     xlim = c(0.01,0.1), xlab = "a", ylab= "density")
#lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$a), col = "blue", lwd = 2)
lines(density(post.invmamm.tsrw$BUGSoutput$sims.list$a.m), col = "red", lwd = 2)
#slight deviation from prior

plot(density(post.misha.pc2p$BUGSoutput$sims.list$c[,1]), col = "black", lwd = 2, type="l",
     xlim = c(0,0.1), xlab = "c", ylab= "density")
#lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$c), col = "blue", lwd = 2)
lines(density(post.invmamm.tsrw$BUGSoutput$sims.list$c.m), col = "red", lwd = 2)
#slight deviation from prior

subset(post.invmamm.tsrw$BUGSoutput$summary,
       rownames(post.invmamm.tsrw$BUGSoutput$summary)=="a.m")
subset(post.invmamm.tsrw$BUGSoutput$summary,
       rownames(post.invmamm.tsrw$BUGSoutput$summary)=="c.m")

plot(density(post.invmamm.tsrw$BUGSoutput$sims.list$exp.ab))


#do the posterior of a.m and c.m 

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,500), ylim = c(0.706, 0.715), xlab = "days", ylab ="Sr 87/86")
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