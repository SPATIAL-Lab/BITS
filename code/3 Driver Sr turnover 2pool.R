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

##################2 pool turnover with 50 pt average (excursion rmv)################
###prep data###
ind.remv.50 <- which(n.avg.misha.50.dist > 17000 & n.avg.misha.50.dist < 17500 & n.avg.misha.50.sr < 0.709)

index.50.anom.remv1 <- c(1:35)
index.50.anom.remv2 <- c(39:length(n.avg.misha.50.dist))
index.50.anom.remv <- c(1:35,39:length(n.avg.misha.50.dist))
n.avg.misha.50.dist.rmv <- n.avg.misha.50.dist[index.50.anom.remv]
n.avg.misha.50.sr.rmv <- n.avg.misha.50.sr[index.50.anom.remv]
n.sd.misha.50.sr.rmv <- n.sd.misha.50.sr[index.50.anom.remv]

###########version with the "new diet" excluded in the analysis###
##this is the version used in the publication####
ind.50.end.remv <- which(n.avg.misha.50.dist.rmv>10000)

n.sd.misha.50.sr.remv <- n.sd.misha.50.sr.rmv[ind.50.end.remv]
n.avg.misha.50.dist.remv <- n.avg.misha.50.dist.rmv[ind.50.end.remv]
n.avg.misha.50.sr.remv <-n.avg.misha.50.sr.rmv[ind.50.end.remv]

R.sd.mea <- n.sd.misha.50.sr.remv
dist.mea <- n.avg.misha.50.dist.remv
R.mea <- n.avg.misha.50.sr.remv
n.mea = length(n.avg.misha.50.sr.remv)

R0 <- 0.707

#Re is the mean ratio of end value  
Re <- 0.711
s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30
#parameters to save
parameters <- c("a", "b","c", "Ivo.rate", "R1.m","R2.m","Rin","dist.index",
                "Sr.pre", "Re.mean", "R0.mean", "switch","dist",
                "flux.ratio", "pool.ratio","Body.mass")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 700, n.mea = n.mea, 
           R0 = R0, Re = Re, s.intv = s.intv, Ivo.rate.mean = Ivo.rate.mean,
           Ivo.rate.sd = Ivo.rate.sd, max.dist.mea = max.dist.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e3
n.burnin = 2e3
n.thin = 1

#Run it
post.misha.pc2p.erm = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS wo intake model3.R", 
                                             parameters.to.save = parameters, 
                                             data = dat, n.chains=5, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 17 hours

save(post.misha.pc2p.erm, file = "out/post.misha.pc2p.erm.RData")
load("out/post.misha.pc2p.erm.RData")

traplot(post.misha.pc2p.erm,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.pc2p.erm,parms = c("a", "b","c"))

#convergence params are ok, but pool constraints on parameter c
subset(post.misha.pc2p.erm$BUGSoutput$summary,
       rownames(post.misha.pc2p.erm$BUGSoutput$summary)=="a")
subset(post.misha.pc2p.erm$BUGSoutput$summary,
       rownames(post.misha.pc2p.erm$BUGSoutput$summary)=="b")
subset(post.misha.pc2p.erm$BUGSoutput$summary,
       rownames(post.misha.pc2p.erm$BUGSoutput$summary)=="c")

MAP.a <- map_estimate(post.misha.pc2p.erm$BUGSoutput$sims.list$a)
MAP.a[1]
HDI.a <- hdi(post.misha.pc2p.erm$BUGSoutput$sims.list$a,0.89)
HDI.a$CI_low
HDI.a$CI_high

######version with the "new diet" included in the analysis#### 
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
n.iter = 5e3
n.burnin = 2e3
n.thin = floor(n.iter-n.burnin)/800

#Run it
post.misha.pc2p = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS wo intake model3.R", 
                                             parameters.to.save = parameters, 
                                             data = dat, n.chains=5, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 20 hours

save(post.misha.pc2p, file = "out/post.misha.pc2p.RData")
load("out/post.misha.pc2p.RData")

traplot(post.misha.pc2p,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.pc2p,parms = c("a", "b","c"))

#convergence params are not so good!
subset(post.misha.pc2p$BUGSoutput$summary,
       rownames(post.misha.pc2p$BUGSoutput$summary)=="a")
subset(post.misha.pc2p$BUGSoutput$summary,
       rownames(post.misha.pc2p$BUGSoutput$summary)=="b")
subset(post.misha.pc2p$BUGSoutput$summary,
       rownames(post.misha.pc2p$BUGSoutput$summary)=="c")
