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

##################2 pool turnover with 50 pt average (excursion rmv)################
###prep data###
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
n.iter = 5e3
n.burnin = 1e3
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
subset(post.misha.pc2p$BUGSoutput$summary,
       rownames(post.misha.pc2p4$BUGSoutput$summary)=="flux.ratio")
subset(post.misha.pc2p$BUGSoutput$summary,
       rownames(post.misha.pc2p4$BUGSoutput$summary)=="pool.ratio")

summary(post.misha.pc2p$BUGSoutput$sims.list$switch)
summary(post.misha.pc2p$BUGSoutput$sims.list$Ivo.rate)

# make a contour map
contour.flux.pool <- kde2d(post.misha.pc2p$BUGSoutput$sims.list$flux.ratio[,1],
                           post.misha.pc2p$BUGSoutput$sims.list$pool.ratio[,1], n = 64,
                           lims = c(c(0,10), c(0,2.5)))
image(contour.flux.pool, col=viridis(64), xlab="Flux ratio: Fs/Fb", ylab="Pool ratio: Ps/Pb",xlim =c(0,10))
contour(contour.flux.pool,lwd = 1.5, add = TRUE, labcex = 1)

#save MAP estimates and 89% CI
MCMC.CI.a <- MCMC.CI.bound(post.misha.pc2p$BUGSoutput$sims.list$a, 0.89)
MCMC.CI.b <- MCMC.CI.bound(post.misha.pc2p$BUGSoutput$sims.list$b, 0.89)
MCMC.CI.c <- MCMC.CI.bound(post.misha.pc2p$BUGSoutput$sims.list$c, 0.89)

plot(density(post.misha.pc2p$BUGSoutput$sims.list$Ps), type = "l", lwd = 2, xlim = c(0, 12),
     ylim = c(0, 4), col = plot.col[2], xlab = "Pool size (mmol)", main = "")
lines(density(post.misha.pc2p$BUGSoutput$sims.list$Pb), lwd = 2, col = plot.col[6])
legend(8, 3, c("Serum","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

plot(density(post.misha.pc2p$BUGSoutput$sims.list$Fin), type = "l", lwd = 2, xlim = c(0, 0.06),
     ylim = c(0, 200), col = plot.col[2], xlab = "Sr Flux (mmol/day)", main = "")
lines(density(post.misha.pc2p$BUGSoutput$sims.list$Fb), lwd = 2, col = plot.col[6])
legend(0.04, 200, c("Intake","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

#a/b = Fin/Fb
#c/b = Ps/Pb

#plotting modeled serum values mapped onto ivory and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.pc2p$BUGSoutput$sims.list$Rs.m,
               post.misha.pc2p$BUGSoutput$sims.list$dist)
lines(n.avg.misha.50.dist.rmv,n.avg.misha.50.sr.rmv,lwd = 2, col = "#00b4ffff")

#plot bone pool values
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.pc2p$BUGSoutput$sims.list$Rb.m,
               post.misha.pc2p$BUGSoutput$sims.list$dist)
lines(n.avg.misha.50.dist.rmv,n.avg.misha.50.sr.rmv,lwd = 2, col = "#00b4ffff")

####check the Rin values####
plot(0,0, xlim = c(0,t), ylim = c(0.706, 0.713), xlab = "days", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
post.misha.pc2p.Rin.89<- MCMC.CI.bound(post.misha.pc2p$BUGSoutput$sims.list$Rin, 0.89)
lines(1:t, post.misha.pc2p.Rin.89[[1]], lwd = 2, col = "ff00ffff")
lines(1:t, post.misha.pc2p.Rin.89[[2]], lwd = 1, lty = 2, col = "ff00ffff")
lines(1:t, post.misha.pc2p.Rin.89[[3]], lwd = 1, lty = 2, col = "ff00ffff")

#reconstructed Rs.m history
plot(0,0, xlim = c(0,t), ylim = c(0.706, 0.713), xlab = "days", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.tl.plot(post.misha.pc2p$BUGSoutput$sims.list$Rs.m,750)
lines(1:750, post.misha.pc2p.Rs.m.89[[1]], lwd = 1, col = "cyan")
lines(1:750, post.misha.pc2p.Rs.m.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(1:750, post.misha.pc2p.Rs.m.89[[3]], lwd = 1, lty = 2, col = "cyan")

##################2 pool turnover with 50 pt average (excursion rmv) sensitivity test Re =0.7118########
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
proc.time() - t1 #~ 20 hours

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
subset(post.misha.pc2p8$BUGSoutput$summary,
       rownames(post.misha.pc2p4$BUGSoutput$summary)=="flux.ratio")
subset(post.misha.pc2p8$BUGSoutput$summary,
       rownames(post.misha.pc2p4$BUGSoutput$summary)=="pool.ratio")

# make a contour map
contour.flux.pool <- kde2d(post.misha.pc2p8$BUGSoutput$sims.list$flux.ratio[,1],
                           post.misha.pc2p8$BUGSoutput$sims.list$pool.ratio[,1], n = 64,
                           lims = c(c(0,10), c(0,2.5)))
image(contour.flux.pool, col=viridis(64), xlab="Flux ratio: Fs/Fb", ylab="Pool ratio: Ps/Pb",xlim =c(0,10))
contour(contour.flux.pool,lwd = 1.5, add = TRUE, labcex = 1)

#plotting some parameters
plot(density(post.misha.pc2p8$BUGSoutput$sims.list$Ps), type = "l", lwd = 2, xlim = c(0, 12),
     ylim = c(0, 4), col = plot.col[2], xlab = "Pool size (mmol)", main = "")
lines(density(post.misha.pc2p8$BUGSoutput$sims.list$Pb), lwd = 2, col = plot.col[6])
legend(8, 3, c("Serum","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

plot(density(post.misha.pc2p8$BUGSoutput$sims.list$Fin), type = "l", lwd = 2, xlim = c(0, 0.06),
     ylim = c(0, 200), col = plot.col[2], xlab = "Sr Flux (mmol/day)", main = "")
lines(density(post.misha.pc2p8$BUGSoutput$sims.list$Fb), lwd = 2, col = plot.col[6])
legend(0.04, 200, c("Intake","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

#a/b = Fin/Fb
#c/b = Ps/Pb

#plotting modeled serum values mapped onto ivory and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.pc2p8$BUGSoutput$sims.list$Rs.m,
               post.misha.pc2p8$BUGSoutput$sims.list$dist)
lines(n.avg.misha.50.dist.rmv,n.avg.misha.50.sr.rmv,lwd = 2, col = "#00b4ffff")

#plot bone pool values
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.pc2p8$BUGSoutput$sims.list$Rb.m,
               post.misha.pc2p8$BUGSoutput$sims.list$dist)
lines(n.avg.misha.50.dist.rmv,n.avg.misha.50.sr.rmv,lwd = 2, col = "#00b4ffff")

####check the Rin values####
plot(0,0, xlim = c(0,t), ylim = c(0.706, 0.713), xlab = "days", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
post.misha.pc2p4.Rin.89<- MCMC.CI.bound(post.misha.pc2p4$BUGSoutput$sims.list$Rin, 0.89)
lines(1:t, post.misha.pc2p4.Rin.89[[1]], lwd = 2, col = "ff00ffff")
lines(1:t, post.misha.pc2p4.Rin.89[[2]], lwd = 1, lty = 2, col = "ff00ffff")
lines(1:t, post.misha.pc2p4.Rin.89[[3]], lwd = 1, lty = 2, col = "ff00ffff")

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