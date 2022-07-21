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
misha <- read.csv("data/Misha ivory.csv")

R.sd.mea <- misha$sd
dist.mea <- misha$dist
R.mea <- misha$mean
n.mea = length(dist.mea)

#adding the model component of food and water mixture as intake#
intake <- read.csv("data/intake.csv")

Hay <- intake[which(intake$type=="H"),]
Pellet <- intake[which(intake$type=="P"),]
Supplement <- intake[which(intake$type=="S"),]
Water <- intake[which(intake$type=="W"),]

#estimate mean and sd for hay Sr87/86
Sr.hay <- enorm(Hay$X87Sr.86Sr)
Sr.hay.mean <- Sr.hay$parameters[1]
Sr.hay.sd <- Sr.hay$parameters[2]

#estimate mean and sd for pellet Sr87/86
Sr.pel <- enorm(Pellet$X87Sr.86Sr)
Sr.pel.mean <- Sr.pel$parameters[1]
Sr.pel.sd <- Sr.pel$parameters[2]

#estimate mean and sd for alfalfa Sr87/86
Sr.sup <- enorm(Supplement$X87Sr.86Sr)
Sr.sup.mean <- Sr.sup$parameters[1]
Sr.sup.sd <- Sr.sup$parameters[2]

#estimate mean and sd for water Sr87/86
Sr.w <- enorm(Water$X87Sr.86Sr)
Sr.w.mean <- Sr.w$parameters[1]
Sr.w.sd <- Sr.w$parameters[2]

#estimate mean and sd for hay Sr concentration
#log-normal distribution
conc.hay <- elnorm(Hay$Sr_conc)
conc.hay.mean <- conc.hay$parameters[1]
conc.hay.sd <- conc.hay$parameters[2]

conc.sup <- elnorm(Supplement$Sr_conc)
conc.sup.mean <- conc.sup$parameters[1]
conc.sup.sd <- conc.sup$parameters[2]

#estimate mean and sd for pellet concentration
#log-normal distribution
conc.pel <- elnorm(Pellet$Sr_conc)
conc.pel.mean <- conc.pel$parameters[1]
conc.pel.sd <- conc.pel$parameters[2]

#estimate mean and sd for water concentration
#log-normal distribution
conc.w <- elnorm(Water$Sr_conc)
conc.w.mean <- conc.w$parameters[1]
conc.w.sd <- conc.w$parameters[2]

##########one pool turnover model########

#R0 is the mean ratio of initial value 
R0 <- 0.7071

m.feed <- 40 #mixing of x kg of hay, this is to constrain Sr variability of hay

Ivo.rate.mean <- 14.7
samp.interval <- mean(dist.mea[1:(n.mea - 1)] - dist.mea[2:n.mea])
n.days.bef.aft <- trunc(samp.interval/2/Ivo.rate.mean) #this parameter has to be supplied to the model

#parameters to save
parameters <- c("a", "Ivo.rate", "Rs.m","Rin","mod.index","mod.dist","f.h","f.pel","h.l",
                "Sr.pre", "Ps", "Fin", "Re.mean", "R0.mean", "switch", "w.contrib", "h.contrib")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 1050, n.mea = n.mea, 
           R0 = R0, Sr.hay.mean = Sr.hay.mean, Sr.hay.sd = Sr.hay.sd, 
           Sr.pel.mean = Sr.pel.mean, Sr.pel.sd = Sr.pel.sd, 
           Sr.sup.mean = Sr.sup.mean, Sr.sup.sd = Sr.sup.sd,
           Sr.w.mean = Sr.w.mean, Sr.w.sd = Sr.w.sd, m.feed = m.feed,
           conc.hay.mean = conc.hay.mean, conc.hay.sd = conc.hay.sd, 
           conc.pel.mean = conc.pel.mean, conc.pel.sd = conc.pel.sd,
           conc.sup.mean = conc.sup.mean, conc.sup.sd = conc.sup.sd,
           conc.w.mean = conc.w.mean, conc.w.sd = conc.w.sd,
           n.days.bef.aft = n.days.bef.aft)

#Start time
t1 = proc.time()

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 2e3
n.thin = floor(n.iter-n.burnin)/1000

#Run it
post.misha.nb5 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS no bone v5.R", 
                                            parameters.to.save = parameters, 
                                            data = dat, n.chains=3, n.iter = n.iter, 
                                            n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 4 hours

post.misha.nb5$BUGSoutput$summary

save(post.misha.nb5, file = "out/post.misha.nb5.RData")

load("out/post.misha.nb5.RData")

plot(density(post.misha.nb5$BUGSoutput$sims.list$a))

summary(post.misha.nb5$BUGSoutput$sims.list$switch)

plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rs.index.nb5 <- MCMC.ts.dist(post.misha.nb5$BUGSoutput$sims.list$Rs.m, 
                                     post.misha.nb5$BUGSoutput$sims.list$mod.dist,
                                     post.misha.nb5$BUGSoutput$sims.list$mod.index)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
MCMC.ts.RS.89.nb5 <- MCMC.ts(MCMC.ts.Rs.index.nb5)

lines(dist.mea, MCMC.ts.RS.89.nb5[[1]], lwd = 1, col = "cyan")
lines(dist.mea, MCMC.ts.RS.89.nb5[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, MCMC.ts.RS.89.nb5[[3]], lwd = 1, lty = 2, col = "cyan")

plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rin.index.nb5 <- MCMC.ts.dist(post.misha.nb5$BUGSoutput$sims.list$Rin, 
                                      post.misha.nb5$BUGSoutput$sims.list$mod.dist,
                                      post.misha.nb5$BUGSoutput$sims.list$mod.index)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
MCMC.ts.Rin.89.nb5 <- MCMC.ts(MCMC.ts.Rin.index.nb5)

lines(dist.mea, MCMC.ts.Rin.89.nb5[[1]], lwd = 1, col = "cyan")
lines(dist.mea, MCMC.ts.Rin.89.nb5[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, MCMC.ts.Rin.89.nb5[[3]], lwd = 1, lty = 2, col = "cyan")

plot(0,0, xlim = c(20000,8000), ylim = c(0, 1), xlab = "distance", ylab ="r")
MCMC.ts.fh.index.nb5 <- MCMC.ts.dist(post.misha.nb5$BUGSoutput$sims.list$f.h, 
                                     post.misha.nb5$BUGSoutput$sims.list$mod.dist,
                                     post.misha.nb5$BUGSoutput$sims.list$mod.index)

plot(0,0, xlim = c(20000,8000), ylim = c(0, 1), xlab = "distance", ylab ="r")
MCMC.ts.fpel.index.nb5 <- MCMC.ts.dist(post.misha.nb5$BUGSoutput$sims.list$f.pel, 
                                       post.misha.nb5$BUGSoutput$sims.list$mod.dist,
                                       post.misha.nb5$BUGSoutput$sims.list$mod.index)

plot(density(post.misha.nb5$BUGSoutput$sims.list$w.contrib))
plot(density(post.misha.nb5$BUGSoutput$sims.list$h.contrib))
