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

##########two pool turnover model########
#####testing version 4 of the turnover model with evaluation on the mean of neighbouring data points##
#and variable R0 and Re values###
#and variable switch positions##
R.sd.mea <- misha$sd
dist.mea <- misha$dist
R.mea <- misha$mean
n.mea = length(dist.mea)
#R0 is the mean ratio of initial value 
R0 <- 0.7071
#Re is the mean ratio of end value  
Re <- 0.7112

#hist(1/rgamma(1000, shape = 10, rate = 0.2)^0.5)

#parameters to save
parameters <- c("a", "b", "c", "Ivo.rate", "Rs.m", "Rb.m","Rin","mod.index","mod.dist",
                "flux.ratio", "pool.ratio", "Sr.pre", "Ps", "Pb", "Fin", "Fb", "switch")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 1050, n.mea = n.mea, 
           R0 = R0, Re = Re)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e4
n.burnin = 4e3
n.thin = floor(n.iter-n.burnin)/1000

#Run it
post.misha.4 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS v4.R", 
                                        parameters.to.save = parameters, 
                                        data = dat, n.chains=3, n.iter = n.iter, 
                                        n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 5 hours

post.misha.4$BUGSoutput$summary

save(post.misha.4, file = "out/post.misha.4.RData")
load("out/post.misha.4.RData")

traplot(post.misha.4,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.4,parms = c("a", "b","c"))

summary(post.misha.4$BUGSoutput$sims.list$switch)
summary(post.misha.4$BUGSoutput$sims.list$Ivo.rate)

plot(density(post.misha.4$BUGSoutput$sims.list$a))


hist(post.misha.4$BUGSoutput$sims.list$flux.ratio)
hist(post.misha.4$BUGSoutput$sims.list$pool.ratio)
#make a contour map
contour.flux.pool <- kde2d(post.misha.4$BUGSoutput$sims.list$flux.ratio[,1], 
                           post.misha.4$BUGSoutput$sims.list$pool.ratio[,1], n = 64)
image(contour.flux.pool, col=viridis(64), xlab="Flux ratio: Fs/Fb", ylab="Pool ratio: Ps/Pb",xlim =c(0,200))
contour(contour.flux.pool,lwd = 1.5, add = TRUE, labcex = 1)

#plotting some parameters
hist(post.misha.4$BUGSoutput$sims.list$Ivo.rate)
map_estimate(post.misha.4$BUGSoutput$sims.list$Pb) #3.24 mmol
hist(post.misha.4$BUGSoutput$sims.list$Pb)
map_estimate(post.misha.4$BUGSoutput$sims.list$Ps) #1.05 mmol
map_estimate(post.misha.4$BUGSoutput$sims.list$Fin) #0.03 mmol/day
hist(post.misha.4$BUGSoutput$sims.list$Fin)
map_estimate(post.misha.4$BUGSoutput$sims.list$Fb) #0.00736 mmol/day
hist(post.misha.4$BUGSoutput$sims.list$Fb)

plot(density(post.misha.4$BUGSoutput$sims.list$Ps), type = "l", lwd = 2, xlim = c(0, 12),
     ylim = c(0, 4), col = plot.col[2], xlab = "Pool size (mmol)", main = "")
lines(density(post.misha.4$BUGSoutput$sims.list$Pb), lwd = 2, col = plot.col[6])
legend(8, 3, c("Serum","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

plot(density(post.misha.4$BUGSoutput$sims.list$Fin), type = "l", lwd = 2, xlim = c(0, 0.1),
     ylim = c(0, 80), col = plot.col[2], xlab = "Sr Flux (mmol/day)", main = "")
lines(density(post.misha.4$BUGSoutput$sims.list$Fb), lwd = 2, col = plot.col[6])
legend(0.06, 80, c("Intake","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

#a/b = Fin/Fb
#c/b = Ps/Pb

#plotting modeled serum values and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
post.misha.4.Rs.index <- MCMC.ts.dist(post.misha.4$BUGSoutput$sims.list$Rs.m, 
                                 post.misha.4$BUGSoutput$sims.list$mod.dist,
                                 post.misha.4$BUGSoutput$sims.list$mod.index)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
post.misha.4.RS.89 <- MCMC.ts(post.misha.4.Rs.index)

lines(dist.mea, post.misha.4.RS.89[[1]], lwd = 1, col = "cyan")
lines(dist.mea, post.misha.4.RS.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, post.misha.4.RS.89[[3]], lwd = 1, lty = 2, col = "cyan")

#plotting modeled bone values and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
post.misha.4.Rb.index <- MCMC.ts.dist(post.misha.4$BUGSoutput$sims.list$Rb.m, 
                                      post.misha.4$BUGSoutput$sims.list$mod.dist,
                                      post.misha.4$BUGSoutput$sims.list$mod.index)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
post.misha.4.Rb.89 <- MCMC.ts(post.misha.4.Rb.index)

lines(dist.mea, post.misha.4.Rb.89[[1]], lwd = 1, col = "cyan")
lines(dist.mea, post.misha.4.Rb.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, post.misha.4.Rb.89[[3]], lwd = 1, lty = 2, col = "cyan")

####check the Rin values####
plot(0,0, xlim = c(0,1050), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.misha.4.Rin <- MCMC.ts(post.misha.4$BUGSoutput$sims.list$Rin)
lines(1:1050,MCMC.ts.misha.4.Rin[[1]],lwd = 2, col = "black")
lines(1:1050,MCMC.ts.misha.4.Rin[[2]], lwd = 1, lty = 2, col = "black")
lines(1:1050,MCMC.ts.misha.4.Rin[[3]], lwd = 1, lty = 2, col = "black")

plot(density(post.misha.4$BUGSoutput$sims.list$Ps), type = "l", lwd = 2, xlim = c(0, 12),
     ylim = c(0, 3), col = plot.col[2], xlab = "Pool size (mmol)", main = "")
lines(density(post.misha.4$BUGSoutput$sims.list$Pb), lwd = 2, col = plot.col[6])
legend(8, 3, c("Serum","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

plot(density(post.misha.4$BUGSoutput$sims.list$Fin), type = "l", lwd = 2, xlim = c(0, 0.4),
     ylim = c(0, 20), col = plot.col[2], xlab = "Sr Flux (mmol/day)", main = "")
lines(density(post.misha.4$BUGSoutput$sims.list$Fb), lwd = 2, col = plot.col[6])
legend(0.3, 20, c("Intake","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

###Version 4 of one pool model######
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

plot(seq(0, 1, length=100), dbeta(seq(0, 1, length=100), 8, 32))

#daily dry matter intake by adults of about 1.5% of body weight

R0 <- 0.7071

m.feed <- 40 #mixing of x kg of hay, this is to constrain Sr variability of hay

#parameters to save
parameters <- c("a", "Ivo.rate", "Rs.m","Rin","mod.index","mod.dist","f.h","f.pel",
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
           conc.w.mean = conc.w.mean, conc.w.sd = conc.w.sd)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 2e3
n.thin = floor(n.iter-n.burnin)/1000

#Run it
post.misha.nb4 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS no bone v4.R", 
                                            parameters.to.save = parameters, 
                                            data = dat, n.chains=3, n.iter = n.iter, 
                                            n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 4 hours

post.misha.nb4$BUGSoutput$summary

save(post.misha.nb4, file = "out/post.misha.nb4.RData")

load("out/post.misha.nb4.RData")

plot(density(post.misha.nb4$BUGSoutput$sims.list$a))

summary(post.misha.nb4$BUGSoutput$sims.list$switch)

plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rs.index.nb4 <- MCMC.ts.dist(post.misha.nb4$BUGSoutput$sims.list$Rs.m, 
                                     post.misha.nb4$BUGSoutput$sims.list$mod.dist,
                                     post.misha.nb4$BUGSoutput$sims.list$mod.index)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
MCMC.ts.RS.89.nb4 <- MCMC.ts(MCMC.ts.Rs.index.nb4)

lines(dist.mea, MCMC.ts.RS.89.nb4[[1]], lwd = 1, col = "cyan")
lines(dist.mea, MCMC.ts.RS.89.nb4[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, MCMC.ts.RS.89.nb4[[3]], lwd = 1, lty = 2, col = "cyan")

plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Re.index.nb4 <- MCMC.ts.dist(post.misha.nb4$BUGSoutput$sims.list$Re.mean, 
                                     post.misha.nb4$BUGSoutput$sims.list$mod.dist,
                                     post.misha.nb4$BUGSoutput$sims.list$mod.index)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
MCMC.ts.Re.89.nb4 <- MCMC.ts(MCMC.ts.Re.index.nb4)

lines(dist.mea, MCMC.ts.Re.89.nb4[[1]], lwd = 1, col = "cyan")
lines(dist.mea, MCMC.ts.Re.89.nb4[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, MCMC.ts.Re.89.nb4[[3]], lwd = 1, lty = 2, col = "cyan")

plot(0,0, xlim = c(20000,8000), ylim = c(0, 1), xlab = "distance", ylab ="r")
MCMC.ts.fr.index.nb4 <- MCMC.ts.dist(post.misha.nb4$BUGSoutput$sims.list$f.h, 
                                     post.misha.nb4$BUGSoutput$sims.list$mod.dist,
                                     post.misha.nb4$BUGSoutput$sims.list$mod.index)

plot(density(post.misha.nb4$BUGSoutput$sims.list$f.h))
plot(density(post.misha.nb4$BUGSoutput$sims.list$f.pel))

#####ver 5 no bone with variable food ratios but fixed Sr values
R0 <- 0.7071

m.feed <- 40 #mixing of x kg of hay, this is to constrain Sr variability of hay

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
           conc.w.mean = conc.w.mean, conc.w.sd = conc.w.sd)

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

#####no bone but with time averaging parameter
R0 <- 0.7071

m.feed <- 40 #mixing of x kg of hay, this is to constrain Sr variability of hay

#calculate n.days before and after the day of evaluation for time averaging
#this number is truncated for more conservative evaluation

#calculate average sampling interval
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

set.seed(t1[3])
n.iter = 5e3
n.burnin = 1e3
n.thin = floor(n.iter-n.burnin)/1000

#Run it
post.misha.nbm = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS no bone micromill.R", 
                                            parameters.to.save = parameters, 
                                            data = dat, n.chains=3, n.iter = n.iter, 
                                            n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 4 hours

post.misha.nbm$BUGSoutput$summary

save(post.misha.nbm, file = "out/post.misha.nbm.RData")

load("out/post.misha.nbm.RData")

plot(density(post.misha.nbm$BUGSoutput$sims.list$a))

summary(post.misha.nbm$BUGSoutput$sims.list$switch)

plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rs.index.nbm <- MCMC.ts.dist(post.misha.nbm$BUGSoutput$sims.list$Rs.m, 
                                     post.misha.nbm$BUGSoutput$sims.list$mod.dist,
                                     post.misha.nbm$BUGSoutput$sims.list$mod.index)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
MCMC.ts.RS.89.nbm <- MCMC.ts(MCMC.ts.Rs.index.nbm)

lines(dist.mea, MCMC.ts.RS.89.nbm[[1]], lwd = 1, col = "cyan")
lines(dist.mea, MCMC.ts.RS.89.nbm[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, MCMC.ts.RS.89.nbm[[3]], lwd = 1, lty = 2, col = "cyan")

plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rin.index.nbm <- MCMC.ts.dist(post.misha.nbm$BUGSoutput$sims.list$Rin, 
                                      post.misha.nbm$BUGSoutput$sims.list$mod.dist,
                                      post.misha.nbm$BUGSoutput$sims.list$mod.index)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
MCMC.ts.Rin.89.nbm <- MCMC.ts(MCMC.ts.Rin.index.nbm)

lines(dist.mea, MCMC.ts.Rin.89.nbm[[1]], lwd = 1, col = "cyan")
lines(dist.mea, MCMC.ts.Rin.89.nbm[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, MCMC.ts.Rin.89.nbm[[3]], lwd = 1, lty = 2, col = "cyan")

########inversion model#######
#####v3 inversion for Rin using posterior distributions of a, b and c####
#assign input values before and after the switch, but allows some variation
R0 <- 0.70706

#Re is the mean ratio of end value  
Re <- 0.71179

switch <- 75

#parameters to save
parameters <- c("Ivo.rate", "Rs.cal", "Rb.cal", "mod.index.cal","mod.dist.cal","Rin.cal","a","b","c",
                "Rs.m", "Rb.m","Rin.m","mod.index","mod.dist", "a.m", "b.m", "c.m","Rin.m.eps.ac")

##Data to pass to the model
#compared to the turnover model that is essentially the .cal part here 
#the inversion takes the measured value of potentially a different ivory series
dat = list(R.cal = R.mea, dist.cal = dist.mea, R.sd.cal = R.sd.mea, t.cal = 1200, n.cal = n.mea, 
           switch = switch, R0 = R0, Re = Re, 
           R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 1200, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 2e3
n.thin = floor(n.iter-n.burnin)/600

#Run it
post.misha.inversion3 = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS v3.R", 
                                                   parameters.to.save = parameters, 
                                                   data = dat, n.chains=5, n.iter = n.iter, 
                                                   n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 8 hours

save(post.misha.inversion3, file = "out/post.misha.inversion3.RData")

post.misha.inversion3$BUGSoutput$summary

post.misha.inversion3 <- load("out/post.misha.inversion3.RData")

#plotting modeled serum values and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rs.cal.index.inversion3 <- MCMC.ts.dist(post.misha.inversion3$BUGSoutput$sims.list$Rs.cal, 
                                                post.misha.inversion3$BUGSoutput$sims.list$mod.dist.cal,
                                                post.misha.inversion3$BUGSoutput$sims.list$mod.index.cal)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
MCMC.ts.Rs.89.inversion3 <- MCMC.ts(MCMC.ts.Rs.cal.index.inversion3)

lines(dist.mea, MCMC.ts.Rs.89.inversion3[[1]], lwd = 1, col = "cyan")
lines(dist.mea, MCMC.ts.Rs.89.inversion3[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, MCMC.ts.Rs.89.inversion3[[3]], lwd = 1, lty = 2, col = "cyan")
legend(12000, 0.709, c("measured ivory","modeled serum"),lwd = c(2, 1), col=c("red","cyan"))
#end plot

#plotting reconstructed Rin and compared that with posterior of Rin
plot(0,0, xlim = c(1,1000), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rin.cal.89.3 <- MCMC.ts(post.misha.inversion3$BUGSoutput$sims.list$Rin.cal)
lines(1:1200,MCMC.ts.Rin.cal.89.3[[1]],lwd = 2, col = "black")
lines(1:1200,MCMC.ts.Rin.cal.89.3[[2]], lwd = 1, lty = 2, col = "black")
lines(1:1200,MCMC.ts.Rin.cal.89.3[[3]], lwd = 1, lty = 2, col = "black")

MCMC.ts.Rin.m.89.3 <- MCMC.ts(post.misha.inversion3$BUGSoutput$sims.list$Rin.m)
lines(1:1200,MCMC.ts.Rin.m.89.3[[1]],lwd = 2, col = "red")
lines(1:1200,MCMC.ts.Rin.m.89.3[[2]], lwd = 1, lty = 2, col = "red")
lines(1:1200,MCMC.ts.Rin.m.89.3[[3]], lwd = 1, lty = 2, col = "red")
legend(0, 0.716, c("PD input","PD reconstructed"),lwd = c(2, 2), col=c("black","red"))
#end plot
