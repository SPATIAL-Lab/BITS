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


###Version 4 of one pool model######


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

#####no bone but with time averaging parameter #this should be used in inversion model
#################################
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

########inversion model#######
#####v4 inversion for Rin using posterior distributions of a, b and c####
#inversion based on micromill results
micromill <- read.csv("data/Misha micromill.csv")

R.mic <- micromill$X87Sr.86Sr
R.sd.mic <- micromill$StdErr
dist.mic <- micromill$dist
n.mic <- length(dist.mic)

#calculate average sampling interval
Ivo.rate.mic <- 20 #this is estimated for M640, according to Uno 2012
samp.interval <- mean(dist.mic[1:(n.mic - 1)] - dist.mic[2:n.mic])
n.days.bef.aft <- trunc(samp.interval/2/Ivo.rate.mic) #this parameter has to be supplied to the model

Ivo.rate.mean <- 14.7
samp.interval <- mean(dist.mea[1:(n.mea - 1)] - dist.mea[2:n.mea])
n.days.bef.aft.cal <- trunc(samp.interval/2/Ivo.rate.mean) #this parameter has to be supplied to the model


R0 <- 0.70706

#parameters to save
parameters <- c("Ivo.rate", "Rs.cal", "Rb.cal", "mod.index.cal","mod.dist.cal","Rin.cal","a","b","c",
                "Rs.m", "Rb.m","Rin.m","mod.index","mod.dist", "a.m", "b.m", "c.m","Rin.m.eps.ac")

##Data to pass to the model
#compared to the turnover model that is essentially the .cal part here 
#the inversion takes the measured value of potentially a different ivory series
dat = list(R.cal = R.mea, dist.cal = dist.mea, R.sd.cal = R.sd.mea, t.cal = 1200, n.cal = n.mea, 
           R0 = R0,  
           Sr.hay.mean = Sr.hay.mean, Sr.hay.sd = Sr.hay.sd, 
           Sr.pel.mean = Sr.pel.mean, Sr.pel.sd = Sr.pel.sd, 
           Sr.sup.mean = Sr.sup.mean, Sr.sup.sd = Sr.sup.sd,
           Sr.w.mean = Sr.w.mean, Sr.w.sd = Sr.w.sd, m.feed = m.feed,
           conc.hay.mean = conc.hay.mean, conc.hay.sd = conc.hay.sd, 
           conc.pel.mean = conc.pel.mean, conc.pel.sd = conc.pel.sd,
           conc.sup.mean = conc.sup.mean, conc.sup.sd = conc.sup.sd,
           conc.w.mean = conc.w.mean, conc.w.sd = conc.w.sd,
           n.days.bef.aft = n.days.bef.aft, n.days.bef.aft.cal = n.days.bef.aft.cal,
           R.mea = R.mic, dist.mea = dist.mic, R.sd.mea = R.sd.mic, t = 800, n.mea = n.mic)

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

#Notes on generating appositional weight sequence
#create vector of distance (in micron), represented 

#Convert into intervals that are rounded into integer
intv.int[t] = ifelse(intv[i] - 0.5 > trunc(intv[i]), round(intv[i]), trunc(intv[i]))
intv.int[1:t-1] = dist.int[2:t] - dist.int[1:t-1]
#rounding distance
for(i in 1:t){
  dist.int[i] = ifelse(dist[i] - 0.5 > trunc(dist[i]), round(dist[i]), trunc(dist[i]))
}

##generating sequence, t = number of weeks
for(i in 2:t){
  dist[i] = dist[i-1] + intv[i] * 7
  intv[i] ~ dnorm(ext.rate.m, ext.rate.pre)
}
dist[1] = intv[1] * 7 #weekly resolution
intv[1] ~ dnorm(ext.rate.m, ext.rate.pre)

en.leng = 1:(l.en + la + lm) #this will have to be supplied to the 

lm.mean = 69800
lm.sd = 4800
ET = 3100

alpha = 3.2/180*pi #appositional angle ~ 3.2 degrees

la.mean = 55600
la.sd = 2000

ext.rate.m = 4.5
ext.rate.pre = 1/0.1^2

fi.mean = 0.65
fi.sd = 0.05

DER = 41

#calculate weights and sum of products 

for(i in 1:n.mea){
  for (j in 1:t){
    ifelse(i == mod.index[j], weights[i,j] = 1/(mod.dist[j] - dist.mea[i])^2, weights[i,j] = 0)
    
    prod.w.Rs.m[i,j] = weights[i,j] * Rs.m[j]
    
  }
  sum.w[i] = sum(weights[i,])
  sum.prod.w.Rs.m[i] = sum(prod.w.Rs.m[i,])
}