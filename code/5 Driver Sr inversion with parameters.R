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

########inversion model with calibration for either the one/two pool model#######
#inversion based on micromill results
micromill <- read.csv("data/Misha micromill.csv")

R.mic <- rev(micromill$X87Sr.86Sr)
R.sd.mic <- rev(micromill$StdErr)
dist.mic <- rev(micromill$dist)
n.mic <- length(dist.mic)

#calculate average sampling interval
Ivo.rate.mic <- 20 #this is estimated for the M640 slab, according to Uno 2012
samp.interval <- mean(dist.mic[1:(n.mic - 1)] - dist.mic[2:n.mic])
n.days.bef.aft <- trunc(samp.interval/2/Ivo.rate.mic) #this parameter has to be supplied to the model

Ivo.rate.mean <- 14.7
cal.interval <- mean(dist.mea[1:(n.mea - 1)] - dist.mea[2:n.mea])
n.days.bef.aft.cal <- trunc(cal.interval/2/Ivo.rate.mean) #this parameter has to be supplied to the model


R0 <- 0.7071

#parameters to save
parameters <- c("Ivo.rate", "Rs.cal", "Rb.cal", "mod.index.cal","mod.dist.cal","Rin.cal","a","b","c",
                "Rs.m", "Rb.m","Rin.m","mod.index","mod.dist", "a.m", "b.m", "c.m","Rin.m.eps.ac",
                "Ps", "Fin","Pb", "Fb", "Re.mean", "switch", "w.contrib", "h.contrib",
                "flux.ratio", "pool.ratio")

##Data to pass to the model
#compared to the turnover model that is essentially the .cal part here 
#the inversion takes the measured value of potentially a different ivory series
dat = list(R.cal = R.mea, dist.cal = dist.mea, R.sd.cal = R.sd.mea, t.cal = 1400, n.cal = n.mea, 
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
n.iter = 5e3
n.burnin = 1e3
n.thin = floor(n.iter-n.burnin)/600

#Run it
post.misha.inversion4 = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS v4.R", 
                                                   parameters.to.save = parameters, 
                                                   data = dat, n.chains=5, n.iter = n.iter, 
                                                   n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 8 hours

save(post.misha.inversion4, file = "out/post.misha.inversion4.RData")

post.misha.inversion4$BUGSoutput$summary

load("out/post.misha.inversion4.RData")

#plotting modeled serum values and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rs.cal.index.inversion4 <- MCMC.ts.dist(post.misha.inversion4$BUGSoutput$sims.list$Rs.cal, 
                                                post.misha.inversion4$BUGSoutput$sims.list$mod.dist.cal,
                                                post.misha.inversion4$BUGSoutput$sims.list$mod.index.cal)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
MCMC.ts.Rs.89.inversion4 <- MCMC.ts(MCMC.ts.Rs.cal.index.inversion4)

lines(dist.mea, MCMC.ts.Rs.89.inversion4[[1]], lwd = 1, col = "cyan")
lines(dist.mea, MCMC.ts.Rs.89.inversion4[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, MCMC.ts.Rs.89.inversion4[[3]], lwd = 1, lty = 2, col = "cyan")
legend(12000, 0.709, c("measured ivory","modeled serum"),lwd = c(2, 1), col=c("red","cyan"))
#end plot

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,t), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to days using rate Ivo.rate.mic
lines(dist.mic/Ivo.rate.mic, R.mic, lwd = 2, col = "blue") #results from micromill

MCMC.ts.Rin.m.89.4 <- MCMC.ts(post.misha.inversion4$BUGSoutput$sims.list$Rin.m)
lines(1:t,MCMC.ts.Rin.m.89.4[[1]],lwd = 2, col = "firebrick4")
lines(1:t,MCMC.ts.Rin.m.89.4[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:t,MCMC.ts.Rin.m.89.4[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))