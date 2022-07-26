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

########inversion model without calibration, parameters taken from the two pool model#######
#inversion based on micromill results

misha <- read.csv("data/Misha ivory.csv")

R.sd.cal <- misha$sd
dist.cal <- misha$dist
R.cal <- misha$mean
n.cal = length(dist.mea)

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

#assign input values before and after the switch, but allows some variation
R0 <- 0.70706

#Re is the mean ratio of end value  
Re <- 0.7112

#parameters to save
parameters <- c("Ivo.rate", "Rs.cal", "Rb.cal", "mod.index.cal","mod.dist.cal","Rin.cal","a","b","c",
                "Rs.m", "Rb.m","Rin.m","mod.index","mod.dist", "a.m", "b.m", "c.m","Rin.m.eps.ac")

##Data to pass to the model
#compared to the turnover model that is essentially the .cal part here 
#the inversion takes the measured value of potentially a different ivory series
dat = list(R.cal = R.cal, dist.cal = dist.cal, R.sd.cal = R.sd.cal, t.cal = 1100, n.cal = n.cal, 
           R0 = R0, Re = Re,
           n.days.bef.aft = n.days.bef.aft, n.days.bef.aft.cal = n.days.bef.aft.cal,
           R.mea = R.mic, dist.mea = dist.mic, R.sd.mea = R.sd.mic, t = 500, n.mea = n.mic)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 2e3
n.thin = floor(n.iter-n.burnin)/600

#Run it
post.misha.inversion.p = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS wo intake model.R", 
                                                   parameters.to.save = parameters, 
                                                   data = dat, n.chains=5, n.iter = n.iter, 
                                                   n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 8 hours

save(post.misha.inversion.p, file = "out/post.misha.inversion.p.RData")

post.misha.inversion.p$BUGSoutput$summary

load("out/post.misha.inversion.p.RData")

plot(density(post.misha.inversion.p$BUGSoutput$sims.list$a))

plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rs.c.index.inversion.p <- MCMC.ts.dist(post.misha.inversion.p$BUGSoutput$sims.list$Rs.cal, 
                                               post.misha.inversion.p$BUGSoutput$sims.list$mod.dist.cal,
                                               post.misha.inversion.p$BUGSoutput$sims.list$mod.index.cal)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
MCMC.ts.Rs.c89.inversion.p <- MCMC.ts(MCMC.ts.Rs.c.index.inversion.p)

lines(misha$dist, MCMC.ts.Rs.c89.inversion.p[[1]], lwd = 1, col = "cyan")
lines(misha$dist, MCMC.ts.Rs.c89.inversion.p[[2]], lwd = 1, lty = 2, col = "cyan")
lines(misha$dist, MCMC.ts.Rs.c89.inversion.p[[3]], lwd = 1, lty = 2, col = "cyan")

#plotting modeled serum values and checking the fit of the data
plot(0,0, xlim = c(10000,0), ylim = c(0.709, 0.714), xlab = "distance", ylab ="Sr 87/86")
MCMC.ts.Rs.m.index.inversion.p <- MCMC.ts.dist(post.misha.inversion.p$BUGSoutput$sims.list$Rs.m, 
                                                post.misha.inversion.p$BUGSoutput$sims.list$mod.dist,
                                                post.misha.inversion.p$BUGSoutput$sims.list$mod.index)
lines(dist.mic,R.mic,lwd = 2, col = "red")
MCMC.ts.Rs.m89.inversion.p <- MCMC.ts(MCMC.ts.Rs.m.index.inversion.p)

lines(dist.mea, MCMC.ts.Rs.m89.inversion.p[[1]], lwd = 1, col = "cyan")
lines(dist.mea, MCMC.ts.Rs.m89.inversion.p[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, MCMC.ts.Rs.m89.inversion.p[[3]], lwd = 1, lty = 2, col = "cyan")
legend(12000, 0.709, c("measured ivory","modeled serum"),lwd = c(2, 1), col=c("red","cyan"))
#end plot

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,t), ylim = c(0.709, 0.714), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to days using rate Ivo.rate.mic
lines(dist.mic/Ivo.rate.mic, R.mic, lwd = 2, col = "blue") #results from micromill

MCMC.ts.Rin.m.89.p <- MCMC.ts(post.misha.inversion.p$BUGSoutput$sims.list$Rin.m)
lines(1:t,MCMC.ts.Rin.m.89.p[[1]],lwd = 2, col = "firebrick4")
lines(1:t,MCMC.ts.Rin.m.89.p[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:t,MCMC.ts.Rin.m.89.p[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))