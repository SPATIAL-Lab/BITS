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

########inversion model with calibration for either the one/two pool model#######
#inversion based on micromill results
micromill <- read.csv("data/Misha micromill.csv")

R.mic <- rev(micromill$X87Sr.86Sr)
R.sd.mic <- rev(micromill$StdErr)
dist.mic <- rev(micromill$dist)
n.mic <- length(dist.mic)

Ivo.rate.mic <- 20

#sample interval
s.intv <- 400

#parameters to save
parameters <- c("Body.mass", "Body.mass.m","Rin.m.cps",
                "Rs.m", "Rb.m","Rin.m","dist", "a.m", "b.m", "c.m","Rin.m.cps.ac")

##Data to pass to the model
#omitting the turnover model here simplifies the model run time
#the inversion takes the measured value of potentially a different ivory series
dat = list(s.intv = s.intv, 
           params.mu = turnover.params.mu, params.vcov = turnover.params.vcov,
           R.mea = R.mic, dist.mea = dist.mic, R.sd.mea = R.sd.mic, t = 450, n.mea = n.mic)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 4e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/600

#Run it
post.misha.inv.p = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS with parameters.R", 
                                                   parameters.to.save = parameters, 
                                                   data = dat, n.chains=5, n.iter = n.iter, 
                                                   n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 8 hours

save(post.misha.inv.p, file = "out/post.misha.inv.p.RData")

post.misha.inv.p$BUGSoutput$summary

load("out/post.misha.inv.p.RData")

plot(density(post.misha.inv.p$BUGSoutput$sims.list$Rin.m.cps.ac))
plot(density(post.misha.inv.p$BUGSoutput$sims.list$Rin.m.cps))
plot(density(post.misha.inv.p$BUGSoutput$sims.list$Body.mass.m))

#plotting reconstructed Rin.m history
plot(0,0, xlim = c(1,t), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines((max(dist.mic+200) - dist.mic)/Ivo.rate.mic + 1, R.mic, lwd = 2, col = "blue") #approximate results from micromill

MCMC.inv.p.Rin.m.50 <- MCMC.CI.bound(post.misha.inv.p$BUGSoutput$sims.list$Rin.m, 0.5)
lines(1:t,MCMC.inv.p.Rin.m.50[[1]],lwd = 2, col = "firebrick4")
lines(1:t,MCMC.inv.p.Rin.m.50[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:t,MCMC.inv.p.Rin.m.50[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))

####brownian motion per step model####
micromill <- read.csv("data/Misha micromill.csv")

R.mic <- rev(micromill$X87Sr.86Sr)
R.sd.mic <- rev(micromill$StdErr)
dist.mic <- rev(micromill$dist)
n.mic <- length(dist.mic)

Ivo.rate.mic <- 20

#sample interval
s.intv <- 400

#parameters to save
parameters <- c("Body.mass", "Body.mass.m","Rin.m.cps",
                "Rs.m", "Rb.m","Rin.m","dist", "a.m", "b.m", "c.m","Rin.m.cps.ac")

##Data to pass to the model
#omitting the turnover model here simplifies the model run time
#the inversion takes the measured value of potentially a different ivory series
dat = list(s.intv = s.intv, 
           params.mu = turnover.params.mu, params.vcov = turnover.params.vcov,
           R.mea = R.mic, dist.mea = dist.mic, R.sd.mea = R.sd.mic, t = 450, n.mea = n.mic)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 4e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/600

#Run it
post.misha.inv.bm = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS BM with parameters.R", 
                                              parameters.to.save = parameters, 
                                              data = dat, n.chains=5, n.iter = n.iter, 
                                              n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 8 hours

save(post.misha.inv.bm, file = "out/post.misha.inv.bm.RData")

post.misha.inv.bm$BUGSoutput$summary

load("out/post.misha.inv.bm.RData")

#plotting modeled serum values mapped onto ivory and checking the fit of the data

plot(density(post.misha.inv.bm$BUGSoutput$sims.list$Rin.m.cps.ac))
plot(density(post.misha.inv.bm$BUGSoutput$sims.list$Rin.m.cps))
plot(density(post.misha.inv.bm$BUGSoutput$sims.list$Body.mass.m))

#plotting reconstructed Rin.m history
plot(0,0, xlim = c(1,t), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines((max(dist.mic+200) - dist.mic)/Ivo.rate.mic + 1, R.mic, lwd = 2, col = "blue") #approximate results from micromill

MCMC.inv.bm.Rin.m.89 <- MCMC.CI.bound(post.misha.inv.bm$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:t,MCMC.inv.bm.Rin.m.89[[1]],lwd = 2, col = "firebrick4")
lines(1:t,MCMC.inv.bm.Rin.m.89[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:t,MCMC.inv.bm.Rin.m.89[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))