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

########inversion model with calibration for either the one/two pool model#######
#inversion based on micromill results
micromill <- read.csv("data/Misha micromill.csv")

R.mic <- rev(micromill$X87Sr.86Sr)
R.sd.mic <- rev(micromill$StdErr)
dist.mic <- rev(micromill$dist)
n.mic <- length(dist.mic)

#calibration interval and sampling interval
cal.intv <- 100

s.intv <- 400

R0 <- 0.7071

#Re is the mean ratio of end value  
Re <- 0.7112

#parameters to save
parameters <- c("Ivo.rate", "Rs.cal", "Rb.cal", "mod.index.cal","mod.dist.cal","Rin.cal","a","b","c",
                "Rs.m", "Rb.m","Rin.m","mod.index","mod.dist", "a.m", "b.m", "c.m","Rin.m.eps.ac",
                "Ps", "Fin","Pb", "Fb", "Re.mean", "switch", "w.contrib", "h.contrib",
                "flux.ratio", "pool.ratio")

##Data to pass to the model
#compared to the turnover model that is essentially the .cal part here 
#the inversion takes the measured value of potentially a different ivory series
dat = list(R.cal = R.mea, dist.cal = dist.mea, R.sd.cal = R.sd.mea, t.cal = 800, n.cal = n.mea, 
           R0 = R0,  Re = Re, cal.intv = cal.intv, s.intv = s.intv, 
           params.mu = turnover.params.mu, params.vcov = turnover.params.vcov,
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