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
n.cal = length(dist.cal)

micromill <- read.csv("data/Misha micromill.csv")

R.mic <- rev(micromill$X87Sr.86Sr)
R.sd.mic <- rev(micromill$StdErr)
dist.mic <- rev(micromill$dist)
n.mic <- length(dist.mic)

#assign input values before and after the switch, but allows some variation
R0 <- 0.70706

#Re is the mean ratio of end value  
Re <- 0.7112

s.intv <- 400

cal.intv <- 100

#parameters to save
parameters <- c("Ivo.rate", "Rs.cal", "Rb.cal", "mod.index.cal","mod.dist.cal","Rin.cal","a","b","c",
                "Rs.m", "Rb.m","Rin.m","mod.index","mod.dist", "a.m", "b.m", "c.m","Rin.m.cps.ac")

##Data to pass to the model
#compared to the turnover model that is essentially the .cal part here 
#the inversion takes the measured value of potentially a different ivory series
dat = list(R.cal = R.cal, dist.cal = dist.cal, R.sd.cal = R.sd.cal, t.cal = 750, n.cal = n.cal, 
           R0 = R0, Re = Re, cal.intv = cal.intv, s.intv = s.intv,
           R.mea = R.mic, dist.mea = dist.mic, R.sd.mea = R.sd.mic, t = 450, n.mea = n.mic)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 4e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/600

#Run it
post.misha.inversion.woi = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS wo intake model.R", 
                                                   parameters.to.save = parameters, 
                                                   data = dat, n.chains=5, n.iter = n.iter, 
                                                   n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 8 hours

save(post.misha.inversion.woi, file = "out/post.misha.inversion.woi.RData")

post.misha.inversion.woi$BUGSoutput$summary

load("out/post.misha.inversion.woi.RData")

plot(density(post.misha.inversion.woi$BUGSoutput$sims.list$a))

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,t), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to days using rate Ivo.rate.mic
lines((max(dist.mic+200) - dist.mic)/Ivo.rate.mic + 1, R.mic, lwd = 2, col = "blue")#results from micromill

MCMC.ts.Rin.m.89.woi <- MCMC.ts(post.misha.inversion.woi$BUGSoutput$sims.list$Rin.m)
lines(1:t,MCMC.ts.Rin.m.89.woi[[1]],lwd = 2, col = "firebrick4")
lines(1:t,MCMC.ts.Rin.m.89.woi[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:t,MCMC.ts.Rin.m.89.woi[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Micromill","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))

