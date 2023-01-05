library(coda)
library(rjags)
library(R2jags)
library(mcmcplots)
library(bayestestR)
library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)
library(zoo)

plot.col<-viridis(7)

setwd("C:/Users/ydmag/Google Drive/U of U/Elephant movement/Sr-in-ivory")

####supplementary figures######

#######Fig S1 Position of cracks and how they have little effect on measured Sr 87/86######

#######Fig S2 Log-normal goodness of fit for parameter a#######
#check q-q plot for fit
qqPlot(post.misha.fdnb1pr$BUGSoutput$sims.list$a[,1],distribution = "lnorm",
       param.list=list(mean=a.param.nb1pr$parameters[1],sd=a.param.nb1pr$parameters[2]),add.line=T)

#######Fig S3 results of intake model###########
#use a different JAGS file for evaluation!

#######Fig S4 sensitivity test 1 Sensitivity to excursion
#excursion makes rate parameter smaller
#########1 pool turnover model with full 50 pt average dataset###################
###prep data###
R.sd.mea <- n.sd.misha.50.sr
dist.mea <- n.avg.misha.50.dist
R.mea <- n.avg.misha.50.sr
n.mea = length(n.avg.misha.50.sr)

R0 <- 0.707

#Re is the mean ratio of end value  
Re <- 0.711
s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30
#parameters to save
parameters <- c("a", "Ivo.rate", "r1.m","Rin","dist","h.l",
                "Sr.pre", "Re.mean", "R0.mean", "switch")
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
post.misha.2p50 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS no bone woi.R", 
                                             parameters.to.save = parameters, 
                                             data = dat, n.chains=5, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 3.5 hours

post.misha.2p50$BUGSoutput$summary

save(post.misha.2p50, file = "out/post.misha.2p50.RData")
load("out/post.misha.2p50.RData")

traplot(post.misha.2p50,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.2p50,parms = c("a"))

summary(post.misha.2p50$BUGSoutput$sims.list$switch)
summary(post.misha.2p50$BUGSoutput$sims.list$Ivo.rate)

#plotting some parameters
plot(density(post.misha.2p50$BUGSoutput$sims.list$Ivo.rate))
map_estimate(post.misha.2p50$BUGSoutput$sims.list$h.l)
plot(density(post.misha.2p50$BUGSoutput$sims.list$h.l))

#plotting modeled serum values mapped onto ivory and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.2p50$BUGSoutput$sims.list$R1.m,
               post.misha.2p50$BUGSoutput$sims.list$dist)
lines(n.avg.misha.50.dist,n.avg.misha.50.sr,lwd = 2, col = "red")

#######Fig S5 sensitivity test 2 Sensitivity to switch date
#how do you determine switch date?

######Fig S6: Sensitivity test 3 Sensitivity to precision terms?


