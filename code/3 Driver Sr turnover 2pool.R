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

##########two pool turnover model########
###with intake model component###

#R0 is the mean ratio of initial value 

R0 <- 0.7071

m.feed <- 40 #mixing of x kg of hay, this is to constrain Sr variability of hay

s.intv <- 100

#parameters to save
parameters <- c("a", "b","c", "Ivo.rate", "Rs.m","Rb.m","Rin","mod.index","mod.dist","f.h","f.pel","h.l.s",
                "Sr.pre", "Ps", "Fin","Pb", "Fb", "Re.mean", "switch", "w.contrib", "h.contrib",
                "flux.ratio", "pool.ratio")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 800, n.mea = n.mea, 
           R0 = R0, Sr.hay.mean = Sr.hay.mean, Sr.hay.sd = Sr.hay.sd, 
           Sr.pel.mean = Sr.pel.mean, Sr.pel.sd = Sr.pel.sd, 
           Sr.sup.mean = Sr.sup.mean, Sr.sup.sd = Sr.sup.sd,
           Sr.w.mean = Sr.w.mean, Sr.w.sd = Sr.w.sd, m.feed = m.feed,
           conc.hay.mean = conc.hay.mean, conc.hay.sd = conc.hay.sd, 
           conc.pel.mean = conc.pel.mean, conc.pel.sd = conc.pel.sd,
           conc.sup.mean = conc.sup.mean, conc.sup.sd = conc.sup.sd,
           conc.w.mean = conc.w.mean, conc.w.sd = conc.w.sd, s.intv = s.intv)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e4
n.burnin = 2e3
n.thin = floor(n.iter-n.burnin)/600

#Run it
post.misha.5 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS v6.R", 
                                          parameters.to.save = parameters, 
                                          data = dat, n.chains=5, n.iter = n.iter, 
                                          n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 5 hours

post.misha.5$BUGSoutput$summary

save(post.misha.5, file = "out/post.misha.5.RData")
load("out/post.misha.5.RData")

traplot(post.misha.5,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.5,parms = c("a", "b","c"))

summary(post.misha.5$BUGSoutput$sims.list$switch)
summary(post.misha.5$BUGSoutput$sims.list$Ivo.rate)

#plotting some parameters
plot(density(post.misha.5$BUGSoutput$sims.list$Ivo.rate))
map_estimate(post.misha.5$BUGSoutput$sims.list$Pb) #1.22 mmol
plot(density(post.misha.5$BUGSoutput$sims.list$Pb))
map_estimate(post.misha.5$BUGSoutput$sims.list$Ps) #1.03 mmol
map_estimate(post.misha.5$BUGSoutput$sims.list$Fin) #0.03 mmol/day
plot(density(post.misha.5$BUGSoutput$sims.list$Fin))
map_estimate(post.misha.5$BUGSoutput$sims.list$Fb) #0.00561 mmol/day
plot(density(post.misha.5$BUGSoutput$sims.list$Fb))

plot(density(post.misha.5$BUGSoutput$sims.list$Ps), type = "l", lwd = 2, xlim = c(0, 12),
     ylim = c(0, 4), col = plot.col[2], xlab = "Pool size (mmol)", main = "")
lines(density(post.misha.5$BUGSoutput$sims.list$Pb), lwd = 2, col = plot.col[6])
legend(8, 3, c("Serum","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

plot(density(post.misha.5$BUGSoutput$sims.list$Fin), type = "l", lwd = 2, xlim = c(0, 0.1),
     ylim = c(0, 80), col = plot.col[2], xlab = "Sr Flux (mmol/day)", main = "")
lines(density(post.misha.5$BUGSoutput$sims.list$Fb), lwd = 2, col = plot.col[6])
legend(0.06, 80, c("Intake","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

#a/b = Fin/Fb
#c/b = Ps/Pb

#plotting modeled serum values mapped onto ivory and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.5$BUGSoutput$sims.list$Rs.m,
               post.misha.5$BUGSoutput$sims.list$dist)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
post.misha.5.Rs.m.89<- MCMC.CI.bound(post.misha.5$BUGSoutput$sims.list$Rs.m, 0.89)
#extract median distance from MCMC(t) results
med.dist.misha.5<- MCMC.dist.median(post.misha.5$BUGSoutput$sims.list$dist)
lines(med.dist.misha.5, post.misha.5.Rs.m.89[[1]], lwd = 1, col = "cyan")
lines(med.dist.misha.5, post.misha.5.Rs.m.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(med.dist.misha.5, post.misha.5.Rs.m.89[[3]], lwd = 1, lty = 2, col = "cyan")

####check the Rin values####
plot(0,0, xlim = c(0,t), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.tl.plot(post.misha.5$BUGSoutput$sims.list$Rin,t)
post.misha.5.Rin.89<- MCMC.CI.bound(post.misha.5$BUGSoutput$sims.list$Rin, 0.89)
lines(1:t, post.misha.5.Rin.89[[1]], lwd = 1, col = "cyan")
lines(1:t, post.misha.5.Rin.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(1:t, post.misha.5.Rin.89[[3]], lwd = 1, lty = 2, col = "cyan")

#reconstructed Rs.m history
plot(0,0, xlim = c(0,t), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.tl.plot(post.misha.woint3$BUGSoutput$sims.list$Rs.m,t)
lines(1:t, post.misha.5.Rs.m.89[[1]], lwd = 1, col = "cyan")
lines(1:t, post.misha.5.Rs.m.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(1:t, post.misha.5.Rs.m.89[[3]], lwd = 1, lty = 2, col = "cyan")

plot(density(post.misha.4$BUGSoutput$sims.list$Ps), type = "l", lwd = 2, xlim = c(0, 12),
     ylim = c(0, 3), col = plot.col[2], xlab = "Pool size (mmol)", main = "")
lines(density(post.misha.4$BUGSoutput$sims.list$Pb), lwd = 2, col = plot.col[6])
legend(8, 3, c("Serum","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

plot(density(post.misha.4$BUGSoutput$sims.list$Fin), type = "l", lwd = 2, xlim = c(0, 0.4),
     ylim = c(0, 20), col = plot.col[2], xlab = "Sr Flux (mmol/day)", main = "")
lines(density(post.misha.4$BUGSoutput$sims.list$Fb), lwd = 2, col = plot.col[6])
legend(0.3, 20, c("Intake","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

#####without intake model try randomly generated distance sequence###
#this also works with variable growth rates#

R.sd.mea <- misha$sd
dist.mea <- misha$dist
R.mea <- misha$mean
n.mea = length(dist.mea)
#R0 is the mean ratio of initial value 
R0 <- 0.7071

#Re is the mean ratio of end value  
Re <- 0.7112
s.intv <- 100

#parameters to save
parameters <- c("a", "b","c", "Ivo.rate", "Rs.m","Rb.m","Rin","dist.index",
                "Sr.pre", "Ps", "Fin","Pb", "Fb", "Re.mean", "switch","dist",
                "flux.ratio", "pool.ratio","Body.mass")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 750, n.mea = n.mea, 
           R0 = R0, Re = Re, s.intv = s.intv)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e3
n.burnin = 2e2
n.thin = floor(n.iter-n.burnin)/600

#Run it
post.misha.woint3 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS wo intake model3.R", 
                                              parameters.to.save = parameters, 
                                              data = dat, n.chains=5, n.iter = n.iter, 
                                              n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 20 hours

post.misha.woint3$BUGSoutput$summary

save(post.misha.woint3, file = "out/post.misha.woint3.RData")
load("out/post.misha.woint3.RData")

traplot(post.misha.woint3,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.woint3,parms = c("a", "b","c"))

summary(post.misha.woint3$BUGSoutput$sims.list$switch)
summary(post.misha.woint3$BUGSoutput$sims.list$Ivo.rate)

#make a contour map
# contour.flux.pool <- kde2d(post.misha.5$BUGSoutput$sims.list$flux.ratio[,1], 
#                            post.misha.5$BUGSoutput$sims.list$pool.ratio[,1], n = 64)
# image(contour.flux.pool, col=viridis(64), xlab="Flux ratio: Fs/Fb", ylab="Pool ratio: Ps/Pb",xlim =c(0,200))
# contour(contour.flux.pool,lwd = 1.5, add = TRUE, labcex = 1)

#plotting some parameters
plot(density(post.misha.woint3$BUGSoutput$sims.list$Ivo.rate))
map_estimate(post.misha.woint3$BUGSoutput$sims.list$Pb) #1.27 mmol
plot(density(post.misha.woint3$BUGSoutput$sims.list$Pb))
map_estimate(post.misha.woint3$BUGSoutput$sims.list$Ps) #1.05 mmol
map_estimate(post.misha.woint3$BUGSoutput$sims.list$Fin) #0.02 mmol/day
plot(density(post.misha.woint3$BUGSoutput$sims.list$Fin))
map_estimate(post.misha.woint3$BUGSoutput$sims.list$Fb) #0.00520 mmol/day
plot(density(post.misha.woint3$BUGSoutput$sims.list$Fb))

plot(density(post.misha.woint3$BUGSoutput$sims.list$Ps), type = "l", lwd = 2, xlim = c(0, 12),
     ylim = c(0, 4), col = plot.col[2], xlab = "Pool size (mmol)", main = "")
lines(density(post.misha.woint$BUGSoutput3$sims.list$Pb), lwd = 2, col = plot.col[6])
legend(8, 3, c("Serum","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

plot(density(post.misha.woint3$BUGSoutput$sims.list$Fin), type = "l", lwd = 2, xlim = c(0, 0.1),
     ylim = c(0, 80), col = plot.col[2], xlab = "Sr Flux (mmol/day)", main = "")
lines(density(post.misha.woint$BUGSoutput3$sims.list$Fb), lwd = 2, col = plot.col[6])
legend(0.06, 80, c("Intake","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

#a/b = Fin/Fb
#c/b = Ps/Pb

#plotting modeled serum values mapped onto ivory and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.woint3$BUGSoutput$sims.list$Rs.m,
               post.misha.woint3$BUGSoutput$sims.list$dist)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
post.misha.woint3.Rs.m.89<- MCMC.CI.bound(post.misha.woint3$BUGSoutput$sims.list$Rs.m, 0.89)
#extract median distance from MCMC(t) results
med.dist.woint3<- MCMC.dist.median(post.misha.woint3$BUGSoutput$sims.list$dist)
lines(med.dist.woint3, post.misha.woint3.Rs.m.89[[1]], lwd = 1, col = "cyan")
lines(med.dist.woint3, post.misha.woint3.Rs.m.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(med.dist.woint3, post.misha.woint3.Rs.m.89[[3]], lwd = 1, lty = 2, col = "cyan")

####check the Rin values####
plot(0,0, xlim = c(0,t), ylim = c(0.706, 0.713), xlab = "days", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
post.misha.woint3.Rin.89<- MCMC.CI.bound(post.misha.woint3$BUGSoutput$sims.list$Rin, 0.89)
lines(1:t, post.misha.woint3.Rin.89[[1]], lwd = 2, col = "firebrick4")
lines(1:t, post.misha.woint3.Rin.89[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:t, post.misha.woint3.Rin.89[[3]], lwd = 1, lty = 2, col = "firebrick4")

#reconstructed Rs.m history
plot(0,0, xlim = c(0,t), ylim = c(0.706, 0.713), xlab = "days", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.tl.plot(post.misha.woint3$BUGSoutput$sims.list$Rs.m,750)
lines(1:750, post.misha.woint3.Rs.m.89[[1]], lwd = 1, col = "cyan")
lines(1:750, post.misha.woint3.Rs.m.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(1:750, post.misha.woint3.Rs.m.89[[3]], lwd = 1, lty = 2, col = "cyan")

#check parameters a, b, and c
par(mfrow = c(1,3))
plot(density(post.misha.woint3$BUGSoutput$sims.list$a))
plot(density(post.misha.woint3$BUGSoutput$sims.list$b))
plot(density(post.misha.woint3$BUGSoutput$sims.list$c))

###MAP estimates, and 89% CI for parameters a, b, and c
MCMC.CI.a <- MCMC.CI.bound(post.misha.woint3$BUGSoutput$sims.list$a, 0.89)
MCMC.CI.b <- MCMC.CI.bound(post.misha.woint3$BUGSoutput$sims.list$b, 0.89)
MCMC.CI.c <- MCMC.CI.bound(post.misha.woint3$BUGSoutput$sims.list$c, 0.89)

#estimate log-normal parameters for each, but they are correlated! should use a correlated structure!
a.param <- elnorm(post.misha.woint3$BUGSoutput$sims.list$a[,1])
b.param <- elnorm(post.misha.woint3$BUGSoutput$sims.list$b[,1])
c.param <- elnorm(post.misha.woint3$BUGSoutput$sims.list$c[,1])

#log transform posterior distribution
log.a <- log(post.misha.woint3$BUGSoutput$sims.list$a[,1])
log.b <- log(post.misha.woint3$BUGSoutput$sims.list$b[,1])
log.c <- log(post.misha.woint3$BUGSoutput$sims.list$c[,1])
turnover.params.mu <- c(mean(log.a), mean(log.b), mean(log.c))

turnover.params<- data.frame(log.a, log.b, log.c)
#parameters are correlated, calculate v.cov matrix
turnover.params.vcov <- var(turnover.params)

#check the log-normal fit using q-q plots
qqPlot(post.misha.woint3$BUGSoutput$sims.list$a[,1],distribution = "lnorm",
       param.list=list(mean=a.param$parameters[1],sd=a.param$parameters[2]),add.line=T)

qqPlot(post.misha.woint3$BUGSoutput$sims.list$b[,1],distribution = "lnorm",
       param.list=list(mean=b.param$parameters[1],sd=b.param$parameters[2]),add.line=T)

qqPlot(post.misha.woint3$BUGSoutput$sims.list$c[,1],distribution = "lnorm",
       param.list=list(mean=c.param$parameters[1],sd=c.param$parameters[2]),add.line=T)

#plot prior distributions
abc.prior.params <- pri.multi.norm.den(-10,0,turnover.params.mu,turnover.params.vcov)
par(mfrow=c(1,3))
plot(abc.prior.params$x,abc.prior.params$y[1,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "a", ylab= "density")
plot(abc.prior.params$x,abc.prior.params$y[2,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "b", ylab= "density")
plot(abc.prior.params$x,abc.prior.params$y[3,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "c", ylab= "density")