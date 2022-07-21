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

#R0 is the mean ratio of initial value 

R0 <- 0.7071

m.feed <- 40 #mixing of x kg of hay, this is to constrain Sr variability of hay

#parameters to save
parameters <- c("a", "b","c", "Ivo.rate", "Rs.m","Rb.m","Rin","mod.index","mod.dist","f.h","f.pel","h.l.s",
                "Sr.pre", "Ps", "Fin","Pb", "Fb", "Re.mean", "switch", "w.contrib", "h.contrib",
                "flux.ratio", "pool.ratio")
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
n.iter = 5e3
n.burnin = 1e3
n.thin = floor(n.iter-n.burnin)/1000

#Run it
post.misha.5 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS v5.R", 
                                          parameters.to.save = parameters, 
                                          data = dat, n.chains=3, n.iter = n.iter, 
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

plot(density(post.misha.5$BUGSoutput$sims.list$a))

#make a contour map
contour.flux.pool <- kde2d(post.misha.5$BUGSoutput$sims.list$flux.ratio[,1], 
                           post.misha.5$BUGSoutput$sims.list$pool.ratio[,1], n = 64)
image(contour.flux.pool, col=viridis(64), xlab="Flux ratio: Fs/Fb", ylab="Pool ratio: Ps/Pb",xlim =c(0,200))
contour(contour.flux.pool,lwd = 1.5, add = TRUE, labcex = 1)

#plotting some parameters
hist(post.misha.5$BUGSoutput$sims.list$Ivo.rate)
map_estimate(post.misha.5$BUGSoutput$sims.list$Pb) #1.22 mmol
hist(post.misha.5$BUGSoutput$sims.list$Pb)
map_estimate(post.misha.5$BUGSoutput$sims.list$Ps) #1.03 mmol
map_estimate(post.misha.5$BUGSoutput$sims.list$Fin) #0.03 mmol/day
hist(post.misha.5$BUGSoutput$sims.list$Fin)
map_estimate(post.misha.5$BUGSoutput$sims.list$Fb) #0.00561 mmol/day
hist(post.misha.5$BUGSoutput$sims.list$Fb)

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

#plotting modeled serum values and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
post.misha.5.Rs.index <- MCMC.ts.dist(post.misha.5$BUGSoutput$sims.list$Rs.m, 
                                      post.misha.5$BUGSoutput$sims.list$mod.dist,
                                      post.misha.5$BUGSoutput$sims.list$mod.index)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
post.misha.5.RS.89 <- MCMC.ts(post.misha.5.Rs.index)

lines(dist.mea, post.misha.5.RS.89[[1]], lwd = 1, col = "cyan")
lines(dist.mea, post.misha.5.RS.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, post.misha.5.RS.89[[3]], lwd = 1, lty = 2, col = "cyan")

#plotting modeled bone values and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
post.misha.5.Rb.index <- MCMC.ts.dist(post.misha.5$BUGSoutput$sims.list$Rb.m, 
                                      post.misha.5$BUGSoutput$sims.list$mod.dist,
                                      post.misha.5$BUGSoutput$sims.list$mod.index)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
post.misha.5.Rb.89 <- MCMC.ts(post.misha.5.Rb.index)

lines(dist.mea, post.misha.5.Rb.89[[1]], lwd = 1, col = "cyan")
lines(dist.mea, post.misha.5.Rb.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, post.misha.5.Rb.89[[3]], lwd = 1, lty = 2, col = "cyan")

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
