library(coda)
library(rjags)
library(R2jags)
library(mcmcplots)
library(bayestestR)
library(scales)
library(MASS)
library(viridisLite)

plot.col<-viridis(7)

setwd("C:/Users/ydmag/Google Drive/U of U/Elephant movement/Sr-in-ivory")
misha <- read.csv("data/Misha ivory.csv")

#helper functions

MCMC.ts <- function (MCMC.res){
  require(KernSmooth)
  require(bayestestR)
  dim.MCMC <- dim(MCMC.res)
  #typically, the first element is # of interation
  #the second element is the time series
  map.res <- rep(0, dim.MCMC[2])
  hdi.high <- rep(0, dim.MCMC[2])
  hdi.low <- rep(0, dim.MCMC[2])
  
  for(i in 1:dim.MCMC[2]){
    map.res[i] <- map_estimate(MCMC.res[,i], method = "KernSmooth")
    hdi.89 <- hdi(MCMC.res[,i], ci = 0.89)
    hdi.low[i] <- hdi.89$CI_low
    hdi.high[i] <- hdi.89$CI_high
  }
  
  return (list(map.res, hdi.low, hdi.high))
}

MCMC.ts.dist <- function(MCMC.res, MCMC.dist, MCMC.index){
  require(scales)
  
  #typically, the first element is # of interation
  #the second element is the time series
  
  dim.MCMC <- dim(MCMC.res) #this is 800
  n <- dim.MCMC[1] #number of interation
  dim.MCMC.dist <- dim(MCMC.dist) #this is 100
  MCMC.res.index <- array(0, c(n, dim.MCMC.dist[2]))
  for(i in 1:n){ #for each iteration, extract MCMC index results
    MCMC.res.index[i,] <- MCMC.res[i, MCMC.index[i,]]
  }
  
  for(i in 1:n){
    lines(MCMC.dist[i,], MCMC.res.index[i,], col = alpha("black", 0.02))
  }

  return(MCMC.res.index)
}

MCMC.dist.median <- function(MCMC.res){
  dim.MCMC <- dim(MCMC.res)
  n <- dim.MCMC[2] #number of data points
  MCMC.res.med <- rep(0, n)
  MCMC.res.max <- rep(0, n)
  MCMC.res.min <- rep(0, n)
  for(i in 1:n){ #for each iteration, extract MCMC index results
    MCMC.res.med[i] <- median (MCMC.res[,i])
    MCMC.res.min[i] <- min (MCMC.res[,i])
    MCMC.res.max[i] <- median (MCMC.res[,i])
  }
  return(MCMC.res.med)
}

R.sd.mea <- misha$sd
dist.mea <- misha$dist
R.mea <- misha$mean
n.mea = length(dist.mea)
# #R0 is the mean ratio of initial value 
# R0 <- 0.7070
# #Re is the mean ratio of end value  
# Re <- 0.711
# Rin.mean <- c(rep(R0,81),rep(Re,769))
# 
# #parameters to save
# parameters <- c("a", "b", "c", "Ivo.rate", "Rs.m", "Rb.m","Rin","mod.index","mod.dist",
#                 "flux.ratio", "pool.ratio", "Sr.pre", "Ps", "Pb", "Fin", "Fb")
# ##Data to pass to the model
# dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 850, n.mea = n.mea, 
#            Rin.mean = Rin.mean, R0 = R0, Re = Re)
# 
# #Start time
# t1 = proc.time()
# 
# set.seed(t1[3])
# n.iter = 1e4
# n.burnin = 2e3
# n.thin = floor(n.iter-n.burnin)/1000
# 
# #Run it
# post.misha = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS v1.R", 
#                                         parameters.to.save = parameters, 
#                                         data = dat, n.chains=3, n.iter = n.iter, 
#                                         n.burnin = n.burnin, n.thin = n.thin))
# 
# #Time taken
# proc.time() - t1 #~ 2 hours
# 
# post.misha$BUGSoutput$summary
# traplot(post.misha,parms = c("flux.ratio", "pool.ratio"))
# traplot(post.misha,parms = c("a", "b","c"))
# 
# hist(post.misha$BUGSoutput$sims.list$flux.ratio)
# hist(post.misha$BUGSoutput$sims.list$pool.ratio)
# #make a contour map
# contour.flux.pool <- kde2d(post.misha$BUGSoutput$sims.list$flux.ratio[,1], 
#       post.misha$BUGSoutput$sims.list$pool.ratio[,1], n = 64)
# image(contour.flux.pool, col=viridis(64), xlab="Flux ratio: Fs/Fb", ylab="Pool ratio: Ps/Pb",xlim =c(0,200))
# contour(contour.flux.pool,lwd = 1.5, add = TRUE, labcex = 1)
# 
# hist(post.misha$BUGSoutput$sims.list$Ivo.rate)
# map_estimate(post.misha$BUGSoutput$sims.list$Pb) #3.24 mmol
# hist(post.misha$BUGSoutput$sims.list$Pb)
# map_estimate(post.misha$BUGSoutput$sims.list$Ps) #1.05 mmol
# map_estimate(post.misha$BUGSoutput$sims.list$Fin) #0.03 mmol/day
# hist(post.misha$BUGSoutput$sims.list$Fin)
# map_estimate(post.misha$BUGSoutput$sims.list$Fb) #0.00736 mmol/day
# hist(post.misha$BUGSoutput$sims.list$Fb)
# 
# 
# plot(density(post.misha$BUGSoutput$sims.list$Ps), type = "l", lwd = 2, xlim = c(0, 12),
#      ylim = c(0, 4), col = plot.col[2], xlab = "Pool size (mmol)", main = "")
# lines(density(post.misha$BUGSoutput$sims.list$Pb), lwd = 2, col = plot.col[6])
# legend(8, 3, c("Serum","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])
# 
# plot(density(post.misha$BUGSoutput$sims.list$Fin), type = "l", lwd = 2, xlim = c(0, 0.1),
#      ylim = c(0, 80), col = plot.col[2], xlab = "Sr Flux (mmol/day)", main = "")
# lines(density(post.misha$BUGSoutput$sims.list$Fb), lwd = 2, col = plot.col[6])
# legend(0.06, 80, c("Intake","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])
# 
# #a/b = Fin/Fb
# #c/b = Ps/Pb
# 
# plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
# abline(h = R0, lwd = 2, lty = 2)
# abline(h = Re, lwd = 2, lty = 2)
# MCMC.ts.Rs.index <- MCMC.ts.dist(post.misha$BUGSoutput$sims.list$Rs.m, 
#                                  post.misha$BUGSoutput$sims.list$mod.dist,
#                                 post.misha$BUGSoutput$sims.list$mod.index)
# lines(misha$dist,misha$mean,lwd = 2, col = "red")
# MCMC.ts.RS.89 <- MCMC.ts(MCMC.ts.Rs.index)
# #Misha.med.dist <- MCMC.dist.median(post.misha$BUGSoutput$sims.list$mod.dist)
# 
# lines(dist.mea, MCMC.ts.RS.89[[1]], lwd = 1, col = "cyan")
# lines(dist.mea, MCMC.ts.RS.89[[2]], lwd = 1, lty = 2, col = "cyan")
# lines(dist.mea, MCMC.ts.RS.89[[3]], lwd = 1, lty = 2, col = "cyan")
# 
# plot(density(post.misha$BUGSoutput$sims.list$Ps), type = "l", lwd = 2, xlim = c(0, 12),
#      ylim = c(0, 3), col = plot.col[2], xlab = "Pool size (mmol)", main = "")
# lines(density(post.misha$BUGSoutput$sims.list$Pb), lwd = 2, col = plot.col[6])
# legend(8, 3, c("Serum","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])
# 
# plot(density(post.misha$BUGSoutput$sims.list$Fin), type = "l", lwd = 2, xlim = c(0, 0.4),
#      ylim = c(0, 20), col = plot.col[2], xlab = "Sr Flux (mmol/day)", main = "")
# lines(density(post.misha$BUGSoutput$sims.list$Fb), lwd = 2, col = plot.col[6])
# legend(0.3, 20, c("Intake","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

#####testing version 2 of the turnover model with evaluation on the mean of neighbouring data points##
R.sd.mea <- misha$sd
dist.mea <- misha$dist
R.mea <- misha$mean
n.mea = length(dist.mea)
#R0 is the mean ratio of initial value 
R0 <- 0.7070
#Re is the mean ratio of end value  
Re <- 0.711
Rin.mean <- c(rep(R0,81),rep(Re,769))

hist(1/rgamma(1000, shape = 10, rate = 0.2)^0.5)

#parameters to save
parameters <- c("a", "b", "c", "Ivo.rate", "Rs.m", "Rb.m","Rin","mod.index","mod.dist",
                "flux.ratio", "pool.ratio", "Sr.pre", "Ps", "Pb", "Fin", "Fb")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 850, n.mea = n.mea, 
           Rin.mean = Rin.mean, R0 = R0, Re = Re)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 2e3
n.thin = floor(n.iter-n.burnin)/1000

#Run it
post.misha = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS v2.R", 
                                        parameters.to.save = parameters, 
                                        data = dat, n.chains=3, n.iter = n.iter, 
                                        n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 2 hours

post.misha$BUGSoutput$summary
traplot(post.misha,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha,parms = c("a", "b","c"))

hist(post.misha$BUGSoutput$sims.list$flux.ratio)
hist(post.misha$BUGSoutput$sims.list$pool.ratio)
#make a contour map
contour.flux.pool <- kde2d(post.misha$BUGSoutput$sims.list$flux.ratio[,1], 
                           post.misha$BUGSoutput$sims.list$pool.ratio[,1], n = 64)
image(contour.flux.pool, col=viridis(64), xlab="Flux ratio: Fs/Fb", ylab="Pool ratio: Ps/Pb",xlim =c(0,200))
contour(contour.flux.pool,lwd = 1.5, add = TRUE, labcex = 1)

#plotting some parameters
hist(post.misha$BUGSoutput$sims.list$Ivo.rate)
map_estimate(post.misha$BUGSoutput$sims.list$Pb) #3.24 mmol
hist(post.misha$BUGSoutput$sims.list$Pb)
map_estimate(post.misha$BUGSoutput$sims.list$Ps) #1.05 mmol
map_estimate(post.misha$BUGSoutput$sims.list$Fin) #0.03 mmol/day
hist(post.misha$BUGSoutput$sims.list$Fin)
map_estimate(post.misha$BUGSoutput$sims.list$Fb) #0.00736 mmol/day
hist(post.misha$BUGSoutput$sims.list$Fb)

plot(density(post.misha$BUGSoutput$sims.list$Ps), type = "l", lwd = 2, xlim = c(0, 12),
     ylim = c(0, 4), col = plot.col[2], xlab = "Pool size (mmol)", main = "")
lines(density(post.misha$BUGSoutput$sims.list$Pb), lwd = 2, col = plot.col[6])
legend(8, 3, c("Serum","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

plot(density(post.misha$BUGSoutput$sims.list$Fin), type = "l", lwd = 2, xlim = c(0, 0.1),
     ylim = c(0, 80), col = plot.col[2], xlab = "Sr Flux (mmol/day)", main = "")
lines(density(post.misha$BUGSoutput$sims.list$Fb), lwd = 2, col = plot.col[6])
legend(0.06, 80, c("Intake","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

#a/b = Fin/Fb
#c/b = Ps/Pb

#plotting modeled serum values and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rs.index <- MCMC.ts.dist(post.misha$BUGSoutput$sims.list$Rs.m, 
                                 post.misha$BUGSoutput$sims.list$mod.dist,
                                 post.misha$BUGSoutput$sims.list$mod.index)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
MCMC.ts.RS.89 <- MCMC.ts(MCMC.ts.Rs.index)
#Misha.med.dist <- MCMC.dist.median(post.misha$BUGSoutput$sims.list$mod.dist)

lines(dist.mea, MCMC.ts.RS.89[[1]], lwd = 1, col = "cyan")
lines(dist.mea, MCMC.ts.RS.89[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, MCMC.ts.RS.89[[3]], lwd = 1, lty = 2, col = "cyan")

plot(density(post.misha$BUGSoutput$sims.list$Ps), type = "l", lwd = 2, xlim = c(0, 12),
     ylim = c(0, 3), col = plot.col[2], xlab = "Pool size (mmol)", main = "")
lines(density(post.misha$BUGSoutput$sims.list$Pb), lwd = 2, col = plot.col[6])
legend(8, 3, c("Serum","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

plot(density(post.misha$BUGSoutput$sims.list$Fin), type = "l", lwd = 2, xlim = c(0, 0.4),
     ylim = c(0, 20), col = plot.col[2], xlab = "Sr Flux (mmol/day)", main = "")
lines(density(post.misha$BUGSoutput$sims.list$Fb), lwd = 2, col = plot.col[6])
legend(0.3, 20, c("Intake","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

####v2 model with 0.7119 end value#####
#R0 is the mean ratio of initial value 
R0 <- 0.70706

#Re is the mean ratio of end value  
Re <- 0.71179

Rin.mean <- c(rep(R0,81),rep(Re,819))

#parameters to save
parameters <- c("a", "b", "c", "Ivo.rate", "Rs.m", "Rb.m","Rin","mod.index","mod.dist",
                "flux.ratio", "pool.ratio", "Sr.pre", "Ps", "Pb", "Fin", "Fb")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 900, n.mea = n.mea, 
           Rin.mean = Rin.mean, R0 = R0, Re = Re)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 2e3
n.thin = floor(n.iter-n.burnin)/1000

#Run it
post.misha.19 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS v2.R", 
                                        parameters.to.save = parameters, 
                                        data = dat, n.chains=3, n.iter = n.iter, 
                                        n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 2 hours
save(post.misha.19, file = "out/post_misha_19.RData")

post.misha.19$BUGSoutput$summary
traplot(post.misha.19,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.19,parms = c("a", "b","c"))

hist(post.misha.19$BUGSoutput$sims.list$flux.ratio)
hist(post.misha.19$BUGSoutput$sims.list$pool.ratio)
plot(density(post.misha.19$BUGSoutput$sims.list$c))
plot(density(post.misha.19$BUGSoutput$sims.list$a))
plot(density(post.misha.19$BUGSoutput$sims.list$b))

map_estimate(post.misha.19$BUGSoutput$sims.list$Pb) #3.14 mmol
map_estimate(post.misha.19$BUGSoutput$sims.list$Ps) #1.02 mmol
map_estimate(post.misha.19$BUGSoutput$sims.list$Fin) #0.03 mmol/day
map_estimate(post.misha.19$BUGSoutput$sims.list$Fb) #0.007 mmol/day
map_estimate(post.misha.19$BUGSoutput$sims.list$a) #0.03
map_estimate(post.misha.19$BUGSoutput$sims.list$c) #0.003

#make a contour map
contour.flux.pool.19 <- kde2d(post.misha.19$BUGSoutput$sims.list$flux.ratio[,1], 
                           post.misha.19$BUGSoutput$sims.list$pool.ratio[,1], n = 64)
image(contour.flux.pool.19, col=viridis(64), xlab="Flux ratio: Fin/Fb", ylab="Pool ratio: Ps/Pb")
contour(contour.flux.pool.19,lwd = 1.5, add = TRUE, labcex = 1)

####R serum ratio###
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rs.index.19 <- MCMC.ts.dist(post.misha.19$BUGSoutput$sims.list$Rs.m, 
                                 post.misha.19$BUGSoutput$sims.list$mod.dist,
                                 post.misha.19$BUGSoutput$sims.list$mod.index)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
MCMC.ts.RS.89.19 <- MCMC.ts(MCMC.ts.Rs.index.19)
#Misha.med.dist <- MCMC.dist.median(post.misha$BUGSoutput$sims.list$mod.dist)

lines(dist.mea, MCMC.ts.RS.89.19[[1]], lwd = 1, col = "cyan")
lines(dist.mea, MCMC.ts.RS.89.19[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, MCMC.ts.RS.89.19[[3]], lwd = 1, lty = 2, col = "cyan")

#R bone ratio##
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rb.index.19 <- MCMC.ts.dist(post.misha.19$BUGSoutput$sims.list$Rb.m, 
                                    post.misha.19$BUGSoutput$sims.list$mod.dist,
                                    post.misha.19$BUGSoutput$sims.list$mod.index)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
MCMC.ts.Rb.89.19 <- MCMC.ts(MCMC.ts.Rb.index.19)
#Misha.med.dist <- MCMC.dist.median(post.misha$BUGSoutput$sims.list$mod.dist)

lines(dist.mea, MCMC.ts.Rb.89.19[[1]], lwd = 1, col = "cyan")
lines(dist.mea, MCMC.ts.Rb.89.19[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, MCMC.ts.Rb.89.19[[3]], lwd = 1, lty = 2, col = "cyan")


plot(density(post.misha.19$BUGSoutput$sims.list$Ps), type = "l", lwd = 2, xlim = c(0, 12),
     ylim = c(0, 3), col = plot.col[2], xlab = "Pool size (mmol)", main = "")
lines(density(post.misha.19$BUGSoutput$sims.list$Pb), lwd = 2, col = plot.col[6])
legend(8, 3, c("Serum","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

plot(density(post.misha.19$BUGSoutput$sims.list$Fin), type = "l", lwd = 2, xlim = c(0, 0.1),
     ylim = c(0, 80), col = plot.col[2], xlab = "Sr Flux (mmol/day)", main = "")
lines(density(post.misha.19$BUGSoutput$sims.list$Fb), lwd = 2, col = plot.col[6])
legend(0.3, 20, c("Intake","Bone"),lwd = c(2, 2), col = plot.col[c(2, 6)])

#####inversion for Rin using posterior distributions of a, b and c####
# R0 <- 0.70706
# 
# #Re is the mean ratio of end value  
# Re <- 0.71179
# 
# Rin.mean <- c(rep(R0,75),rep(Re,1125))
# #parameters to save
# parameters <- c("Ivo.rate", "Rs.cal", "Rb.cal", "mod.index.cal","mod.dist.cal","Rin.cal","a","b","c",
#                 "Rs.m", "Rb.m","Rin.m","mod.index","mod.dist", "a.m", "b.m", "c.m","Rin.m.eps.ac")
# ##Data to pass to the model
# dat = list(R.cal = R.mea, dist.cal = dist.mea, R.sd.cal = R.sd.mea, t.cal = 1200, n.cal = n.mea, 
#            Rin.mean = Rin.mean, R0 = R0, Re = Re, 
#            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 1200, n.mea = n.mea)
# 
# #Start time
# t1 = proc.time()
# 
# set.seed(t1[3])
# n.iter = 1e4
# n.burnin = 2e3
# n.thin = floor(n.iter-n.burnin)/600
# 
# #Run it
# post.misha.inversion = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS v1.R", 
#                                            parameters.to.save = parameters, 
#                                            data = dat, n.chains=5, n.iter = n.iter, 
#                                            n.burnin = n.burnin, n.thin = n.thin))
# 
# #Time taken
# proc.time() - t1 #~ 4 hours
# 
# #save(post.misha.inversion, file = "out/post.misha.inversion.RData")
# 
# post.misha.inversion$BUGSoutput$summary
# 
# plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
# abline(h = R0, lwd = 2, lty = 2)
# abline(h = Re, lwd = 2, lty = 2)
# MCMC.ts.Rs.cal.index.inversion <- MCMC.ts.dist(post.misha.inversion$BUGSoutput$sims.list$Rs.cal, 
#                                         post.misha.inversion$BUGSoutput$sims.list$mod.dist.cal,
#                                         post.misha.inversion$BUGSoutput$sims.list$mod.index.cal)
# lines(misha$dist,misha$mean,lwd = 2, col = "red")
# MCMC.ts.Rs.89.inversion <- MCMC.ts(MCMC.ts.Rs.cal.index.inversion)
# 
# lines(dist.mea, MCMC.ts.Rs.89.inversion[[1]], lwd = 1, col = "cyan")
# lines(dist.mea, MCMC.ts.Rs.89.inversion[[2]], lwd = 1, lty = 2, col = "cyan")
# lines(dist.mea, MCMC.ts.Rs.89.inversion[[3]], lwd = 1, lty = 2, col = "cyan")
# 
# #reconstructed Rin
# plot(0,0, xlim = c(1,1000), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
# abline(h = R0, lwd = 2, lty = 2)
# abline(h = Re, lwd = 2, lty = 2)
# MCMC.ts.Rin.cal.89 <- MCMC.ts(post.misha.inversion$BUGSoutput$sims.list$Rin.cal)
# lines(1:1200,MCMC.ts.Rin.cal.89[[1]],lwd = 2, col = "black")
# lines(1:1200,MCMC.ts.Rin.cal.89[[2]], lwd = 1, lty = 2, col = "black")
# lines(1:1200,MCMC.ts.Rin.cal.89[[3]], lwd = 1, lty = 2, col = "black")
# 
# MCMC.ts.Rin.m.89 <- MCMC.ts(post.misha.inversion$BUGSoutput$sims.list$Rin.m)
# lines(1:1200,MCMC.ts.Rin.m.89[[1]],lwd = 2, col = "red")
# lines(1:1200,MCMC.ts.Rin.m.89[[2]], lwd = 1, lty = 2, col = "red")
# lines(1:1200,MCMC.ts.Rin.m.89[[3]], lwd = 1, lty = 2, col = "red")
# 
# plot(density(post.misha.inversion$BUGSoutput$sims.list$Rin.m.eps.ac),type="l")
# plot(density(post.misha.inversion$BUGSoutput$sims.list$Ivo.rate),type="l")
# plot(density(post.misha.inversion$BUGSoutput$sims.list$a.m),type="l")
# plot(density(post.misha.inversion$BUGSoutput$sims.list$b.m),type="l")
# plot(density(post.misha.inversion$BUGSoutput$sims.list$c.m),type="l")

#####v2 inversion for Rin using posterior distributions of a, b and c####
#assign input values before and after the switch
R0 <- 0.70706

#Re is the mean ratio of end value  
Re <- 0.71179

#assign the mean values of the input sequence. The switch is at a fixed location
Rin.mean <- c(rep(R0,75),rep(Re,1125))

#parameters to save
parameters <- c("Ivo.rate", "Rs.cal", "Rb.cal", "mod.index.cal","mod.dist.cal","Rin.cal","a","b","c",
                "Rs.m", "Rb.m","Rin.m","mod.index","mod.dist", "a.m", "b.m", "c.m","Rin.m.eps.ac")

##Data to pass to the model
#compared to the turnover model that is essentially the .cal part here 
#the inversion takes the measured value of potentially a different ivory series
dat = list(R.cal = R.mea, dist.cal = dist.mea, R.sd.cal = R.sd.mea, t.cal = 1200, n.cal = n.mea, 
           Rin.mean = Rin.mean, R0 = R0, Re = Re, 
           R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 1200, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 2e3
n.thin = floor(n.iter-n.burnin)/600

#Run it
post.misha.inversion2 = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS v2.R", 
                                                  parameters.to.save = parameters, 
                                                  data = dat, n.chains=5, n.iter = n.iter, 
                                                  n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 8 hours

save(post.misha.inversion2, file = "out/post.misha.inversion2.RData")

post.misha.inversion2$BUGSoutput$summary

#plotting modeled serum values and checking the fit of the data
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rs.cal.index.inversion2 <- MCMC.ts.dist(post.misha.inversion2$BUGSoutput$sims.list$Rs.cal, 
                                               post.misha.inversion2$BUGSoutput$sims.list$mod.dist.cal,
                                               post.misha.inversion2$BUGSoutput$sims.list$mod.index.cal)
lines(misha$dist,misha$mean,lwd = 2, col = "red")
MCMC.ts.Rs.89.inversion2 <- MCMC.ts(MCMC.ts.Rs.cal.index.inversion2)

lines(dist.mea, MCMC.ts.Rs.89.inversion2[[1]], lwd = 1, col = "cyan")
lines(dist.mea, MCMC.ts.Rs.89.inversion2[[2]], lwd = 1, lty = 2, col = "cyan")
lines(dist.mea, MCMC.ts.Rs.89.inversion2[[3]], lwd = 1, lty = 2, col = "cyan")
legend(12000, 0.709, c("measured ivory","modeled serum"),lwd = c(2, 1), col=c("red","cyan"))
#end plot

#plotting reconstructed Rin and compared that with posterior of Rin
plot(0,0, xlim = c(1,1000), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
MCMC.ts.Rin.cal.89.2 <- MCMC.ts(post.misha.inversion2$BUGSoutput$sims.list$Rin.cal)
lines(1:1200,MCMC.ts.Rin.cal.89.2[[1]],lwd = 2, col = "black")
lines(1:1200,MCMC.ts.Rin.cal.89.2[[2]], lwd = 1, lty = 2, col = "black")
lines(1:1200,MCMC.ts.Rin.cal.89.2[[3]], lwd = 1, lty = 2, col = "black")

MCMC.ts.Rin.m.89.2 <- MCMC.ts(post.misha.inversion2$BUGSoutput$sims.list$Rin.m)
lines(1:1200,MCMC.ts.Rin.m.89.2[[1]],lwd = 2, col = "red")
lines(1:1200,MCMC.ts.Rin.m.89.2[[2]], lwd = 1, lty = 2, col = "red")
lines(1:1200,MCMC.ts.Rin.m.89.2[[3]], lwd = 1, lty = 2, col = "red")
legend(0, 0.716, c("PD input","PD reconstructed"),lwd = c(2, 2), col=c("black","red"))
#end plot

