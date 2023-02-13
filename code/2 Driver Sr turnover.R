library(coda)
library(rjags)
library(R2jags)
library(mcmcplots)
library(bayestestR)
library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)

source("code/1 Helper functions.R")

plot.col<-viridis(7)

setwd("C:/Users/ydmag/Google Drive/U of U/Elephant movement/Sr-in-ivory")

####The one pool model is only used in sensitivity test in the supporting information#####

misha.raw <- read.csv("data/Misha raw.csv")

misha.raw.Sr <- misha.raw$Corr..87Sr.86Sr
misha.raw.dist <- misha.raw$dist

##calculate 50 point average to reduce data amount
n.avg <- 50 #now many data points to average
n <- floor(length(misha.raw$Corr..87Sr.86Sr)/n.avg) #99 data points, discard the last <50 data points

#initiate vectors
n.avg.misha.50.sr <- rep(0, n)
n.sd.misha.50.sr <- rep(0, n)

n.avg.misha.50.dist <- rep(0, n)

for(i in 1:n){
  x <- ((i-1)*n.avg + 1):(i*n.avg)
  n.avg.misha.50.sr[i] <- mean(misha.raw.Sr[x])
  n.sd.misha.50.sr[i] <- sd(misha.raw.Sr[x])
  
  n.avg.misha.50.dist[i] <- mean(misha.raw.dist[x])#mean and median are the same
}

#check data range
plot(n.avg.misha.50.dist, n.avg.misha.50.sr, main="50 pt average",
     xlim=c(max(n.avg.misha.50.dist),min(n.avg.misha.50.dist)))
lines(n.avg.misha.50.dist, n.avg.misha.50.sr)
hist(n.sd.misha.50.sr)

##################2 pool turnover with 50 pt average (excursion rmv)################
###prep data###
ind.remv.50 <- which(n.avg.misha.50.dist > 17000 & n.avg.misha.50.dist < 17500 & n.avg.misha.50.sr < 0.7075)
ind.remv.50 #36 37 38 39
index.50.anom.remv1 <- c(1:35)
index.50.anom.remv2 <- c(40:length(n.avg.misha.50.dist))
index.50.anom.remv <- c(1:35,40:length(n.avg.misha.50.dist))
n.avg.misha.50.dist.rmv <- n.avg.misha.50.dist[index.50.anom.remv]
n.avg.misha.50.sr.rmv <- n.avg.misha.50.sr[index.50.anom.remv]
n.sd.misha.50.sr.rmv <- n.sd.misha.50.sr[index.50.anom.remv]


#######This is the version used in the article, constraint on parameter b######
R.sd.mea <- n.sd.misha.50.sr.rmv
dist.mea <- n.avg.misha.50.dist.rmv
R.mea <- n.avg.misha.50.sr.rmv
n.mea = length(n.avg.misha.50.sr.rmv)

Rpri <- 0.706

#Re is the mean ratio of end value  
Raft <- 0.711
s.intv <- 52.4

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.8
max.dist.mea <- max(n.avg.misha.50.dist) + 30
#parameters to save
parameters <- c("a", "b","c", "Ivo.rate", "R1.m","R2.m","Rin","dist.index",
                "Sr.pre", "Raft.mean", "Rpri.mean", "switch","dist",
                "flux.ratio", "pool.ratio")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 750, n.mea = n.mea, 
           Rpri = Rpri, Raft = Raft, s.intv = s.intv, Ivo.rate.mean = Ivo.rate.mean,
           Ivo.rate.sd = Ivo.rate.sd, max.dist.mea = max.dist.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 4e3
n.thin = 1

#Run it
post.misha.pc2p3 = do.call(jags.parallel,list(model.file = "code/Sr turnover JAGS wo intake model3.R", 
                                             parameters.to.save = parameters, 
                                             data = dat, n.chains=5, n.iter = n.iter, 
                                             n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 50 hours

save(post.misha.pc2p3, file = "out/post.misha.pc2p3.RData")
load("out/post.misha.pc2p3.RData")

traplot(post.misha.pc2p3,parms = c("flux.ratio", "pool.ratio"))
traplot(post.misha.pc2p3,parms = c("a", "b","c"))

plot(density(post.misha.pc2p3$BUGSoutput$sims.list$Ivo.rate))
abline(v=14.7)

#convergence params are ok
subset(post.misha.pc2p3$BUGSoutput$summary,
       rownames(post.misha.pc2p3$BUGSoutput$summary)=="a")
subset(post.misha.pc2p3$BUGSoutput$summary,
       rownames(post.misha.pc2p3$BUGSoutput$summary)=="b")
subset(post.misha.pc2p3$BUGSoutput$summary,
       rownames(post.misha.pc2p3$BUGSoutput$summary)=="c")

#save MAPEs and CIs
MAP.a <- map_estimate(post.misha.pc2p3$BUGSoutput$sims.list$a)
MAP.a[1]
log(2)/MAP.a[1]
MAP.b <- map_estimate(post.misha.pc2p3$BUGSoutput$sims.list$b)
MAP.b[1]

MAP.c <- map_estimate(post.misha.pc2p3$BUGSoutput$sims.list$c)
MAP.c[1]
log(2)/MAP.c[1]

MCMC.CI.a <- hdi(post.misha.pc2p3$BUGSoutput$sims.list$a,0.89)
MCMC.CI.a$CI_low
MCMC.CI.a$CI_high
MCMC.CI.b <- hdi(post.misha.pc2p3$BUGSoutput$sims.list$b,0.89)
MCMC.CI.b$CI_low
MCMC.CI.b$CI_high
MCMC.CI.c <- hdi(post.misha.pc2p3$BUGSoutput$sims.list$c,0.89)
MCMC.CI.c$CI_low
MCMC.CI.c$CI_high
