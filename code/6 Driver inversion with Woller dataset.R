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

####invterting laser abalition data from Wooller et al., 2021####
Wooller <- read.csv("data/Wooller_Data_S3.csv")

#dist is in cm, convert to mm
wooller.micron <- Wooller$Dist_Seg01*10000

#foward model simulating micromill results (500 micron band)
mm.bwidth <- 500 #microns
index.wooller.dist<- ceiling(wooller.micron/mm.bwidth) #this is about the same as averaging per 100 data points!

#number of micromill simulations, discarding the last bit of data in the sequence
mm.sim.n <- max(index.wooller.dist) - 1 

mm.sim.avg.Wooller.sr <- rep(0,mm.sim.n)#initiate vectors
mm.sim.sd.Wooller.sr <- rep(0,mm.sim.n)
mm.sim.avg.Wooller.dist <- rep(0,mm.sim.n)

for(i in 1:mm.sim.n){
  temp.sr <- subset(Wooller.sr,index.wooller.dist==i)
  temp.dist <- subset(wooller.micron,index.wooller.dist==i)

  mm.sim.avg.Wooller.sr[i] <- mean(temp.sr)
  mm.sim.sd.Wooller.sr[i] <- sd(temp.sr)
  
  mm.sim.avg.Wooller.dist[i] <- median(temp.dist)#median is less sensitive to potential data gaps
}#takes ~30s 

plot(mm.sim.avg.Wooller.dist, mm.sim.avg.Wooller.sr,type="l")

######ivory extension rate Wooller et al 2021####
wooller.COSr<-read.csv("data/Wooller_isotope_data.csv")

wooller.rate <- rep(NA,27)

#first 2 years of life is neonate (Wooller et al 2021)
#last year is close to death, so these data are omitted
for(i in 3:27){
  test.sub.last <- subset(wooller.COSr$d,wooller.COSr$year==(i-1))
  test.sub <- subset(wooller.COSr$d,wooller.COSr$year==i)
  #d is in an increasing order, so get the first element of last year
  test.comb <- rbind(test.sub,test.sub.next[1])
  wooller.rate[i] <- 10000*(max(test.comb)- min(test.comb))/365 #distance in cm, 365 days in a year
  #results are in microns/day
}

wooller.rate<- na.omit(wooller.rate)
mean.wooller.rate <- mean(wooller.rate)
sd.wooller.rate <- sd(wooller.rate)
hist(wooller.rate[3:27])
plot(3:27,wooller.rate[3:27])

mean.wooller.rate*365*28 #this is close to the total length of the tusk at 1.7 meters

####subseting the entire dataset (150 data points)
sub <- 851:1000
sub.mm.sim.avg.dist <- rev(mm.sim.avg.Wooller.dist[sub])
sub.mm.sim.avg.sr <- rev(mm.sim.avg.Wooller.sr[sub])
sub.mm.sim.sd.sr <- rev(mm.sim.sd.Wooller.sr[sub])

plot(sub.mm.sim.avg.dist, sub.mm.sim.avg.sr, type="l", 
     xlim=c(max(sub.mm.sim.avg.dist),min(sub.mm.sim.avg.dist)))
#back to raw data, which is in the dist
Wooller.sub.raw <- subset(Wooller, (wooller.micron> min(sub.mm.sim.avg.dist -250)) & (wooller.micron< 250 + max(sub.mm.sim.avg.dist)))

########################inversion by sampling posterior of 2-pool model########
######inversion of Mammoth record with one pool and parameter#######
Ivo.rate.mean <- mean.wooller.rate #microns per day
Ivo.rate.sd <- sd.wooller.rate

R.sd.mea <- sub.mm.sim.sd.sr
dist.mea <- sub.mm.sim.avg.dist
R.mea <- sub.mm.sim.avg.sr
n.mea = length(sub.mm.sim.avg.sr)

s.intv <- mm.bwidth

max.dist.mea <- max(sub.mm.sim.avg.dist)+ 800 #add some distance before the simulation

#posterior samples of parameters from Misha calibration

a.post <- post.misha.pc2p3$BUGSoutput$sims.list$a[,1]
b.post <- post.misha.pc2p3$BUGSoutput$sims.list$b[,1]
c.post <- post.misha.pc2p3$BUGSoutput$sims.list$c[,1]
post.leng <- length(a.post)

parameters <- c("Ivo.rate","dist", "R1.m","Rin.m", "R2.m","a","b","c","exp.a","exp.bc",
               "Rin.m.pre","a.m","b.m","c.m","Body.mass.m", "Body.mass")

dat = list( s.intv = s.intv, max.dist.mea = max.dist.mea, post.leng=post.leng, 
            a.post=a.post, b.post=b.post, c.post=c.post,
            Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 480, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e3
n.burnin = 2e3
n.thin = 1

#Run it
post.misha.invmamm.param = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param mamm.R", 
                                                      parameters.to.save = parameters, 
                                                      data = dat, n.chains=5, n.iter = n.iter, 
                                                      n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 10 hours

save(post.misha.invmamm.param, file = "out/post.misha.invmamm.param.RData")

post.misha.invmamm.param$BUGSoutput$summary

load("out/post.misha.invmamm.param.RData")

plot(density(post.misha.invmamm.param$BUGSoutput$sims.list$Ivo.rate))

#check prior vs posterior parameters
plot(density(post.misha.pc2p3$BUGSoutput$sims.list$a[,1]), col = "black", lwd = 2, type="l",
     xlim = c(0.01,0.1), xlab = "a", ylab= "density")
#lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$a), col = "blue", lwd = 2)
lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$a.m), col = "red", lwd = 2)
#slight deviation from prior

plot(density(post.misha.pc2p3$BUGSoutput$sims.list$c[,1]), col = "black", lwd = 2, type="l",
     xlim = c(0,0.1), xlab = "c", ylab= "density")
#lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$c), col = "blue", lwd = 2)
lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$c.m), col = "red", lwd = 2)
#slight deviation from prior

subset(post.misha.invmamm.param$BUGSoutput$summary,
       rownames(post.misha.invmamm.param$BUGSoutput$summary)=="a.m")
subset(post.misha.invmamm.param$BUGSoutput$summary,
       rownames(post.misha.invmamm.param$BUGSoutput$summary)=="c.m")

plot(density(post.misha.invmamm.param$BUGSoutput$sims.list$exp.bc))


#do the posterior of a.m and c.m 

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,450), ylim = c(0.706, 0.715), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(sub.mm.sim.avg.dist)+ 800-sub.mm.sim.avg.dist)/mean.wooller.rate,
       sub.mm.sim.avg.sr, pch= 18, col="#00b4ffff")
lines((max(sub.mm.sim.avg.dist)+ 800-sub.mm.sim.avg.dist)/mean.wooller.rate,
      sub.mm.sim.avg.sr, lwd= 1.5, col="#00b4ffff")
#estimated input series
MCMC.ts.Rin.m.invmamm.param.89<- MCMC.CI.bound(post.misha.invmamm.param$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:480,MCMC.ts.Rin.m.invmamm.param.89[[1]],lwd = 2, col = "magenta")
lines(1:480,MCMC.ts.Rin.m.invmamm.param.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:480,MCMC.ts.Rin.m.invmamm.param.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(300, 0.710, c("Measured ivory","Reconstructed input"),lwd = c(1.5, 2), col=c("#00b3ffff","magenta"))


#use this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
########################inversion by sampling posterior of 2-pool model########
######inversion of Mammoth record with remv post and parameter#######
Ivo.rate.mean <- mean.wooller.rate #microns per day
Ivo.rate.sd <- sd.wooller.rate

R.sd.mea <- sub.mm.sim.sd.sr
dist.mea <- sub.mm.sim.avg.dist
R.mea <- sub.mm.sim.avg.sr
n.mea = length(sub.mm.sim.avg.sr)

s.intv <- mm.bwidth

max.dist.mea <- max(sub.mm.sim.avg.dist)+ 800 #add some distance before the simulation

#posterior samples of parameters from Misha calibration

a.post <- post.misha.pc2p.erm$BUGSoutput$sims.list$a[,1]
b.post <- post.misha.pc2p.erm$BUGSoutput$sims.list$b[,1]
c.post <- post.misha.pc2p.erm$BUGSoutput$sims.list$c[,1]
post.leng <- length(a.post)

parameters <- c("Ivo.rate","dist", "R1.m","Rin.m", "R2.m","a","b","c","exp.a","exp.bc",
                "Rin.m.pre","a.m","b.m","c.m","Body.mass.m", "Body.mass")

dat = list( s.intv = s.intv, max.dist.mea = max.dist.mea, post.leng=post.leng, 
            a.post=a.post, b.post=b.post, c.post=c.post,
            Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 480, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 8e3
n.burnin = 4e3
n.thin = floor(n.iter-n.burnin)/500

#Run it
post.misha.invmamm.param.erm = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param mamm.R", 
                                                      parameters.to.save = parameters, 
                                                      data = dat, n.chains=5, n.iter = n.iter, 
                                                      n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 10 hours

save(post.misha.invmamm.param.erm, file = "out/post.misha.invmamm.param.erm.RData")

post.misha.invmamm.param.erm$BUGSoutput$summary

load("out/post.misha.invmamm.param.erm.RData")

plot(density(post.misha.invmamm.param.erm$BUGSoutput$sims.list$Ivo.rate))

#check prior vs posterior parameters
plot(density(post.misha.pc2p.erm$BUGSoutput$sims.list$a[,1]), col = "black", lwd = 2, type="l",
     xlim = c(0.01,0.06), xlab = "a", ylab= "density")
#lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$a), col = "blue", lwd = 2)
lines(density(post.misha.invmamm.param.erm$BUGSoutput$sims.list$a.m), col = "red", lwd = 2)
#slight deviation from prior

plot(density(post.misha.pc2p.erm$BUGSoutput$sims.list$c[,1]), col = "black", lwd = 2, type="l",
     xlim = c(0,1), xlab = "c", ylab= "density")
#lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$c), col = "blue", lwd = 2)
lines(density(post.misha.invmamm.param.erm$BUGSoutput$sims.list$c.m), col = "red", lwd = 2)
#slight deviation from prior

subset(post.misha.invmamm.param.erm$BUGSoutput$summary,
       rownames(post.misha.invmamm.param.erm$BUGSoutput$summary)=="a.m")
subset(post.misha.invmamm.param.erm$BUGSoutput$summary,
       rownames(post.misha.invmamm.param.erm$BUGSoutput$summary)=="c.m")

plot(density(post.misha.invmamm.param.erm$BUGSoutput$sims.list$exp.a))
plot(density(post.misha.invmamm.param.erm$BUGSoutput$sims.list$exp.bc))


#do the posterior of a.m and c.m 

#plotting reconstructed Rin history
#precision term 1e-7, and constrained c#
plot(0,0, xlim = c(1,450), ylim = c(0.706, 0.715), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(sub.mm.sim.avg.dist)+ 800-sub.mm.sim.avg.dist)/mean.wooller.rate,
       sub.mm.sim.avg.sr, pch= 18, col="#00b4ffff")
lines((max(sub.mm.sim.avg.dist)+ 800-sub.mm.sim.avg.dist)/mean.wooller.rate,
      sub.mm.sim.avg.sr, lwd= 1.5, col="#00b4ffff")
#estimated input series
MCMC.ts.Rin.m.invmamm.param.erm.89<- MCMC.CI.bound(post.misha.invmamm.param.erm$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:480,MCMC.ts.Rin.m.invmamm.param.erm.89[[1]],lwd = 2, col = "magenta")
lines(1:480,MCMC.ts.Rin.m.invmamm.param.erm.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:480,MCMC.ts.Rin.m.invmamm.param.erm.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Measured ivory","Reconstructed input"),lwd = c(1.5, 2), col=c("#00b3ffff","magenta"))

########################inversion by sampling posterior of 2-pool model########
######inversion of Mammoth record with remv post and parameter with auto correlation#######
Ivo.rate.mean <- mean.wooller.rate #microns per day
Ivo.rate.sd <- sd.wooller.rate

R.sd.mea <- sub.mm.sim.sd.sr
dist.mea <- sub.mm.sim.avg.dist
R.mea <- sub.mm.sim.avg.sr
n.mea = length(sub.mm.sim.avg.sr)

s.intv <- mm.bwidth

max.dist.mea <- max(sub.mm.sim.avg.dist)+ 800 #add some distance before the simulation

#posterior samples of parameters from Misha calibration

a.post <- post.misha.pc2p.erm$BUGSoutput$sims.list$a[,1]
b.post <- post.misha.pc2p.erm$BUGSoutput$sims.list$b[,1]
c.post <- post.misha.pc2p.erm$BUGSoutput$sims.list$c[,1]
post.leng <- length(a.post)

parameters <- c("Ivo.rate","dist", "R1.m","Rin.m", "R2.m","a","b","c","exp.ab",
                "Rin.m.pre","a.m","b.m","c.m","Body.mass.m", "Body.mass", "Rin.m.cps.ac")

dat = list( s.intv = s.intv, max.dist.mea = max.dist.mea, post.leng=post.leng, 
            a.post=a.post, b.post=b.post, c.post=c.post,
            Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 480, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e3
n.burnin = 2e3
n.thin = 1#floor(n.iter-n.burnin)/500

#Run it
post.misha.invmamm.tsrwca.erm = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param mamm tsrwca.R", 
                                                          parameters.to.save = parameters, 
                                                          data = dat, n.chains=5, n.iter = n.iter, 
                                                          n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 13.5 hours

save(post.misha.invmamm.tsrwca.erm, file = "out/post.misha.invmamm.tsrwca.erm.RData")

post.misha.invmamm.tsrwca.erm$BUGSoutput$summary

load("out/post.misha.invmamm.tsrwca.erm.RData")

plot(density(post.misha.invmamm.tsrwca.erm$BUGSoutput$sims.list$Ivo.rate))

#check prior vs posterior parameters
plot(density(post.misha.pc2p.erm$BUGSoutput$sims.list$a[,1]), col = "black", lwd = 2, type="l",
     xlim = c(0.01,0.1), xlab = "a", ylab= "density")
#lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$a), col = "blue", lwd = 2)
lines(density(post.misha.invmamm.tsrwca.erm$BUGSoutput$sims.list$a.m), col = "red", lwd = 2)
#slight deviation from prior


subset(post.misha.invmamm.tsrwca.erm$BUGSoutput$summary,
       rownames(post.misha.invmamm.tsrwca.erm$BUGSoutput$summary)=="a.m")#
subset(post.misha.invmamm.tsrwca.erm$BUGSoutput$summary,
       rownames(post.misha.invmamm.tsrwca.erm$BUGSoutput$summary)=="Rin.m.cps.ac")#poor convergence on autocorrelation

plot(density(post.misha.invmamm.tsrwca.erm$BUGSoutput$sims.list$exp.ab))
plot(density(post.misha.invmamm.tsrwca.erm$BUGSoutput$sims.list$Rin.m.cps.ac))

#do the posterior of a.m and c.m 

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,450), ylim = c(0.706, 0.715), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(sub.mm.sim.avg.dist)+ 800-sub.mm.sim.avg.dist)/mean.wooller.rate,
       sub.mm.sim.avg.sr, pch= 18, col="#00b4ffff")
lines((max(sub.mm.sim.avg.dist)+ 800-sub.mm.sim.avg.dist)/mean.wooller.rate,
      sub.mm.sim.avg.sr, lwd= 1.5, col="#00b4ffff")
#estimated input series
MCMC.ts.Rin.m.invmamm.tsrwca.erm.89<- MCMC.CI.bound(post.misha.invmamm.tsrwca.erm$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:480,MCMC.ts.Rin.m.invmamm.tsrwca.erm.89[[1]],lwd = 2, col = "magenta")
lines(1:480,MCMC.ts.Rin.m.invmamm.tsrwca.erm.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:480,MCMC.ts.Rin.m.invmamm.tsrwca.erm.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Measured ivory","Reconstructed input"),lwd = c(1.5, 2), col=c("#00b3ffff","magenta"))

plot(density(post.misha.pc2p.erm$BUGSoutput$sims.list$a, from = 0), xlim = c(0.01,0.05),ylim= c(0,120),
     lwd = 2, col = "red",main="a", xlab="parameter estimate")
lines(density(post.misha.invmamm.tsrwca.erm$BUGSoutput$sims.list$a.m, from = 0),lwd = 2,col="blue")

legend(0.04,120, c("Calibration","Fidelity test"),lwd = c(2, 2), col=c("red","blue"))

