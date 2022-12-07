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

a.post <- post.misha.pc2p$BUGSoutput$sims.list$a[,1]
b.post <- post.misha.pc2p$BUGSoutput$sims.list$b[,1]
c.post <- post.misha.pc2p$BUGSoutput$sims.list$c[,1]
post.leng <- length(a.post)

parameters <- c("Ivo.rate","dist", "R1.m","Rin.m", "R2.m","a","b","c","exp.ab",
               "Rin.m.pre","a.m","b.m","c.m","Body.mass.m", "Body.mass")

dat = list( s.intv = s.intv, max.dist.mea = max.dist.mea, post.leng=post.leng, 
            a.post=a.post, b.post=b.post, c.post=c.post,
            Ivo.rate.mean = Ivo.rate.mean, Ivo.rate.sd = Ivo.rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 480, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 3e3
n.burnin = 1e3
n.thin = floor(n.iter-n.burnin)/500

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
plot(density(post.misha.pc2p$BUGSoutput$sims.list$a[,1]), col = "black", lwd = 2, type="l",
     xlim = c(0.01,0.1), xlab = "a", ylab= "density")
#lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$a), col = "blue", lwd = 2)
lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$a.m), col = "red", lwd = 2)
#slight deviation from prior

plot(density(post.misha.pc2p$BUGSoutput$sims.list$c[,1]), col = "black", lwd = 2, type="l",
     xlim = c(0,0.1), xlab = "c", ylab= "density")
#lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$c), col = "blue", lwd = 2)
lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$c.m), col = "red", lwd = 2)
#slight deviation from prior

subset(post.misha.invmamm.param$BUGSoutput$summary,
       rownames(post.misha.invmamm.param$BUGSoutput$summary)=="a.m")
subset(post.misha.invmamm.param$BUGSoutput$summary,
       rownames(post.misha.invmamm.param$BUGSoutput$summary)=="c.m")

plot(density(post.misha.invmamm.param$BUGSoutput$sims.list$exp.ab))


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
legend(0, 0.715, c("Measured ivory","Reconstructed input"),lwd = c(1.5, 2), col=c("#00b3ffff","magenta"))

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

parameters <- c("Ivo.rate","dist", "R1.m","Rin.m", "R2.m","a","b","c","exp.ab",
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
     xlim = c(0.01,0.1), xlab = "a", ylab= "density")
#lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$a), col = "blue", lwd = 2)
lines(density(post.misha.invmamm.param.erm$BUGSoutput$sims.list$a.m), col = "red", lwd = 2)
#slight deviation from prior

plot(density(post.misha.pc2p.erm$BUGSoutput$sims.list$c[,1]), col = "black", lwd = 2, type="l",
     xlim = c(0,0.1), xlab = "c", ylab= "density")
#lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$c), col = "blue", lwd = 2)
lines(density(post.misha.invmamm.param.erm$BUGSoutput$sims.list$c.m), col = "red", lwd = 2)
#slight deviation from prior

subset(post.misha.invmamm.param.erm$BUGSoutput$summary,
       rownames(post.misha.invmamm.param.erm$BUGSoutput$summary)=="a.m")
subset(post.misha.invmamm.param.erm$BUGSoutput$summary,
       rownames(post.misha.invmamm.param.erm$BUGSoutput$summary)=="c.m")

plot(density(post.misha.invmamm.param.erm$BUGSoutput$sims.list$exp.ab))


#do the posterior of a.m and c.m 

#plotting reconstructed Rin history
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

#####Inversion using a subset of data from Wooller et al., 2021, with full calibration##############
#assign input values before and after the switch, but allows some variation

s.intv <- 26.2

max.dist.mea <- max(n.avg.misha.25.dist) + 30
R.sd.cal <- n.sd.misha.25.sr.rmv
dist.cal <- n.avg.misha.25.dist.rmv
R.cal <- n.avg.misha.25.sr.rmv
n.cal = length(n.avg.misha.25.sr.rmv)

R0 <- 0.707

#Re is the mean ratio of end value  
Re <- 0.711

max.dist.cal <- max(n.avg.misha.25.dist) + 30
cal.intv <- 26.2

#calibration rate
Ivo.rate.cal.mean <- 14.7 #microns per day
Ivo.rate.cal.sd <- 0.6

s.intv <- mm.bwidth

max.dist.mea <- max(sub.mm.sim.avg.dist)+ 800 #add some distance before the simulation

n.mea <- length(sub.mm.sim.avg.sr)

#parameters to save
parameters <- c("Ivo.rate", "R1.cal", "dist.cal.m","dist","Rin.cal","a",
                "R1.m","Rin.m", "a.m","Rin.m.cps.ac","Ivo.rate.cal","Body.mass.m", "Body.mass")

##Data to pass to the model
#compared to the turnover model that is essentially the .cal part here 
#the inversion takes the measured value of potentially a different ivory series
dat = list(R.cal = R.cal, dist.cal = dist.cal, R.sd.cal = R.sd.cal, t.cal = 750, n.cal = n.cal, 
           R0 = R0, Re = Re, cal.intv = cal.intv, s.intv = s.intv,
           max.dist.cal = max.dist.cal, 
           Ivo.rate.cal.mean = Ivo.rate.cal.mean, Ivo.rate.cal.sd = Ivo.rate.cal.sd,
           Ivo.rate.mean = mean.wooller.rate, Ivo.rate.sd = sd.wooller.rate,
           max.dist.mea = max.dist.mea,
           R.mea = sub.mm.sim.avg.sr, dist.mea = sub.mm.sim.avg.dist, 
           R.sd.mea = sub.mm.sim.sd.sr, t = 500, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/400

post.mamm.inv.woi = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS wo intake mamm.R", 
                                               parameters.to.save = parameters, 
                                               data = dat, n.chains=5, n.iter = n.iter, 
                                               n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 61 hours!

save(post.mamm.inv.woi, file = "out/post.mamm.inv.woi.RData")

post.mamm.inv.woi$BUGSoutput$summary

load("out/post.mamm.inv.woi.RData")

plot(density(post.mamm.inv.woi$BUGSoutput$sims.list$a.m)) #check parameters
plot(density(post.mamm.inv.bm$BUGSoutput$sims.list$Body.mass.m))

#plot calibration curve
par(mfrow=c(1,1))
plot(0,0, xlim = c(20000,8000), ylim = c(0.706, 0.713), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.mamm.inv.woi$BUGSoutput$sims.list$R1.cal,
               post.mamm.inv.woi$BUGSoutput$sims.list$dist.cal.m)
lines(misha$dist,misha$mean,lwd = 2, col = "red")

post.mamm.inv.woi.R1.cal.89<- MCMC.CI.bound(post.mamm.inv.woi$BUGSoutput$sims.list$R1.cal, 0.89)

#plotting reconstructed Rin.m history
plot(0,0, xlim = c(1,550), ylim = c(0.706, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines((max.dist.mea - sub.mm.sim.avg.dist)/mean.wooller.rate + 1, sub.mm.sim.avg.sr, lwd = 2, col = "blue") #approximate results from micromill

MCMC.inv.mamm.woi.Rin.m.89 <- MCMC.CI.bound(post.mamm.inv.woi$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:550,MCMC.inv.mamm.woi.Rin.m.89[[1]],lwd = 2, col = "firebrick4")
lines(1:550,MCMC.inv.mamm.woi.Rin.m.89[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:550,MCMC.inv.mamm.woi.Rin.m.89[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Mammoth LA 500 micron averaged","Reconstructed intake"),lwd = c(2, 2), col=c("blue","firebrick4"))

###plot possible bone history###
plot(0,0, xlim = c(1,550), ylim = c(0.706, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines((max.dist.mea - sub.mm.sim.avg.dist)/mean.wooller.rate + 1, sub.mm.sim.avg.sr, lwd = 2, col = "blue") #approximate results from micromill

MCMC.inv.mamm.woi.R2.m.89 <- MCMC.CI.bound(post.mamm.inv.woi$BUGSoutput$sims.list$R2.m, 0.89)
lines(1:550,MCMC.inv.mamm.woi.R2.m.89[[1]],lwd = 2, col = "firebrick2")
lines(1:550,MCMC.inv.mamm.woi.R2.m.89[[2]], lwd = 1, lty = 2, col = "firebrick2")
lines(1:550,MCMC.inv.mamm.woi.R2.m.89[[3]], lwd = 1, lty = 2, col = "firebrick2")
legend(0, 0.716, c("Mammoth LA 500 micron averaged","Reconstructed bone pool"),lwd = c(2, 2), col=c("blue","firebrick2"))

par(mfrow=c(1,3))

plot(abc.prior.params$x,abc.prior.params$y[1,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "a", ylab= "density")
lines(density(log(post.mamm.inv.woi$BUGSoutput$sims.list$a)), col = "red", lwd = 2)

plot(abc.prior.params$x,abc.prior.params$y[2,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "b", ylab= "density")
lines(density(log(post.mamm.inv.woi$BUGSoutput$sims.list$b)), col = "red", lwd = 2)

plot(abc.prior.params$x,abc.prior.params$y[3,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "c", ylab= "density")
lines(density(log(post.mamm.inv.woi$BUGSoutput$sims.list$c)), col = "red", lwd = 2)