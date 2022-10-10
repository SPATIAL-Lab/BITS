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

# n.avg <- 100 #now many data points to average
# n <- floor(length(Wooller$Sr_Seg01)/n.avg) #4287 data points, discard the last <100 data points
# Wooller.sr <- Wooller$Sr_Seg01
# #initiate vectors
# n.avg.Wooller.sr <- rep(0, n)
# n.sd.Wooller.sr <- rep(0, n) 
# 
# n.avg.Wooller.dist <- rep(0, n) 
# n.med.Wooller.dist <- rep(0, n) 
#   
# for(i in 1:n){
#   x <- ((i-1)*n.avg + 1):(i*n.avg)
#   n.avg.Wooller.sr[i] <- mean(Wooller.sr[x])
#   n.sd.Wooller.sr[i] <- sd(Wooller.sr[x])
#   
#   n.avg.Wooller.dist[i] <- mean(wooller.micron[x])#mean and median are the same
# }
# 
# plot(n.avg.Wooller.dist, n.avg.Wooller.sr,type="l")
# lines(n.med.Wooller.dist, n.avg.Wooller.sr, col="red")
# 
# ####subseting the entire dataset
# sub <- 1001:1300
# sub.n.avg.dist <- rev(n.avg.Wooller.dist[sub])
# sub.n.avg.sr <- rev(n.avg.Wooller.sr[sub])
# sub.n.sd.sr <- rev(n.sd.Wooller.sr[sub])
# 
# plot(sub.n.avg.dist, sub.n.avg.sr, type="l", xlim=c(max(sub.n.avg.dist),min(sub.n.avg.dist)))

#foward model simulating micromill results (400 micron band)
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

####subseting the entire dataset (300 is too many, try 150)
sub <- 851:1000
sub.mm.sim.avg.dist <- rev(mm.sim.avg.Wooller.dist[sub])
sub.mm.sim.avg.sr <- rev(mm.sim.avg.Wooller.sr[sub])
sub.mm.sim.sd.sr <- rev(mm.sim.sd.Wooller.sr[sub])

plot(sub.mm.sim.avg.dist, sub.mm.sim.avg.sr, type="l", 
     xlim=c(max(sub.mm.sim.avg.dist),min(sub.mm.sim.avg.dist)))

##body mass scaling is embedded in the JAGS model##
####scaling parameters a, b, c to body mass of the subject#####
# adjusting a, b and c to the body mass of the elephant investigated
# rate ~ scale with e3/4 body mass (basal matabolic rate)
# pool ~ scale with 1 body mass
# a, b and c are rate/pool, so they should scale with -1/4 body mass
# a.m <- a * ((Body.mass.ratio)^-0.25)


#####Inversion using a subset of data from Wooller et al., 2021##############
#####Inversion using params didn't work, so try withfull calibration####
#sample interval
s.intv <- mm.bwidth

max.dist.mea <- max(sub.mm.sim.avg.dist)+ 3000 #add some distance before the simulation

n.mea <- length(sub.mm.sim.avg.sr)

#parameters to save
parameters <- c("Body.mass", "Body.mass.m","Rin.m.cps",
                "Rs.m", "Rb.m","Rin.m","dist", "a.m", "b.m", "c.m","Rin.m.cps.ac")

##Data to pass to the model
#omitting the turnover model here simplifies the model run time
#the inversion takes the measured value of potentially a different ivory series
dat = list(s.intv = s.intv, 
           params.mu = turnover.params.mu, params.vcov = turnover.params.vcov,
           Ivo.rate.mean = mean.wooller.rate, Ivo.rate.sd = sd.wooller.rate,
           max.dist.mea = max.dist.mea,
           R.mea = sub.mm.sim.avg.sr, dist.mea = sub.mm.sim.avg.dist, 
           R.sd.mea = sub.mm.sim.sd.sr, t = 550, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 4e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/400

#Run it
post.mamm.inv.bm = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS BM with parameters Mamm.R", 
                                               parameters.to.save = parameters, 
                                               data = dat, n.chains=4, n.iter = n.iter, 
                                               n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 8 hours

save(post.mamm.inv.bm, file = "out/post.mamm.inv.bm.RData")

post.mamm.inv.bm$BUGSoutput$summary

load("out/post.mamm.inv.bm.RData")

#plotting modeled serum values mapped onto ivory and checking the fit of the data

plot(density(post.mamm.inv.bm$BUGSoutput$sims.list$Rin.m.cps))
plot(density(post.mamm.inv.bm$BUGSoutput$sims.list$Body.mass.m))

#plotting reconstructed Rin.m history
plot(0,0, xlim = c(1,500), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines((max.dist.mea - sub.mm.sim.avg.dist)/mean.wooller.rate + 1, sub.mm.sim.avg.sr, lwd = 2, col = "blue") #approximate results from micromill

MCMC.inv.mamm.Rin.m.89 <- MCMC.CI.bound(post.mamm.inv.bm$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:500,MCMC.inv.mamm.Rin.m.89[[1]],lwd = 2, col = "firebrick4")
lines(1:500,MCMC.inv.mamm.Rin.m.89[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:500,MCMC.inv.mamm.Rin.m.89[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.714, c("Mammoth LA 500 micron Averaged","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))

par(mfrow=c(1,3))

plot(abc.prior.params$x,abc.prior.params$y[1,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "a", ylab= "density")
lines(density(log(post.mamm.inv.bm$BUGSoutput$sims.list$a.m)), col = "red", lwd = 2)

plot(abc.prior.params$x,abc.prior.params$y[2,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "b", ylab= "density")
lines(density(log(post.mamm.inv.bm$BUGSoutput$sims.list$b.m)), col = "red", lwd = 2)

plot(abc.prior.params$x,abc.prior.params$y[3,], col = "blue", lwd = 2, type="l",
     xlim = c(-10,0), xlab = "c", ylab= "density")
lines(density(log(post.mamm.inv.bm$BUGSoutput$sims.list$c.m)), col = "red", lwd = 2)

#####Inversion using a subset of data from Wooller et al., 2021, with full calibration##############
#assign input values before and after the switch, but allows some variation
misha <- read.csv("data/Misha ivory.csv")

R.sd.cal <- misha$sd
dist.cal <- misha$dist
R.cal <- misha$mean
n.cal = length(dist.cal)

R0 <- 0.70706

#Re is the mean ratio of end value  
Re <- 0.7112

max.dist.cal <- 19200
cal.intv <- 100

#calibration rate
Ivo.rate.cal.mean <- 14.7 #microns per day
Ivo.rate.cal.sd <- 0.6

s.intv <- mm.bwidth

max.dist.mea <- max(sub.mm.sim.avg.dist)+ 8000 #add some distance before the simulation

n.mea <- length(sub.mm.sim.avg.sr)

#parameters to save
parameters <- c("Ivo.rate", "Rs.cal", "Rb.cal", "dist.cal.m","dist","Rin.cal","a","b","c",
                "Rs.m", "Rb.m","Rin.m", "a.m", "b.m", "c.m","Rin.m.cps.ac","Ivo.rate.cal")

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
           R.sd.mea = sub.mm.sim.sd.sr, t = 550, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 2e3
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

MCMC.dist.plot(post.mamm.inv.woi$BUGSoutput$sims.list$Rs.cal,
               post.mamm.inv.woi$BUGSoutput$sims.list$dist.cal.m)
lines(misha$dist,misha$mean,lwd = 2, col = "red")

post.mamm.inv.woi.Rs.cal.89<- MCMC.CI.bound(post.mamm.inv.woi$BUGSoutput$sims.list$Rs.cal, 0.89)

#plotting reconstructed Rin.m history
plot(0,0, xlim = c(1,500), ylim = c(0.706, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines((max.dist.mea - sub.mm.sim.avg.dist)/mean.wooller.rate + 1, sub.mm.sim.avg.sr, lwd = 2, col = "blue") #approximate results from micromill

MCMC.inv.mamm.woi.Rin.m.89 <- MCMC.CI.bound(post.mamm.inv.woi$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:550,MCMC.inv.mamm.woi.Rin.m.89[[1]],lwd = 2, col = "firebrick4")
lines(1:550,MCMC.inv.mamm.woi.Rin.m.89[[2]], lwd = 1, lty = 2, col = "firebrick4")
lines(1:550,MCMC.inv.mamm.woi.Rin.m.89[[3]], lwd = 1, lty = 2, col = "firebrick4")
legend(0, 0.716, c("Mammoth LA 500 micron Averaged","Reconstructed input"),lwd = c(2, 2), col=c("blue","firebrick4"))
