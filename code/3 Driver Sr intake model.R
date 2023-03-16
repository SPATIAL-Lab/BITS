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

#adding the model component of food and water mixture as intake#
intake <- read.csv("data/intake.csv")
intake$Sr_conc

Hay <- subset(intake, intake$type=="H")
Pellet <- subset(intake, intake$type=="P")
Supplement <- subset(intake, intake$type=="S")
Water <- subset(intake, intake$type=="W")

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

#estimate mean and sd of posterior distribution of Re from the two-pool model
#normal distribution
R.intake <- enorm(post.misha.pc2p$BUGSoutput$sims.list$Re.mean[,1])
R.intake.mean <- R.intake$parameters[1]
R.intake.sd <- R.intake$parameters[2]

####alternative mean and sd estimates
R.intake.mean <- intake.af
R.intake.sd <- 0.0005
########inversion model with calibration for either the one/two pool model#######
#mixing of 10 kg of hay in diet
m.feed <- 10

#parameters to save
parameters <- c("w.contrib","h.contrib", "w.food", "w.hay", "w.pel", "w.sup", "w.water",
                "f.h", "f.pel", "Rin")

##Data to pass to the model
dat = list(Sr.hay.mean = Sr.hay.mean, Sr.hay.sd = Sr.hay.sd, 
           Sr.pel.mean = Sr.pel.mean, Sr.pel.sd = Sr.pel.sd, 
           Sr.sup.mean = Sr.sup.mean, Sr.sup.sd = Sr.sup.sd,
           Sr.w.mean = Sr.w.mean, Sr.w.sd = Sr.w.sd, m.feed = m.feed,
           conc.hay.mean = conc.hay.mean, conc.hay.sd = conc.hay.sd, 
           conc.pel.mean = conc.pel.mean, conc.pel.sd = conc.pel.sd,
           conc.sup.mean = conc.sup.mean, conc.sup.sd = conc.sup.sd,
           conc.w.mean = conc.w.mean, conc.w.sd = conc.w.sd,
           R.intake.mean = R.intake.mean, R.intake.sd = R.intake.sd)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 5e3
n.thin = 1

#Run it
post.misha.intake = do.call(jags.parallel,list(model.file = "code/Sr intake JAGS.R", 
                                                   parameters.to.save = parameters, 
                                                   data = dat, n.chains=5, n.iter = n.iter, 
                                                   n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #instantaneous

save(post.misha.intake, file = "out/post.misha.intake.RData")

post.misha.intake$BUGSoutput$summary

w.intake.map <- map_estimate(post.misha.intake$BUGSoutput$sims.list$w.contrib[,1])
w.intake.hid <- hdi(post.misha.intake$BUGSoutput$sims.list$w.contrib[,1],0.89)
w.intake.hid[2]
w.intake.hid[3]
