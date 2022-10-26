install.packages("Cairo")
library(Cairo)
library(viridisLite)

plot.col.6 <- inferno(6)
########Figure 2 ##############
#######plotting 25 pt average#######
#600 * 400
plot(n.avg.misha.25.dist, n.avg.misha.25.sr,col="#00b4ffff",type = "l",lwd=2,
     xlim=c(20000,8000),ylim=c(0.706,0.712),main="LA-ICP-MS 25 pt average",
     xlab="distance (micron) from pulp cavity", ylab="Sr 87/86")

########Figure 3 ##############
#######plotting 25 pt average#######
#Panel a: 800 * 560

svg(filename = "out/Fig 3 raw.svg",
         width = 8, height = 5.6, pointsize = 12)
plot(misha.raw$dist, misha.raw$X87Sr.86Sr, col=alpha("black",0.25),pch=16,
     xlim=c(20000,8000),ylim=c(0.705,0.714),main="LA-ICP-MS raw data & 25 pt average",
     xlab="distance (micron) from pulp cavity", ylab="Sr 87/86")#plot all raw data points
lines(n.avg.misha.25.dist, n.avg.misha.25.sr,col="#00b4ffff",lwd=1.5)
points(n.avg.misha.25.dist, n.avg.misha.25.sr, col="#00b4ffff",pch=18)
dev.off()
#Panel b: range of measured water and 
intake <- read.csv("data/intake.csv")


stripchart(X87Sr.86Sr ~ type, data = intake, ylim=c(0.705,0.714), vertical=TRUE, method = "stack", pch=19,
           main = "", xlab = "Intake Type", ylab = "Sr 87/86")
## Then compute the group-wise medians
int.med <- tapply(intake[,"X87Sr.86Sr"], intake[,"type"], median)
## Now add line segments corresponding to the group-wise medians
loc <- 1:length(int.med)
segments(loc-0.3, int.med, loc+0.3, int.med, col="red", lwd=3)


boxplot(data=intake, X87Sr.86Sr ~ type)


########Figure 4 ##############
#calibration curves of 2 pool vs 1 pool model
#posterior distribution of parameters 
#Two pool first
plot(0,0, xlim = c(20000,13000), ylim = c(0.7066, 0.712), xlab = "distance", ylab ="Sr 87/86",
     main="Two-pool model")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.fdnbhr$BUGSoutput$sims.list$Rs.m,
               post.misha.fdnbhr$BUGSoutput$sims.list$dist)
points(n.avg.misha.25.dist[index.25.anom.remv1], n.avg.misha.25.sr[index.25.anom.remv1],
       pch=18, col = "#00b4ffff")
points(n.avg.misha.25.dist[index.25.anom.remv2], n.avg.misha.25.sr[index.25.anom.remv2],
       pch=18, col = "#00b4ffff")

lines(n.avg.misha.25.dist[index.25.anom.remv1], n.avg.misha.25.sr[index.25.anom.remv1],
      lwd=1.5, col = "#00b4ffff")
lines(n.avg.misha.25.dist[index.25.anom.remv2], n.avg.misha.25.sr[index.25.anom.remv2],
      lwd=1.5, col = "#00b4ffff")

# then one pool
plot(0,0, xlim = c(20000,13000), ylim = c(0.7066, 0.712), xlab = "distance", ylab ="Sr 87/86",
     main="One-pool model")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.fdnb1pr$BUGSoutput$sims.list$Rs.m,
               post.misha.fdnb1pr$BUGSoutput$sims.list$dist)
points(n.avg.misha.25.dist[index.25.anom.remv1], n.avg.misha.25.sr[index.25.anom.remv1],
       pch=18, col = "#00b4ffff")
points(n.avg.misha.25.dist[index.25.anom.remv2], n.avg.misha.25.sr[index.25.anom.remv2],
       pch=18, col = "#00b4ffff")

lines(n.avg.misha.25.dist[index.25.anom.remv1], n.avg.misha.25.sr[index.25.anom.remv1],
      lwd=1.5, col = "#00b4ffff")
lines(n.avg.misha.25.dist[index.25.anom.remv2], n.avg.misha.25.sr[index.25.anom.remv2],
      lwd=1.5, col = "#00b4ffff")

#check parameters a, b, and c, also use 25 pt with excursion removed to estimate these 3 params.
#first 3 params of the 2 pool model
plot(density(post.misha.fdnbhr$BUGSoutput$sims.list$a, from = 0), xlim = c(0,0.08),ylim= c(0,120),
     lwd = 2, col = plot.col.6[3])
lines(density(post.misha.fdnbhr$BUGSoutput$sims.list$b, from = 0),lwd = 2, col = plot.col.6[4])
lines(density(post.misha.fdnbhr$BUGSoutput$sims.list$c, from = 0),lwd = 2, col = plot.col.6[5])
#second, parameter a of the 1 pool model
lines(density(post.misha.fdnb1pr$BUGSoutput$sims.list$a),lwd = 2, col = plot.col.6[2])

###MAP estimates, and 89% CI for parameters a, b, and c, report in table 1
MCMC.CI.a <- MCMC.CI.bound(post.misha.woint3$BUGSoutput$sims.list$a, 0.89)
MCMC.CI.b <- MCMC.CI.bound(post.misha.woint3$BUGSoutput$sims.list$b, 0.89)
MCMC.CI.c <- MCMC.CI.bound(post.misha.woint3$BUGSoutput$sims.list$c, 0.89)

########Figure 5 ##############
######panel a, plotting calibration curve, one pool model#####
plot(0,0, xlim = c(20000,13000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.fdnb1pr$BUGSoutput$sims.list$Rs.m,
               post.misha.fdnb1pr$BUGSoutput$sims.list$dist)
lines(n.avg.misha.25.dist[index.25.anom.remv1], n.avg.misha.25.sr[index.25.anom.remv1],
      lwd = 2, col = "#00b4ffff")
lines(n.avg.misha.25.dist[index.25.anom.remv2], n.avg.misha.25.sr[index.25.anom.remv2],
      lwd = 2, col = "#00b4ffff")

##panel b, plot rate estimate and half life
plot(density(post.misha.fdnb1p50$BUGSoutput$sims.list$a))
plot(density(post.misha.fdnb1p50$BUGSoutput$sims.list$h.l))

##panel c, prior and posterior for R0 and Re
plot(density(post.misha.fdnb1p$BUGSoutput$sims.list$R0.mean),col="red",xlim=c(0.7065,0.7075))#posterior
lines(seq(0.706,0.708,0.00001),dnorm(seq(0.706,0.708,0.00001),mean = R0, sd= 1/sqrt(100/2e-6)))#prior

plot(density(post.misha.fdnb1p$BUGSoutput$sims.list$Re.mean),col="red",xlim=c(0.7095,0.7125))#posterior
lines(seq(0.7095,0.7125,0.00001),dnorm(seq(0.7095,0.7125,0.00001),mean = Re, sd= 1/sqrt(100/2e-5)))#prior

######Figure 6###########
#Results reconstructed input using Misha LA-ICP-MS series (with excursion, 25 pt?)
#days in x-axis, posterior of input (in calibration) vs reconstructed input (with excursion)
#plotting reconstructed Rin history
plot(0,0, xlim = c(1,400), ylim = c(0.706, 0.715), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(max.dist.mea)-n.avg.misha.25.dist.rmv)/14.7,n.avg.misha.25.sr.rmv, pch= 18, col="#00b3ffff")
#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.fdnb1pr.Rin.89<- MCMC.CI.bound(post.misha.fdnb1pr$BUGSoutput$sims.list$Rin, 0.89)
lines(1:400,post.misha.fdnb1pr.Rin.89[[1]],lwd = 2, col = "black")
lines(1:400,post.misha.fdnb1pr.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:400,post.misha.fdnb1pr.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

#estimated input series
MCMC.ts.Rin.m.inv1phr.param.89<- MCMC.CI.bound(post.misha.inv1phr.param$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:400,MCMC.ts.Rin.m.inv1phr.param.89[[1]],lwd = 2, col = "magenta")
lines(1:400,MCMC.ts.Rin.m.inv1phr.param.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:400,MCMC.ts.Rin.m.inv1phr.param.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Model input","Reconstructed input"),lwd = c(2, 2), col=c("black","magenta"))

######Figure 7###########
#Results reconstructed input using forward model simulated ivory (serum) data
#panel a, synthetic input data with 120 day and 60 day intermediate values
#panel b, synthetic input data with 60 day shorter intervals
#panel c, synthetic input data with 30 day shortest intervals


########Figure 8 ##############
##panel a, raw data and 500 micron average data###
svg(filename = "out/Fig 8 raw.svg",
    width = 8, height = 5.6, pointsize = 12)

plot(Wooller.sub.raw$wooller.micron, Wooller.sub.raw$Sr_Seg01, col=alpha("black",0.25),pch=16,
     xlim=c(500000,400000),ylim=c(0.705,0.714),main="LA-ICP-MS raw data & 25 pt average",
     xlab="distance (micron) from pulp cavity", ylab="Sr 87/86")#plot all raw data points
lines(sub.mm.sim.avg.dist, sub.mm.sim.avg.sr,col="#00b4ffff",lwd=1.5)
points(sub.mm.sim.avg.dist, sub.mm.sim.avg.sr, col="#00b4ffff",pch=18)
dev.off()

##panel b, reconstructed series with days in the x axis

##panel c, posterior distribution of parameters 1) BM, 2) a, 3) half-life, compared to Misha