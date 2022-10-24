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
#800 * 560

svg(filename = "out/Fig 3 raw.svg",
         width = 8, height = 5.6, pointsize = 12)
plot(misha.raw$dist, misha.raw$X87Sr.86Sr, col=alpha("black",0.25),pch=16,
     xlim=c(20000,8000),ylim=c(0.705,0.714),main="LA-ICP-MS raw data & 25 pt average",
     xlab="distance (micron) from pulp cavity", ylab="Sr 87/86")#plot all raw data points
lines(n.avg.misha.25.dist, n.avg.misha.25.sr,col="#00b4ffff",lwd=1.5)
points(n.avg.misha.25.dist, n.avg.misha.25.sr, col="#00b4ffff",pch=18)
dev.off()
########Figure 4 ##############
######plotting calibration curve#####
plot(0,0, xlim = c(20000,13000), ylim = c(0.706, 0.712), xlab = "distance", ylab ="Sr 87/86")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

MCMC.dist.plot(post.misha.fdnb1pr$BUGSoutput$sims.list$Rs.m,
               post.misha.fdnb1pr$BUGSoutput$sims.list$dist)
lines(n.avg.misha.25.dist[index.25.anom.remv1], n.avg.misha.25.sr[index.25.anom.remv1],
      lwd = 2, col = "#00b4ffff")
lines(n.avg.misha.25.dist[index.25.anom.remv2], n.avg.misha.25.sr[index.25.anom.remv2],
      lwd = 2, col = "#00b4ffff")

#plot rate estimate and half life
plot(density(post.misha.fdnb1p50$BUGSoutput$sims.list$a))
plot(density(post.misha.fdnb1p50$BUGSoutput$sims.list$h.l))

#prior and posterior for R0 and Re
plot(density(post.misha.fdnb1p$BUGSoutput$sims.list$R0.mean),col="red",xlim=c(0.7065,0.7075))#posterior
lines(seq(0.706,0.708,0.00001),dnorm(seq(0.706,0.708,0.00001),mean = R0, sd= 1/sqrt(100/2e-6)))#prior

plot(density(post.misha.fdnb1p$BUGSoutput$sims.list$Re.mean),col="red",xlim=c(0.7095,0.7125))#posterior
lines(seq(0.7095,0.7125,0.00001),dnorm(seq(0.7095,0.7125,0.00001),mean = Re, sd= 1/sqrt(100/2e-5)))#prior
