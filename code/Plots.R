library(Cairo)
library(scales)
library(viridisLite)
library(vioplot)

source("code/1 Helper functions.R")

plot.col.6 <- inferno(6)
########Figure 1(c) ##############
#######plotting 50 pt average#######
#540 * 360
plot(n.avg.misha.50.dist, n.avg.misha.50.sr,col="#00b4ffff",type = "l",lwd=2,
     xlim=c(20000,8000),ylim=c(0.705,0.711),main="LA-ICP-MS 50 pt average",
     xlab="distance (micron) from pulp cavity", ylab="Sr 87/86")

########Figure 3 ##############
#######plotting 50 pt average#######
#Panel a: 800 * 560

svg(filename = "out/Fig 3 raw50.svg",
    width = 8, height = 5.6, pointsize = 12)
plot(misha.raw.dist, misha.raw.Sr, col=alpha("black",0.2),pch=16,
     xlim=c(20000,8000),ylim=c(0.704,0.713),main="LA-ICP-MS raw data & 50 pt average",
     xlab="distance (micron) from pulp cavity", ylab="Sr 87/86")#plot all raw data points
lines(n.avg.misha.50.dist, n.avg.misha.50.sr,col="#00b4ffff",lwd=1.5)
points(n.avg.misha.50.dist, n.avg.misha.50.sr, col="#00b4ffff",pch=18)
dev.off()

#Panel b: range of measured water and food values
load("out/post.misha.intake.RData")

post.intake <- cbind(post.misha.intake$BUGSoutput$sims.list$w.hay[,1],
                     post.misha.intake$BUGSoutput$sims.list$w.pel[,1], 
                     post.misha.intake$BUGSoutput$sims.list$w.sup[,1], 
                     post.misha.intake$BUGSoutput$sims.list$w.water[,1])
colnames(post.intake) <- c("Hay", "Pellet", "Supp.","Water")

##### dot plot for measured food and water values######
#536h * 330w
par(mar = c(5.1, 4.1, 4.1, 4.1))
stripchart(X87Sr.86Sr ~ type, data = intake, ylim=c(0.705,0.714), vertical=TRUE, method = "jitter", pch=19,
           cex=1, col=plot.col.6[4],main = "", xlab = "Intake Type", ylab = "Sr 87/86")
## Then compute the group-wise medians
int.med <- tapply(intake[,"X87Sr.86Sr"], intake[,"type"], median)
## Now add line segments corresponding to the group-wise medians
loc <- 1:length(int.med)
segments(loc-0.25, int.med, loc+0.25, int.med, col=plot.col.6[4], lwd=3)

vioplot(post.intake, ylab="Sr intake (mg/day)", main ="Daily Sr intake",col=c(plot.col[5:2]),ylim=c(0,1200))
axis(4,at=c(0, 200, 400, 600) )

#panel c
load("out/post.misha.pc2p3.RData")

svg(filename = "out/Fig 4 2pcal.svg",
    width = 8, height = 5.6, pointsize = 12)
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(0,0, xlim = c(20000,8000), ylim = c(0.705, 0.711), xlab = "distance", ylab ="Sr 87/86",
     main="Calibration")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
#4000 lines are too many
#further thinning to 2000 lines
ind.pc2p<- sample(dim(post.misha.pc2p3$BUGSoutput$sims.list$R1.m)[1],500,replace = F)
MCMC.dist.plot(post.misha.pc2p3$BUGSoutput$sims.list$R1.m[ind.pc2p,],
               post.misha.pc2p3$BUGSoutput$sims.list$dist[ind.pc2p,])
points(n.avg.misha.50.dist[index.50.anom.remv1], n.avg.misha.50.sr[index.50.anom.remv1],
       pch=18, col = "#00b4ffff")
points(n.avg.misha.50.dist[index.50.anom.remv2], n.avg.misha.50.sr[index.50.anom.remv2],
       pch=18, col = "#00b4ffff")

lines(n.avg.misha.50.dist[index.50.anom.remv1], n.avg.misha.50.sr[index.50.anom.remv1],
      lwd=1.5, col = "#00b4ffff")
lines(n.avg.misha.50.dist[index.50.anom.remv2], n.avg.misha.50.sr[index.50.anom.remv2],
      lwd=1.5, col = "#00b4ffff")
dev.off()

######panel d#####
#536h * 330w
#panel d
par(mfrow=c(3,1))
par(mar = c(4.1, 4.1, 2.1, 4.1))
plot(density(post.misha.pc2p3$BUGSoutput$sims.list$a, from = 0), xlim = c(0,0.04),ylim= c(0,100),
     lwd = 2, col = plot.col.6[3],main="a & b", xlab="Parameter estimates")
lines(density(post.misha.pc2p3$BUGSoutput$sims.list$b, from = 0),
     lwd = 2, col = plot.col.6[4])

legend(0, 100, c("a","b"),
       lwd = rep(2, 2), col=c(plot.col.6[3:4]))

plot(density(post.misha.pc2p3$BUGSoutput$sims.list$c, from = 0, to = 0.01), xlim = c(0.001,0.009),ylim= c(0,420),
     lwd = 2, col = plot.col.6[5],main="c", xlab="Parameter estimate")

######panel e#####
# make a contour map
par(mfrow=c(1,1))
contour.flux.pool <- kde2d(post.misha.pc2p3$BUGSoutput$sims.list$flux.ratio[,1],
                           post.misha.pc2p3$BUGSoutput$sims.list$pool.ratio[,1], n = 64,
                           lims = c(c(0.9,1.8), c(0.1,0.5)))
image(contour.flux.pool, col=viridis(64), xlab="Flux ratio: FI/FII", ylab="Pool ratio: PI/PII")
contour(contour.flux.pool,lwd = 1, add = TRUE, labcex = 0.8)


########Figure 4 ##############
load("out/post.misha.inv2perm.param.RData")
#Results reconstructed input using Misha LA-ICP-MS series (without excursion, 50 pt)
#days in x-axis, posterior of input (in calibration) vs reconstructed input (with excursion)
#plotting reconstructed Rin history
#750*500
plot(0,0, xlim = c(1,700), ylim = c(0.705, 0.714), xlab = "days", ylab ="Sr 87/86",
     main="Fidelity test: model input series vs. estimated input series")
#converting misha distance to days using rate Ivo.rate

#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.pc2p3.Rin.89<- MCMC.CI.bound(post.misha.pc2p3$BUGSoutput$sims.list$Rin, 0.89)
lines(1:750,post.misha.pc2p3.Rin.89[[1]],lwd = 2, col = "black")
lines(1:750,post.misha.pc2p3.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:750,post.misha.pc2p3.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

points((max(n.avg.misha.50.dist) + 30-n.avg.misha.50.dist.rmv)/14.7,n.avg.misha.50.sr.rmv, pch= 18, col="#00b3ffff")

#estimated input series
#estimated input series
MCMC.misha.inv2p3.param.Rin.m.89 <- MCMC.CI.bound(post.misha.inv2p3.param$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:750,MCMC.misha.inv2p3.param.Rin.m.89[[1]],lwd = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2p3.param.Rin.m.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2p3.param.Rin.m.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(400, 0.708, c("Model input","Estimated input"),lwd = c(2, 2), col=c("black","magenta"))

##panel b adding prior vs posterior density plots#####
plot(density(post.misha.pc2p3$BUGSoutput$sims.list$a, from = 0), xlim = c(0,0.04),ylim= c(0,140),
     lwd = 2, col = "red",main="a", xlab="parameter estimate")
lines(density(post.misha.inv2p3.param$BUGSoutput$sims.list$a, from = 0),lwd = 2,col="blue")

legend(0.03,100, c("Calibration","Fidelity test"),lwd = c(2, 2), col=c("red","blue"))

########Figure 5 ##############
load("out/post.misha.invmamm.param.erm.RData")
##panel a, raw data and 500 micron average data###
svg(filename = "out/Fig 5ab.svg",
    width = 8, height = 8, pointsize = 12)
par(mfrow=c(2,1))
par(mar = c(4.1, 4.1, 3.1, 2.1))
plot(Wooller.sub.raw$Dist_Seg01 *10000, Wooller.sub.raw$Sr_Seg01, col=alpha("black",0.2),pch=16,
     xlim=c(500000,425000),ylim=c(0.708,0.715),main="LA-ICP-MS raw data & 500 micron average",
     xlab="distance (micron) from pulp cavity", ylab="Sr 87/86")#plot all raw data points
lines(sub.mm.sim.avg.dist, sub.mm.sim.avg.sr,col="#00b4ffff",lwd=1.5)
points(sub.mm.sim.avg.dist, sub.mm.sim.avg.sr, col="#00b4ffff",pch=18)

##panel b, reconstructed series with days in the x axis
#plotting reconstructed Rin history
plot(0,0, xlim = c(1,440), ylim = c(0.705, 0.714), xlab = "days", ylab ="Sr 87/86",
     main="Estimated input series from 500 micron average")
#converting misha distance to days using rate Ivo.rate
points((max(sub.mm.sim.avg.dist)+ 800-sub.mm.sim.avg.dist)/mean.wooller.rate,
       sub.mm.sim.avg.sr, pch= 18, col="#00b4ffff")
# lines((max(sub.mm.sim.avg.dist)+ 800-sub.mm.sim.avg.dist)/mean.wooller.rate,
#       sub.mm.sim.avg.sr, lwd= 1.5, col="#00b4ffff")
#estimated input series
MCMC.ts.Rin.m.invmamm.param.89<- MCMC.CI.bound(post.misha.invmamm.param$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:480,MCMC.ts.Rin.m.invmamm.param.89[[1]],lwd = 2, col = "magenta")
lines(1:480,MCMC.ts.Rin.m.invmamm.param.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:480,MCMC.ts.Rin.m.invmamm.param.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(300, 0.709, c("500 micron average","Estimated input"),lwd = c(1.5, 2), col=c("#00b3ffff","magenta"))
dev.off()
##panel c, posterior distribution of parameters 1) BM, 2) a, 3) half-life, compared to Misha

svg(filename = "out/Fig 5c.svg",
    width = 8, height = 8, pointsize = 12)
par(mfrow=c(2,3))
par(mar = c(4.1, 4.1, 3.1, 1.1))

#a
plot(density(post.misha.invmamm.param.erm$BUGSoutput$sims.list$a.m), xlim = c(0.01, 0.05),ylim = c(0,120),
     lwd=2,col="red", main ="Posterior densities: a", xlab="Parameter estimate")
lines(density(post.misha.pc2p.erm$BUGSoutput$sims.list$a),
      lwd=2, col="blue")
legend(0.04,100, c("Calibration","Case study"),lwd = c(2, 2), col=c("blue","red"))

#b
plot(density(post.misha.invmamm.param.erm$BUGSoutput$sims.list$b.m, from = 0), xlim = c(0, 0.05),ylim = c(0,60),
     lwd=2,col="red", main ="Posterior densities: b", xlab="Parameter estimate")
lines(density(post.misha.pc2p.erm$BUGSoutput$sims.list$b, from = 0),
      lwd=2, col="blue")
legend(0.04,100, c("Calibration","Case study"),lwd = c(2, 2), col=c("blue","red"))

#c
plot(density(post.misha.invmamm.param.erm$BUGSoutput$sims.list$c.m, from = 0, to = 1), xlim = c(0, 1),ylim = c(0,1.8),
     lwd=2,col="red", main ="Posterior densities: c", xlab="Parameter estimate")
lines(density(post.misha.pc2p.erm$BUGSoutput$sims.list$c, from = 0, to = 1),
      lwd=2, col="blue")
legend(1.8,0.6, c("Calibration","Case study"),lwd = c(2, 2), col=c("blue","red"))
dev.off()


################supp. one-pool vs two-pool model comparison#######
#parameter a#
plot(density(post.misha.1p50r$BUGSoutput$sims.list$a),xlim = c(0.01,0.05),ylim= c(0,120),
     lwd = 2, col = "red",main="Parameter a comparison", xlab="Parameter estimate")
lines(density(post.misha.pc2p.erm$BUGSoutput$sims.list$a, from = 0),
     lwd = 2, col = "blue")


##panel c, prior and posterior for Re between the models

svg(filename = "out/Fig 4 post2p.svg",
    width = 8, height = 5.6, pointsize = 12)
par(mar = c(5.1, 4.1, 4.1, 4.1))
par(mfrow=c(1,2))
plot(density(post.misha.pc2p$BUGSoutput$sims.list$Re.mean),col="red",lwd=2,
     xlim=c(0.7095,0.7125),ylim=c(0,6e3), xlab="parameter estimate", main="Two-pool model")#posterior
lines(seq(0.7095,0.7125,0.00001),dnorm(seq(0.7095,0.7125,0.00001),mean = 0.711, sd= 1/sqrt(100/2e-5)))#prior

plot(density(post.misha.1p50r$BUGSoutput$sims.list$Re.mean),col="red",lwd=2,
     xlim=c(0.7095,0.7125),ylim=c(0,6e3), xlab="parameter estimate", main="One-pool model")#posterior
lines(seq(0.7095,0.7125,0.00001),dnorm(seq(0.7095,0.7125,0.00001),mean = 0.711, sd= 1/sqrt(100/2e-5)))#prior
dev.off()

######Figure 6 Results forward model###########
#Results reconstructed input using forward model simulated ivory (serum) data
svg(filename = "out/Fig 6 fwd.svg",
    width = 6, height = 8, pointsize = 12)
par(mfrow=c(4,1))
par(mar = c(4.1, 4.1, 3.1, 2.1))
#panel a, synthetic input data with 150-day long intervals and 30-day short intervals (~85% of the variation)
plot(0,0, xlim = c(1,720), ylim = c(syn.low, syn.high), xlab = "days", ylab ="Sr 87/86",
     main="150 + 30-day annual switch")
lines(1:length(syn.input.150), syn.input.150)
lines(1:length(syn.input.150), Se.bone.res150[[1]],lwd=2, col = "#00b4ffff")
#panel b, synthetic input data with 120-day long intervals and 60-day short intervals (~80% of the variation)
plot(0,0, xlim = c(1,720), ylim = c(syn.low, syn.high), xlab = "days", ylab ="Sr 87/86",
     main="120 + 60-day annual switch")
lines(1:length(syn.input.120), syn.input.120)
lines(1:length(syn.input.120), Se.bone.res120[[1]],lwd=2, col = "#00b4ffff")
#panel c, synthetic input data with 60-day long intervals and 30-day short intervals semiannual (~67%)
plot(0,0, xlim = c(1,720), ylim = c(syn.low, syn.high), 
     xlab = "days", ylab ="Sr 87/86",main="60 + 30-day semiannual switch")
lines(1:length(syn.input.60), syn.input.60)
lines(1:length(syn.input.60), Se.bone.res60[[1]],lwd=2, col = "#00b4ffff")
#panel d, synthetic input data with 90 day intermediate invervals semiannual, to show carry over effect
plot(0,0, xlim = c(1,720), ylim = c(syn.low, syn.high), 
     xlab = "days", ylab ="Sr 87/86",main="90 day semiannual switch with carryover")
lines(1:length(syn.input.90), syn.input.90)
lines(1:length(syn.input.90), Se.bone.res90[[1]],lwd=2, col = "#00b4ffff")
dev.off()


######Figure comparing laser ablation and micromill results######
Misha.mm.4.6 <- read.csv("data/Misha micromill 400-600.csv")
n.misha.4.6 <- dim(Misha.mm.4.6)
par(mfrow=c(1,2))
plot(n.misha.4.6[1]:1, Misha.mm.4.6$X87Sr.86Sr,ylim=c(0.706,0.712), type = "l",  
     main ="Misha micromill", xlab="sample #",ylab="Sr 87/86")

plot(n.avg.misha.50.dist, n.avg.misha.50.sr,col="#00b4ffff",type = "l",lwd=2,
     xlim=c(20000,8000),ylim=c(0.706,0.712),main="LA-ICP-MS 50 pt average",
     xlab="distance (micron) from pulp cavity", ylab="Sr 87/86")