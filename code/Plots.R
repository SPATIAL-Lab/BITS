library(Cairo)
library(scales)
library(viridisLite)
library(vioplot)

plot.col.6 <- inferno(6)
########Figure 2 ##############
#######plotting 25 pt average#######
#600 * 400
plot(n.avg.misha.25.dist, n.avg.misha.25.sr,col="#00b4ffff",type = "l",lwd=2,
     xlim=c(20000,8000),ylim=c(0.706,0.712),main="LA-ICP-MS 25 pt average",
     xlab="distance (micron) from pulp cavity", ylab="Sr 87/86")

#######plotting 50 pt average#######
#540 * 360
plot(n.avg.misha.50.dist, n.avg.misha.50.sr,col="#00b4ffff",type = "l",lwd=2,
     xlim=c(20000,8000),ylim=c(0.706,0.712),main="LA-ICP-MS 50 pt average",
     xlab="distance (micron) from pulp cavity", ylab="Sr 87/86")

########Figure 3 ##############
#######plotting 25 pt average#######
#Panel a: 800 * 560

svg(filename = "out/Fig 3 raw.svg",
         width = 8, height = 5.6, pointsize = 12)
plot(misha.raw$dist, misha.raw$X87Sr.86Sr, col=alpha("black",0.2),pch=16,
     xlim=c(20000,8000),ylim=c(0.705,0.714),main="LA-ICP-MS raw data & 25 pt average",
     xlab="distance (micron) from pulp cavity", ylab="Sr 87/86")#plot all raw data points
lines(n.avg.misha.25.dist, n.avg.misha.25.sr,col="#00b4ffff",lwd=1.5)
points(n.avg.misha.25.dist, n.avg.misha.25.sr, col="#00b4ffff",pch=18)
dev.off()
#Panel b: range of measured water and food values?
intake <- read.csv("data/intake.csv")

#######plotting 50 pt average#######
#Panel a: 800 * 560

svg(filename = "out/Fig 3 raw50.svg",
    width = 8, height = 5.6, pointsize = 12)
plot(misha.raw$dist, misha.raw$X87Sr.86Sr, col=alpha("black",0.2),pch=16,
     xlim=c(20000,8000),ylim=c(0.705,0.714),main="LA-ICP-MS raw data & 50 pt average",
     xlab="distance (micron) from pulp cavity", ylab="Sr 87/86")#plot all raw data points
lines(n.avg.misha.50.dist, n.avg.misha.50.sr,col="#00b4ffff",lwd=1.5)
points(n.avg.misha.50.dist, n.avg.misha.50.sr, col="#00b4ffff",pch=18)
dev.off()
#Panel b: range of measured water and food values?
intake <- read.csv("data/intake.csv")

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
#concentrations in bar plots?
vioplot(post.intake, ylab="Sr intake (mg/day)", main ="Daily Sr intake",col=c(plot.col[5:2]),ylim=c(0,1200))
axis(4,at=c(0, 200, 400, 600) )
# intake.naomit<-na.omit(intake) #for concentration values
# stripchart(Sr_conc ~ type, data = intake.naomit, ylim=c(0,120), vertical=TRUE, method = "jitter", pch=18,
#            cex=1.2,col=plot.col.6[2],main = "", xlab = "", ylab = "",axes=F)
# int.med <- tapply(intake.naomit[,"Sr_conc"], intake.naomit[,"type"], median)
# ## Now add line segments corresponding to the group-wise medians
# loc <- 1:length(int.med)
# segments(loc-0.25, int.med, loc+0.25, int.med, col=plot.col.6[2], lwd=3)
# axis(4,at=c(0, 20, 40, 60) )
########Figure 4 ##############

###########50 pt average curve
load("out/post.misha.pc2p.RData")
load("out/post.misha.1p50r.RData")
#calibration curves of 2 pool vs 1 pool model
#posterior distribution of parameters 
#Two pool first 50 pt average
svg(filename = "out/Fig 4 2pcal.svg",
    width = 8, height = 5.6, pointsize = 12)
par(mar = c(5.1, 4.1, 4.1, 4.1))
plot(0,0, xlim = c(20000,8000), ylim = c(0.7066, 0.712), xlab = "distance", ylab ="Sr 87/86",
     main="Two-pool model")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)
#4000 lines are too many
#further thinning to 2000 lines
ind.pc2p<- sample(dim(post.misha.pc2p$BUGSoutput$sims.list$Rs.m)[1],500,replace = F)
MCMC.dist.plot(post.misha.pc2p$BUGSoutput$sims.list$Rs.m[ind.pc2p,],
               post.misha.pc2p$BUGSoutput$sims.list$dist[ind.pc2p,])
points(n.avg.misha.50.dist[index.50.anom.remv1], n.avg.misha.50.sr[index.50.anom.remv1],
       pch=18, col = "#00b4ffff")
points(n.avg.misha.50.dist[index.50.anom.remv2], n.avg.misha.50.sr[index.50.anom.remv2],
       pch=18, col = "#00b4ffff")

lines(n.avg.misha.50.dist[index.50.anom.remv1], n.avg.misha.50.sr[index.50.anom.remv1],
      lwd=1.5, col = "#00b4ffff")
lines(n.avg.misha.50.dist[index.50.anom.remv2], n.avg.misha.50.sr[index.50.anom.remv2],
      lwd=1.5, col = "#00b4ffff")
dev.off()
# then one pool, also use 50 pt average, run line 600
svg(filename = "out/Fig 4 1pcal.svg",
    width = 8, height = 5.6, pointsize = 12)
par(mar = c(5.1, 4.1, 4.1, 4.1))
plot(0,0, xlim = c(20000,8000), ylim = c(0.7066, 0.712), xlab = "distance", ylab ="Sr 87/86",
     main="One-pool model")
abline(h = R0, lwd = 2, lty = 2)
abline(h = Re, lwd = 2, lty = 2)

ind.1p50r<- sample(dim(post.misha.1p50r$BUGSoutput$sims.list$Rs.m)[1],500,replace = F)
MCMC.dist.plot(post.misha.1p50r$BUGSoutput$sims.list$Rs.m[ind.1p50r,],
               post.misha.1p50r$BUGSoutput$sims.list$dist[ind.1p50r,])
points(n.avg.misha.50.dist[index.50.anom.remv1], n.avg.misha.50.sr[index.50.anom.remv1],
       pch=18, col = "#00b4ffff")
points(n.avg.misha.50.dist[index.50.anom.remv2], n.avg.misha.50.sr[index.50.anom.remv2],
       pch=18, col = "#00b4ffff")

lines(n.avg.misha.50.dist[index.50.anom.remv1], n.avg.misha.50.sr[index.50.anom.remv1],
      lwd=1.5, col = "#00b4ffff")
lines(n.avg.misha.50.dist[index.50.anom.remv2], n.avg.misha.50.sr[index.50.anom.remv2],
      lwd=1.5, col = "#00b4ffff")
dev.off()
#Panel b: compare parameters a, b, and c, also use 50 pt with excursion removed to estimate these 3 params.
#first 3 params of the 2 pool model
svg(filename = "out/Fig 4 params.svg",
    width = 8, height = 5.6, pointsize = 12)
par(mar = c(5.1, 4.1, 4.1, 4.1))
par(mfrow=c(1,2))
plot(density(post.misha.pc2p$BUGSoutput$sims.list$a, from = 0), xlim = c(0,0.06),ylim= c(0,150),
     lwd = 2, col = plot.col.6[3],main="Two-pool model", xlab="parameter estimate")
lines(density(post.misha.pc2p$BUGSoutput$sims.list$b, from = 0),lwd = 2, col = plot.col.6[4])
lines(density(post.misha.pc2p$BUGSoutput$sims.list$c, from = 0),lwd = 2, col = plot.col.6[5])
legend(0.02, 150, c("a, Two-pool","b, Two-pool", "c, Two-pool", "a, One-pool"),
       lwd = rep(2, 4), col=c(plot.col.6[3:5],plot.col.6[2]))
#perhaps use pool size, instead
#second, parameter a of the 1 pool model
plot(density(post.misha.1p50r$BUGSoutput$sims.list$a),xlim = c(0,0.06),ylim= c(0,150),
     lwd = 2, col = plot.col.6[2],main="One-pool model", xlab="parameter estimate")
dev.off()

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
########Figure 5 ##############
#panel a: pool ratio vs flux ratio contour map
svg(filename = "out/Fig 5 bivar.svg",
    width = 6, height = 6, pointsize = 12)
par(mfrow=c(1,1))
par(mar = c(5.1, 5.1, 5.1, 5.1))
# make a contour map of flux ratio and pool ratio
contour.flux.pool <- kde2d(post.misha.pc2p$BUGSoutput$sims.list$flux.ratio[,1],
                           post.misha.pc2p$BUGSoutput$sims.list$pool.ratio[,1], n = 64,
                           lims = c(c(0,8), c(0,1.6)))
image(contour.flux.pool, col=viridis(64), xlab="Flux ratio: Fs/Fb", ylab="Pool ratio: Ps/Pb")
contour(contour.flux.pool,lwd = 1.5, add = TRUE, labcex = 1)
dev.off()

#present half life in a table?
##panel b, plot half life
#calculate half-life for pool 1 and pool 2
hl.p1<- map_estimate(log(2)/post.misha.pc2p$BUGSoutput$sims.list$a)
hl.p2<- map_estimate(log(2)/post.misha.pc2p$BUGSoutput$sims.list$c)

par(mfrow=c(1,2))
par(mar = c(5.1, 4.1, 5.1, 3.1))
plot(density(log(2)/post.misha.pc2p$BUGSoutput$sims.list$a, from = 0, to = 50) ,lwd = 2, col = plot.col.6[3],
     ylim = c(0,0.08),main= "Half-life: pool 1",xlab="days")
abline(v=hl.p1[[1]],lty=2,lwd=2)
plot(density(log(2)/post.misha.pc2p$BUGSoutput$sims.list$c, from = 0, to = 600) ,lwd = 2, col = plot.col.6[5], 
     xlim= c(0,600), ylim = c(0,0.005),main= "Half-life: pool 2 ",xlab="days")
abline(v=hl.p2[[1]],lty=2,lwd=2)


#parameter c is bone pool turnover rate, and associated half-life
plot(density(post.misha.pc2p$BUGSoutput$sims.list$c))
plot(density(log(2)/post.misha.pc2p$BUGSoutput$sims.list$c))

######Figure 6 Results forward model###########
#Results reconstructed input using forward model simulated ivory (serum) data
#panel a, synthetic input data with 120 day long intervals (~80% of the variation)
#panel b, synthetic input data with 60 day shorter intervals
#panel c, synthetic input data with 120 day and 60 day intermediate values, to show carry over effect

#Results reconstructed input using Misha LA-ICP-MS series (without excursion, 50 pt)
#days in x-axis, posterior of input (in calibration) vs reconstructed input (with excursion)
#plotting reconstructed Rin history
plot(0,0, xlim = c(1,750), ylim = c(0.705, 0.714), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate

#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.pc2p.Rin.89<- MCMC.CI.bound(post.misha.pc2p$BUGSoutput$sims.list$Rin, 0.89)
lines(1:750,post.misha.pc2p.Rin.89[[1]],lwd = 2, col = "black")
lines(1:750,post.misha.pc2p.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:750,post.misha.pc2p.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

points((max(n.avg.misha.50.dist) + 30-n.avg.misha.50.dist.rmv)/14.7,n.avg.misha.50.sr.rmv, pch= 18, col="#00b3ffff")

#estimated input series
MCMC.misha.inv2pr.param.Rin.m.89 <- MCMC.CI.bound(post.misha.inv2pr.param$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:750,MCMC.misha.inv2pr.param.Rin.m.89[[1]],lwd = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2pr.param.Rin.m.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2pr.param.Rin.m.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(400, 0.708, c("Model input","Reconstructed input"),lwd = c(2, 2), col=c("black","magenta"))

######Figure 7###########
#refer to file #7 forward model





########Figure 8 ##############
##panel a, raw data and 500 micron average data###
svg(filename = "out/Fig 8 raw.svg",
    width = 8, height = 5.6, pointsize = 12)

plot(Wooller.sub.raw$Dist_Seg01 *10000, Wooller.sub.raw$Sr_Seg01, col=alpha("black",0.2),pch=16,
     xlim=c(500000,420000),ylim=c(0.706,0.715),main="LA-ICP-MS raw data & 25 pt average",
     xlab="distance (micron) from pulp cavity", ylab="Sr 87/86")#plot all raw data points
lines(sub.mm.sim.avg.dist, sub.mm.sim.avg.sr,col="#00b4ffff",lwd=1.5)
points(sub.mm.sim.avg.dist, sub.mm.sim.avg.sr, col="#00b4ffff",pch=18)
dev.off()

##panel b, reconstructed series with days in the x axis
load("out/post.misha.invmamm.param.RData")

#plotting reconstructed Rin history
plot(0,0, xlim = c(1,450), ylim = c(0.706, 0.715), xlab = "days", ylab ="Sr 87/86")
#converting misha distance to days using rate Ivo.rate
points((max(sub.mm.sim.avg.dist)+ 800-sub.mm.sim.avg.dist)/mean.wooller.rate,
       sub.mm.sim.avg.sr, pch= 18, col="#00b4ffff")
lines((max(sub.mm.sim.avg.dist)+ 800-sub.mm.sim.avg.dist)/mean.wooller.rate,
      sub.mm.sim.avg.sr, lwd= 1.5, col="#00b4ffff")
#estimated input series
MCMC.ts.Rin.m.invmamm.param.89<- MCMC.CI.bound(post.misha.invmamm.param$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:500,MCMC.ts.Rin.m.invmamm.param.89[[1]],lwd = 2, col = "magenta")
lines(1:500,MCMC.ts.Rin.m.invmamm.param.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:500,MCMC.ts.Rin.m.invmamm.param.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(0, 0.715, c("Measured ivory","Reconstructed input"),lwd = c(1.5, 2), col=c("#00b3ffff","magenta"))

##panel c, posterior distribution of parameters 1) BM, 2) a, 3) half-life, compared to Misha
#BM
plot(density(post.misha.invmamm.param$BUGSoutput$sims.list$Body.mass.m), 
     xlim = c(2000,9000),ylim= c(0,2.5e-3),lwd=2,col="red")
lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$Body.mass),
      lwd=2, col="blue")
#a
plot(density(post.misha.invmamm.param$BUGSoutput$sims.list$a.m), 
     xlim = c(0.01, 0.06),ylim = c(0,80),lwd=2,col="red")
lines(density(post.misha.pc2p$BUGSoutput$sims.list$a),
      lwd=2, col="blue")
#half-life = ln(2)/rate
plot(density(log(2)/post.misha.invmamm.param$BUGSoutput$sims.list$a.m), 
     xlim = c(20, 60),lwd=2,col="red")
lines(density(log(2)/post.misha.pc2p$BUGSoutput$sims.list$a),
      lwd=2, col="blue")


####supplementary figures######
#######Fig S1 Position of cracks and how they have little effect on measured Sr 87/86######

#######Fig S2 Log-normal goodness of fit for parameter a#######
#check q-q plot for fit
qqPlot(post.misha.fdnb1pr$BUGSoutput$sims.list$a[,1],distribution = "lnorm",
       param.list=list(mean=a.param.nb1pr$parameters[1],sd=a.param.nb1pr$parameters[2]),add.line=T)

#######Fig S3 results of intake model###########
#use a different JAGS file for evaluation!

#######Fig S4 sensitivity test 1 Sensitivity to excursion
#excursion makes rate parameter smaller

#######Fig S5 sensitivity test 2 Sensitivity to switch date
#how do you determine switch date?

######Fig S6: Sensitivity test 3 Sensitivity to precision terms?



