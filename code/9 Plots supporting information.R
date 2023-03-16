library(coda)
library(rjags)
library(R2jags)
library(mcmcplots)
library(bayestestR)
library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)
library(zoo)

plot.col<-viridis(7)

setwd("C:/Users/ydmag/Google Drive/U of U/Elephant movement/Sr-in-ivory")

source("code/1 Helper functions.R")

####supplementary figures######

#################Fig S1 ################
#using Mass 88(Sr) as an indicator of cracks in the ivory slab
misha.elements <- read.csv("data/Misha elements.csv")
#the beginning of the transet needs to be cleaned
misha.elements.cl <- misha.elements[25:nrow(misha.elements),]

#panels a-c
par(mfrow=c(3,1))
plot(n.avg.misha.50.dist, n.avg.misha.50.sr,col="#00b4ffff",type = "l",lwd=2,
     xlim=c(20000,8000),ylim=c(0.705,0.711),main="LA-ICP-MS 50 pt average",
     xlab="Distance (micron) from pulp cavity", ylab="87Sr/86Sr")
abline(v=misha.elements.cl$adj.dist[1058-25],lty=2)
abline(v=misha.elements.cl$adj.dist[2304-25],lty=2)
abline(v=misha.elements.cl$adj.dist[3019-25],lty=2)
plot(misha.raw$dist, misha.raw$X88Sr, type = "l",xlim=c(20000,8000),
     main="Voltage of Mass 88 (Sr)",xlab="Distance (micron) from pulp cavity", ylab="Voltage")
abline(v=misha.elements.cl$adj.dist[1058-25],lty=2)
abline(v=misha.elements.cl$adj.dist[2304-25],lty=2)
abline(v=misha.elements.cl$adj.dist[3019-25],lty=2)
plot(misha.elements.cl$adj.dist, misha.elements.cl$Mass.88,type = "l",
     xlim=c(20000,8000),ylim=c(0,20000),col="red",
     main="Counts of Mass 88 (Sr), Elemental Scan",xlab="Distance (micron) from pulp cavity", ylab="Counts")
abline(v=misha.elements.cl$adj.dist[1058-25],lty=2)
abline(v=misha.elements.cl$adj.dist[2304-25],lty=2)
abline(v=misha.elements.cl$adj.dist[3019-25],lty=2)

#################Fig S2 ################
#comparing solution methods to laser ablation methods
misha.micromill <- read.csv("data/Misha micromill 400-600.csv")

misha.micromill<-na.omit(misha.micromill)

#the micromill was done at 500-micron increment, but a geometric conversion is needed (explained in SI)
#here the conversion is done by multiplying 547 (microns)
misha.micromill.dist <- misha.micromill$Position * 547

misha.micromill.sr <- misha.micromill$corr..87Sr.86Sr
misha.micromill.sr.err <- misha.micromill$comb..Err
par(mfrow=c(1,1))
plot(n.avg.misha.50.dist, n.avg.misha.50.sr,type = "l",col="#00b4ffff",lwd=1.5,
     xlim=c(23000,0),ylim=c(0.705,0.712),xlab="Distance from pulp cavity (microns)",ylab="87Sr/86Sr",
     main="Comparing data from the LA-ICP-MS and the solution methods")
PlotPE(n.avg.misha.50.dist, n.avg.misha.50.sr,n.sd.misha.50.sr,col="#00b4ffff")
points(n.avg.misha.50.dist, n.avg.misha.50.sr, col="#00b4ffff",pch=18,cex=1.1)
PlotPE(misha.micromill.dist,misha.micromill.sr,misha.micromill.sr.err,col="red")
legend(10000, 0.707, c("50-point avg. LA-ICP-MS","Micromill"),
       pch = c(16, 16), col=c("#00b4ffff","red"))

##################Fig S3 sensitivity test Sensitivity to cauchy rate parameter################

par(mfrow=c(3,1))
plot(density(post.misha.pc2p3$BUGSoutput$sims.list$a, from = 0), xlim = c(0.01,0.04),ylim= c(0,120),
     lwd = 2, col = plot.col[2],main="Parameter a, rate = 2e-7", xlab="Parameter estimate")
lines(density(post.misha.invmamm.param$BUGSoutput$sims.list$a.m, from = 0),lwd = 2,col=plot.col[6])
legend(0.03,120, c("Calibration","Mammoth"),lwd = c(2, 2), col=c(plot.col[2],plot.col[6]))

plot(density(post.misha.pc2p3$BUGSoutput$sims.list$a, from = 0), xlim = c(0.01,0.04),ylim= c(0,120),
     lwd = 2, col = plot.col[2],main="Parameter a, rate = 1e-7", xlab="Parameter estimate")
lines(density(post.misha.invmamm.i$BUGSoutput$sims.list$a.m, from = 0),lwd = 2,col=plot.col[6])

plot(density(post.misha.pc2p3$BUGSoutput$sims.list$a, from = 0), xlim = c(0.01,0.04),ylim= c(0,120),
     lwd = 2, col = plot.col[2],main="Parameter a, rate = 3e-8", xlab="Parameter estimate")
lines(density(post.misha.invmamm.s$BUGSoutput$sims.list$a.m, from = 0),lwd = 2,col=plot.col[6])


##################Fig S4 sensitivity test Sensitivity to error structure################
#normal vs cauchy with autocorrelation
par(mfrow=c(2,1))
plot(0,0, xlim = c(1,700), ylim = c(0.705, 0.714), xlab = "days", ylab ="Sr 87/86",
     main = "Fidelity test: Normal error with autocorrelation")
#converting misha distance to days using rate Ivo.rate
points((max(n.avg.misha.50.dist) + 30-n.avg.misha.50.dist.rmv)/14.7,n.avg.misha.50.sr.rmv, pch= 18, col="#00b3ffff")
#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.pc2p3.Rin.89<- MCMC.CI.bound(post.misha.pc2p3$BUGSoutput$sims.list$Rin, 0.89)
lines(1:750,post.misha.pc2p3.Rin.89[[1]],lwd = 2, col = "black")
lines(1:750,post.misha.pc2p3.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:750,post.misha.pc2p3.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

#estimated input series
MCMC.misha.inv2p.tsrw.Rin.m.89 <- MCMC.CI.bound(post.misha.inv2p.tsrw$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:750,MCMC.misha.inv2p.tsrw.Rin.m.89[[1]],lwd = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2p.tsrw.Rin.m.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2p.tsrw.Rin.m.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(400, 0.708, c("Model input","Estimated input"),lwd = c(2, 2), col=c("black","magenta"))

plot(0,0, xlim = c(1,700), ylim = c(0.705, 0.714), xlab = "days", ylab ="Sr 87/86", 
     main = "Fidelity test: Cauchy error with autocorrelation")
#converting misha distance to days using rate Ivo.rate
points((max(n.avg.misha.50.dist) + 30-n.avg.misha.50.dist.rmv)/14.7,n.avg.misha.50.sr.rmv, pch= 18, col="#00b3ffff")
#Use posterior from the calibration run as reference, but use estimates, not posterior
post.misha.pc2p3.Rin.89<- MCMC.CI.bound(post.misha.pc2p3$BUGSoutput$sims.list$Rin, 0.89)
lines(1:750,post.misha.pc2p3.Rin.89[[1]],lwd = 2, col = "black")
lines(1:750,post.misha.pc2p3.Rin.89[[2]], lwd = 1, lty = 2, col = "black")
lines(1:750,post.misha.pc2p3.Rin.89[[3]], lwd = 1, lty = 2, col = "black")

#estimated input series
MCMC.misha.inv2p.tsrwca.Rin.m.89 <- MCMC.CI.bound(post.misha.inv2p.tsrwca$BUGSoutput$sims.list$Rin.m, 0.89)
lines(1:750,MCMC.misha.inv2p.tsrwca.Rin.m.89[[1]],lwd = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2p.tsrwca.Rin.m.89[[2]], lwd = 1, lty = 2, col = "magenta")
lines(1:750,MCMC.misha.inv2p.tsrwca.Rin.m.89[[3]], lwd = 1, lty = 2, col = "magenta")
legend(400, 0.708, c("Model input","Estimated input"),lwd = c(2, 2), col=c("black","magenta"))

################Fig S5: #forward model################
#plot input, simulated ivory with different turnover parameter values
par(mfrow=c(1,2))
plot(0,0, xlim = c(1,900), ylim = c(0.706, 0.711), xlab = "days", ylab ="Sr 87/86",main="Pool ratio")
lines(1:900, input.misha)
lines(1:900, res.misha,lwd=2, col = "#00b4ffff")
lines(1:900, res.prl,lwd=2, col = plot.col.6[1])
lines(1:900, res.prs,lwd=2, col = plot.col.6[3])
legend(300,0.708,c("Intake","Reference (Misha)","Larger pool ratio","Smaller pool ratio"),
       col=c("black","#00b4ffff",plot.col.6[1],plot.col.6[3]),
       lwd = c(1,2,2,2))

plot(0,0, xlim = c(1,900), ylim = c(0.706, 0.711), xlab = "days", ylab ="Sr 87/86",main="Flux ratio")
lines(1:900, input.misha)
lines(1:900, res.misha,lwd=2, col = "#00b4ffff")
lines(1:900, res.frl,lwd=2, col = plot.col.6[1])
lines(1:900, res.frs,lwd=2, col = plot.col.6[3])
legend(300,0.708,c("Intake","Reference (Misha)","Larger flux ratio","Smaller flux ratio"),
       col=c("black","#00b4ffff",plot.col.6[1],plot.col.6[3]),
       lwd = c(1,2,2,2))

#################Fig S6: forward model################
###plots###
par(mfrow=c(3,1))
plot(0,0, xlim = c(1,720), ylim = c(syn.low, syn.high), xlab = "days", ylab ="87Sr/86Sr",
     main="150-day ends & 30-day intermediate switches")
lines(1:length(syn.input.150), syn.input.150)
lines(1:length(syn.input.150), res150[[1]],lwd=2, col = "#00b4ffff")
legend(600,0.711,c("Intake","Ivory"),lwd = c(1,2), col=c("black","#00b4ffff"))
##150 day switch ~65% of the original input

plot(0,0, xlim = c(1,720), ylim = c(syn.low, syn.high), xlab = "days", ylab ="87Sr/86Sr",
     main="120-day ends & 60-day intermediate switches")
lines(1:length(syn.input.120), syn.input.120)
lines(1:length(syn.input.120), res120[[1]],lwd=2, col = "#00b4ffff")
##120 day switch ~60% of the original input

plot(0,0, xlim = c(1,720), ylim = c(syn.low, syn.high), 
     xlab = "days", ylab ="87Sr/86Sr",main="60-day ends & 30-day intermediate switches")
lines(1:length(syn.input.60), syn.input.60)
lines(1:length(res60[[1]]), res60[[1]],lwd=2, col = "#00b4ffff")

##60 day switch ~50% of the original input

###############Fig S7: comparing reaction progress and BITS###########
#reaction progress lines
par(mfrow=c(1,2))
plot(1:length(reac.prog.tk), log(reac.prog.tk),ylab="ln(1-F)",xlab="Days",main="Slow turnover pool")
abline(lm.res.misha.2$coefficients[[1]], lm.res.misha.2$coefficients[[2]], col = "red")
legend(350,0,c("Intercept = -0.654", "Slope = -0.0021","f2 = 0.52"))

plot(c(1:150),reac.prog.res1,ylab="Residual: ln(1-F)",xlab="Days", main = "Fast turnover pool")
abline(lm.res.misha.1$coefficients[[1]], lm.res.misha.1$coefficients[[2]], col = "red")
legend(60,-0.8,c("Intercept = -0.725", "Slope = -0.0338","f1 = 0.48"))
