
###prep data###
R.sd.mea <- n.sd.misha.25.sr.rmv
dist.mea <- n.avg.misha.25.dist.rmv
R.mea <- n.avg.misha.25.sr.rmv
n.mea = length(n.avg.misha.25.sr.rmv)

R0 <- 0.707

#Re is the mean ratio of end value  
Re <- 0.711
s.intv <- 26.2

Ivo.rate.mean <- 14.7 #microns per day
Ivo.rate.sd <- 0.6
max.dist.mea <- max(n.avg.misha.25.dist) + 30 #adding a bit of distance leading into the series

#parameters to save
parameters <- c("a", "b","c", "Ivo.rate", "Rs.m","Rb.m","Rin","dist.index",
                "Sr.pre", "Ps", "Fin","Pb", "Fb", "Re.mean", "switch","dist",
                "flux.ratio", "pool.ratio","Body.mass")
##Data to pass to the model
dat = list(R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = 400, n.mea = n.mea, 
           R0 = R0, Re = Re, s.intv = s.intv, Ivo.rate.mean = Ivo.rate.mean,
           Ivo.rate.sd = Ivo.rate.sd, max.dist.mea = max.dist.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e3
n.burnin = 4e2
n.thin = floor(n.iter-n.burnin)/400

#Run it
post.misha.fdnbhr = do.call(jags.parallel,list(model.file = "code/Gabe JAGS 2 pool.R", 
                                               parameters.to.save = parameters, 
                                               data = dat, n.chains=5, n.iter = n.iter, 
                                               n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 3.5 hours

post.misha.fdnbhr$BUGSoutput$summary
traplot(post.misha.fdnbhr,parms = c("a"))

##plot calibration curves, run Helper functions first
plot(0,0, xlim = c(20000,13000), ylim = c(0.7066, 0.712), xlab = "distance", ylab ="Sr 87/86",
     main="One-pool model")
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

