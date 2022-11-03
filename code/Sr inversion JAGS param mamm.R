model {
  ######inversion######
  for (i in 1:n.mea){
    
    Rs.eva[i] <- inprod(Rs.m, ((dist.mea[i] + s.intv/2) < dist) - ((dist.mea[i] - s.intv/2)  < dist))/sum(((dist.mea[i] + s.intv/2) < dist) - ((dist.mea[i] - s.intv/2)  < dist))
    #evaluate its ratio
    R.mea[i] ~ dnorm(Rs.eva[i], 1/(R.sd.mea[i])^2)
  }
  
  #Data model priors for ivory growth
  for (i in 2:t){
    Ivo.rate[i] ~ dnorm(Ivo.rate.mean, 1/Ivo.rate.sd^2) #ivory growth rate, micron/day
    dist[i] <- dist[i - 1] - Ivo.rate[i] #simulate daily distance increment
  }
  
  dist[1] <- max.dist.mea#maximum distance from the pulp cavity in microns
  
  Ivo.rate[1] ~ dnorm(Ivo.rate.mean, 1/Ivo.rate.sd^2) #ivory growth rate, micron/day

  for (i in 2:t){
    #serum ratio
    Rs.m[i] <- Rs.m[i - 1] + b.m * (Rb.m[i - 1] - Rs.m[i - 1]) + a.m * (Rin.m[i - 1] - Rs.m[i - 1])
    
    #bone ratios
    Rb.m[i] <- Rb.m[i - 1] + c.m * (Rs.m[i - 1] - Rb.m[i - 1])
  }
  
  # assuming the starting value of the two pools are close to the starting Rin value
  #e.g., close to equilibrium with Rin
  Rb.m[1] ~ dnorm(Rin.m[1], Sr.pre.b) #T(0.700, 0.720)#use a different error term for bone
  Rs.m[1] ~ dnorm(Rin.m[1], Sr.pre.s)
  
  #precision for average Sr measurements in bone should be much smaller
  Sr.pre.b ~ dgamma(Sr.pre.shape, Sr.pre.rate.b)  
  Sr.pre.rate.b <- 5e-7
  
  Sr.pre.s ~ dgamma(Sr.pre.shape, Sr.pre.rate.s) 
  
  Sr.pre.shape <- 100
  Sr.pre.rate.s <- 5e-6
  
  #generate null input series
  for (i in 2:t){
    
    Rin.m[i] <- Rin.m[i - 1] + Rin.m.cps[i]
    
    Rin.m.cps[i] ~ dt(0, Rin.m.pre, 1) T(-5e-3, 5e-3) #Brownian motion, cauchy error term

  }
  # initiate the series with an reasonable prior
  Rin.m[1] ~ dnorm(Rin.int, Rin.m.pre) #allowed some variation
  
  Rin.int ~ dnorm(0.710, 1/0.01^2)  #an uninformative initial value
  
  #initial change per step
  Rin.m.cps[1] ~ dt(0, Rin.m.pre, 1) T(-5e-3, 5e-3)
  
  Rin.m.pre ~ dgamma(Rin.m.pre.shp, Rin.m.pre.rate)
  Rin.m.pre.shp = 100
  Rin.m.pre.rate = 5e-8
  
  ####scaling parameters a, b, c to body mass of the subject#####
  a.m <- a.post[indx]* ((Body.mass.m/Body.mass)^exp.ab)
  b.m <- b.post[indx]* ((Body.mass.m/Body.mass)^exp.ab)
  c.m <- c.post[indx]* ((Body.mass.m/Body.mass)^exp.c)
  
  #calculation of exponents 
  exp.c <- exp.r.in - exp.bnbm
  exp.ab <- exp.r.in - exp.smbm
  
  #adjusting a, b and c to the body mass of the elephant investigated
  #rates scale with ~e 3/4 body mass (basal matabolic rate)
  #pools scale with ~e 1 body mas, but bone has a slightly higher slope
  #a, b and c are rate/pool, so it should scale with ~ -1/4 body mass
  exp.r.in ~ dnorm(0.75, 1/0.05^2) #intake rates scales with BM allometrically
  exp.bnbm ~ dnorm(1.05, 1/0.05^2) #bone scales with BM allometrically (Christiansen 2002)
  exp.smbm ~ dnorm(1, 1/0.05^2) #serum scales with BM isometrically
  
  #sampling from parameter posteriors from the calibration
  indx ~ dcat(rep(1, post.leng))
  
  #For example, male Mammuthus primigenius is estimated to be around 7500 +- 500 kg
  #for the purpose of demonstration, here we use the same parameters as in Misha
  Body.mass.m ~ dnorm(Body.mass.m.mean, 1/Body.mass.m.sd^2) T(6000, )
  Body.mass.m.mean <- 7500 # kg
  Body.mass.m.sd <- 500 # kg
  
  #body mass of Misha, used to scale parameters a, b and c
  Body.mass ~ dnorm(Body.mass.mean, 1/Body.mass.sd^2) T(2000, )
  Body.mass.mean <- 3000 # kg
  Body.mass.sd <- 150 # kg
  
}