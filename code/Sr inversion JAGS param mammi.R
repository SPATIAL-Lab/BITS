model {
  ######inversion######
  for (i in 1:n.mea){
    
    R1.eva[i] <- inprod(R1.m, ((dist.mea[i] + s.intv/2) < dist) - ((dist.mea[i] - s.intv/2)  < dist))/sum(((dist.mea[i] + s.intv/2) < dist) - ((dist.mea[i] - s.intv/2)  < dist))
    #evaluate its ratio
    R.mea[i] ~ dnorm(R1.eva[i], 1/(R.sd.mea[i])^2)
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
    R1.m[i] <- R1.m[i - 1] + b.m * (R2.m[i - 1] - R1.m[i - 1]) + a.m * (Rin.m[i - 1] - R1.m[i - 1])
    
    #bone ratios
    R2.m[i] <- R2.m[i - 1] + c.m * (R1.m[i - 1] - R2.m[i - 1])
  }
  
  # assuming the starting value of the two pools are close to the starting Rin value
  #e.g., close to equilibrium with Rin
  R2.m[1] ~ dnorm(Rin.m[1], Sr.pre.2) #T(0.700, 0.720)#use a different error term for bone
  R1.m[1] ~ dnorm(Rin.m[1], Sr.pre.1)
  
  #precision for average Sr measurements in bone should be much smaller
  Sr.pre.2 ~ dgamma(Sr.pre.shape, Sr.pre.rate.2)  
  Sr.pre.rate.2 <- 5e-7
  
  Sr.pre.1 ~ dgamma(Sr.pre.shape, Sr.pre.rate.1) 
  
  Sr.pre.shape <- 100
  Sr.pre.rate.1 <- 5e-6
  
  #generate null input series
  for (i in 2:t){
    
    Rin.m[i] <- Rin.m[i - 1] + Rin.m.cps[i]
    
    Rin.m.cps[i] ~ dt(0, Rin.m.pre, 1) T(-6e-3, 6e-3) #Brownian motion, cauchy error term

  }
  # initiate the series with an reasonable prior
  Rin.m[1] ~ dnorm(Rin.int, 1/0.0005^2) #allowed some variation
  
  Rin.int ~ dnorm(0.710, 1/0.01^2)  #an uninformative initial value
  
  #initial change per step, with the maximum difference crossing a discrete Sr isotope boundary
  Rin.m.cps[1] ~ dt(0, Rin.m.pre, 1) T(-6e-3, 6e-3)
  
  Rin.m.pre ~ dgamma(Rin.m.pre.shp, Rin.m.pre.rate)
  Rin.m.pre.shp = 100
  Rin.m.pre.rate = 1e-7
  
  ####scaling parameters a, b, c to body mass of the subject#####
  a.m <- a * ((Body.mass.m/Body.mass)^exp.a)
  b.m <- b * ((Body.mass.m/Body.mass)^exp.bc)
  c.m <- c * ((Body.mass.m/Body.mass)^exp.bc)
  
  #calculation of exponents 
  #fast turnover pool is not affected by body mass (Thomas & Crowther 2015)
  exp.a ~ dnorm(0, 1/0.05^2)

  exp.bc ~ dnorm(-0.19, 1/0.05^2) #whole body turnover allometry from Thomas & Crowther 2015
  
  a <- a.post[indx]
  b <- b.post[indx]
  c <- c.post[indx]
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