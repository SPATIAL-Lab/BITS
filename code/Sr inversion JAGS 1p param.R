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
    Rs.m[i] <- Rs.m[i - 1] + a.m * (Rin.m[i - 1] - Rs.m[i - 1])
    
  }
  
  # assuming the starting value of the two pools are close to the starting Rin value
  #e.g., close to equilibrium with Rin
  Rs.m[1] ~ dnorm(Rin.m[1], Rin.m.pre)
  
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

  a.m <- exp(a)
  
  a ~ dnorm(a.mean, 1/a.sd^2)
    
}