model {
  ######inversion######

  for (i in 1:n.mea){
    #dist[1:t] is a descending sequence
    
    #Evaluating the mean of neighbouring data points
    #this can accommodate variable growth rate of tusk/enamel
    #the following ">" logical operations create vectors with 1s and 0s,
    #dist.mea is used as reference; if dist falls within the bracket, then 1s will be recorded
    #then inprod()/sum() is used to calculate the mean of all data points within the bracket
    
    Rs.eva[i] <- inprod(Rs.m, ((dist.mea[i] + s.intv/2) < dist) - ((dist.mea[i] - s.intv/2)  < dist))/sum(((dist.mea[i] + s.intv/2) < dist) - ((dist.mea[i] - s.intv/2)  < dist))
    #evaluate its ratio
    R.mea[i] ~ dnorm(Rs.eva[i], 1/R.sd.mea[i]^2)
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
  Rb.m[1] ~ dnorm(Rin.m[1], Sr.pre.b) T(0.700, 0.720)#use a different error term for bone
  Rs.m[1] ~ dnorm(Rin.m[1], Sr.pre.s) T(0.700, 0.720)
  
  #generate null input time series
  for (i in 2:t){
    
    Rin.m[i] <- Rin.m[i - 1] + Rin.m.cps[i] 
    
    #Rin.m.cps[i] ~ dnorm(0, Rin.m.pre) #Brownian motion, gaussian error term
    
    Rin.m.cps[i] ~ dt(0, Rin.m.pre, 1) T(-5e-3, 5e-3) #Brownian motion, cauchy error term
  }
  
  # initiate the series with an reasonable prior
  Rin.m[1] ~ dnorm(Rin.int, Rin.m.pre) #allowed some variation
  
  Rin.int ~ dnorm(0.709, 1e4)  #a reasonable initial value
  #Rin.int ~ dunif(0.700, 0.800)
  
  #initial change per step centered around 0
  Rin.m.cps[1] ~ dt(0, Rin.m.pre, 1) T(-5e-3, 5e-3)
  
  #Rin.m.cps[1] ~ dnorm(0, Rin.m.pre)
  
  Sr.pre.b ~ dgamma(Rin.m.pre.shp, Sr.pre.b.rate)
  Sr.pre.b.rate = 5e-7
  
  Sr.pre.s ~ dgamma(Rin.m.pre.shp, Sr.pre.s.rate)
  Sr.pre.s.rate = 5e-6
  
  Rin.m.pre ~ dgamma(Rin.m.pre.shp, Rin.m.pre.rate)
  Rin.m.pre.shp = 100
  Rin.m.pre.rate = 5e-8

  ####scaling parameters a, b, c to body mass of the subject#####
  # adjusting a, b and c to the body mass of the elephant investigated
  # rate ~ scale with e3/4 body mass (basal matabolic rate)
  # pool ~ scale with 1 body mass
  # a, b and c are rate/pool, so they should scale with -1/4 body mass
  a.m <- a * ((Body.mass.ratio)^-0.25)
  b.m <- b * ((Body.mass.ratio)^-0.25)
  c.m <- c * ((Body.mass.ratio)^-0.25)

  #For example, male Mammuthus primigenius is estimated to be around 7500 +- 500 kg
  #for the purpose of demonstration, here we use the same parameters as in Misha
  Body.mass.ratio <- Body.mass.m/Body.mass
  Body.mass.m ~ dnorm(Body.mass.m.mean, 1/Body.mass.m.sd^2) T(6000, )
  Body.mass.m.mean <- 7500 # kg
  Body.mass.m.sd <- 500 # kg

  #body mass of Misha, used to scale parameters a, b and c
  Body.mass ~ dnorm(Body.mass.mean, 1/Body.mass.sd^2) T(2000, )
  Body.mass.mean <- 3000 # kg
  Body.mass.sd <- 150 # kg
  
  a <- exp(params[1]) 
  
  b <- exp(params[2])
  
  c <- exp(params[3])
  
  #supply the model with estimated parameters
  params ~ dmnorm.vcov(params.mu, params.vcov)
  
}