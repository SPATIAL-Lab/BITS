model {
  ######inversion######

  for (i in 1:n.mea){
    #dist[1:t] is a descending sequence, no negative value is allowed
    
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
    Ivo.rate[i] ~ dnorm(Ivo.rate.mean, Ivo.rate.pre) #ivory growth rate, micron/day
    dist[i] <- dist[i - 1] - Ivo.rate[i] #simulate daily distance increment
  }
  
  dist[1] <- 8000#maximum distance from the pulp cavity in microns
  
  Ivo.rate[1] ~ dnorm(Ivo.rate.mean, Ivo.rate.pre) #ivory growth rate, micron/day
  
  #Parameters for the ivory growth rate of the subject
  Ivo.rate.mean <- 20 #microns per day
  Ivo.rate.pre <- 1/0.6^2 # 1 sd = 0.6 according to Uno 2012
  
  for (i in 2:t){
    #serum ratio
    Rs.m[i] <- Rs.m[i - 1] + b.m * (Rb.m[i - 1] - Rs.m[i - 1]) + a.m * (Rin.m[i - 1] - Rs.m[i - 1])
    
    #bone ratios
    Rb.m[i] <- Rb.m[i - 1] + c.m * (Rs.m[i - 1] - Rb.m[i - 1])
  }
  
  # assuming the starting value of the two pools are close to the starting Rin value
  #e.g., close to equilibrium with Rin
  Rb.m[1] ~ dnorm(Rs.m[1], Sr.pre.b) T(0.706, 0.713)#use a different error term for bone
  Rs.m[1] ~ dnorm(Rin.m[1], Sr.pre.s)
  
  #generate null input series
  for (i in 2:t){
    
    Rin.m[i] <- Rin.m[i - 1] + Rin.m.cps[i]
    
    Rin.m.cps[i] ~ dnorm(Rin.m.cps[i - 1] * Rin.m.cps.ac, Rin.m.pre)
    
    #autocorrelation structure of cps (change per step)
    #Rin.m.cps[i] ~ dnorm(Rin.m.cps[i - 1] * Rin.m.cps.ac[i], Rin.m.pre)
    
    #autocorrelation term is also a time series with an autocorrelation structure
    #it is centered around the previous step with some variation allowed
    #this is used to accommodate the wide range of autocorrelation values in the actual data
    #including a sharp increase in input values during the switch, and steady values before and after
    # Rin.m.cps.ac[i] ~ dnorm(Rin.m.cps.ac[i - 1], Rin.m.cps.ac.pre)
  }
  
  #it can also be modeled as a sequence with independent values
  #the initial value is from an uninformative distribution
  Rin.m.cps.ac ~ dunif(0.01, 1)
  
  # initiate the series with an reasonable prior
  Rin.m[1] ~ dnorm(0.710, 1e10)
  
  #initial change per step centered around 0
  Rin.m.cps[1] ~ dnorm(0, Rin.m.pre)
  
  Sr.pre.b ~ dgamma(Rin.m.pre.shp, Sr.pre.b.rate)
  Sr.pre.b.rate = 1e-8
  
  Sr.pre.s ~ dgamma(Rin.m.pre.shp, Sr.pre.s.rate)
  Sr.pre.s.rate = 1e-8
  
  Rin.m.pre ~ dgamma(Rin.m.pre.shp, Rin.m.pre.rate)
  Rin.m.pre.shp = 100
  Rin.m.pre.rate = 1e-14
  
  # Rin.m.cps.ac.pre ~ dgamma(Rin.m.cps.ac.pre.shp, Rin.m.cps.ac.pre.rate)
  # 
  # Rin.m.cps.ac.pre.shp = 20
  # Rin.m.cps.ac.pre.rate = 0.1
  

  ####scaling parameters a, b, c to body mass of the subject#####
  #adjusting a, b and c to the body mass of the elephant investigated
  #rate ~ scale with e3/4 body mass (basal matabolic rate)
  #pool ~ scale with 1 body mass
  #a, b and c are rate/pool, so it should scale with -1/4 body mass
  # a.m <- a * ((Body.mass.m/Body.mass)^-0.25)
  # b.m <- b * ((Body.mass.m/Body.mass)^-0.25)
  # c.m <- c * ((Body.mass.m/Body.mass)^-0.25)
  # 
  # #For example, Mammuthus primigenius is estimated to be around 9500 +- 500 kg
  # #for the purpose of demonstration, here we use the same parameters as in Misha
  # Body.mass.m ~ dnorm(Body.mass.m.mean, 1/Body.mass.m.sd^2)
  # Body.mass.m.mean <- 4800 # kg
  # Body.mass.m.sd <- 150 # kg
  # 
  # #body mass of Misha, used to scale parameters a, b and c
  # Body.mass ~ dnorm(Body.mass.mean, 1/Body.mass.sd^2)
  # Body.mass.mean <- 4800 # kg
  # Body.mass.sd <- 150 # kg
  # 
  #perhaps use body mass ratio?
  
  a.m <- exp(params[1])
    
  b.m <- exp(params[2])
    
  c.m <- exp(params[3])
  
  #supply the model with estimated parameters
  params ~ dmnorm.vcov(params.mu, params.vcov)


}