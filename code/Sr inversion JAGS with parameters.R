model {
  ######inversion######
  for (i in 1:n.mea){
    
    #evaluate its ratio
    R.mea[i] ~ dnorm(Rs.eva[i], 1/R.sd.mea[i]^2)
  }
  #this should be consistent with the averaging pattern of actual data
  #The data we used here is a 100 point average, at 100 micron interval
  #Evaluating the mean of x micron width of ivory sampled
  for (i in 2:n.mea){ #it is safe to evaluate the sequence at the second measurement
    for(j in 1:n.days.bef.aft){
      Rs.bef[i, j] <- Rs.m[mod.index[i] - j]
      Rs.aft[i, j] <- Rs.m[mod.index[i] + j]
    }
    Rs.eva[i] <- (Rs.m[mod.index[i]] + sum(Rs.bef[i,]) + sum(Rs.aft[i,]))/(2*n.days.bef.aft + 1)
  }
  
  
  Rs.eva[1] <- Rs.m [1]
  
  #also record modeled distance from pulp cavity
  mod.dist <- max.dist.mic - dist.index * Ivo.rate
  
  #converting distance values to a set of indexes (integer)
  #mod.index is measured in days 
  mod.index <- round(dist.index) + 1
  
  dist.index <- (max.dist.mic - dist.mea)/Ivo.rate 
  
  max.dist.mic = 8000 #maximum distance from the pulp cavity in microns
  #Data model of ivory
  #Priors for ivory sampling
  #assuming laser ablation has no averaging effects, which samples daily Rs values
  
  #Parameters for the ivory growth rate of the subject
  Ivo.rate ~ dnorm(Ivo.rate.mean, Ivo.rate.pre) #ivory growth rate, micron/day
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
  Rb.m[1] ~ dnorm(Rs.m[1], Sr.pre.b) #use a different error term for bone
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
  Rin.m.cps.ac ~ dunif(0.01, 0.99)
  
  # initiate the series with an uninformative prior
  Rin.m[1] ~ dnorm(Rin.m.int, Rin.m.pre)
  
  Rin.m.int ~ dunif(min.inv, max.inv) #any value between min and max values
  
  #setting upper and lower bounds of the uninformative prior
  min.inv = 0.707
  max.inv = 0.713
  
  #initial change per step centered around 0
  Rin.m.cps[1] ~ dnorm(0, Rin.m.pre)
  
  Sr.pre.b ~ dgamma(Rin.m.pre.shp, Sr.pre.b.rate)
  Sr.pre.b.rate = 1e-6
  
  Sr.pre.s ~ dgamma(Rin.m.pre.shp, Sr.pre.s.rate)
  Sr.pre.s.rate = 2e-5
  
  Rin.m.pre ~ dgamma(Rin.m.pre.shp, Rin.m.pre.rate)
  Rin.m.pre.shp = 100
  Rin.m.pre.rate = 4e-6
  
  Rin.m.cps.ac.pre ~ dgamma(Rin.m.cps.ac.pre.shp, Rin.m.cps.ac.pre.rate)
  
  Rin.m.cps.ac.pre.shp = 20
  Rin.m.cps.ac.pre.rate = 2
  

  
  ####scaling parameters a, b, c to body mass of the subject#####
  #adjusting a, b and c to the body mass of the elephant investigated
  #rate ~ scale with e3/4 body mass (basal matabolic rate)
  #pool ~ scale with 1 body mass
  #a, b and c are rate/pool, so it should scale with -1/4 body mass
  a.m <- a * ((Body.mass.m/Body.mass)^-0.25)
  b.m <- b * ((Body.mass.m/Body.mass)^-0.25)
  c.m <- c * ((Body.mass.m/Body.mass)^-0.25)
  
  #For example, Mammuthus primigenius is estimated to be around 9500 +- 500 kg
  #for the purpose of demonstration, here we use the same parameters as in Misha
  Body.mass.m ~ dnorm(Body.mass.m.mean, 1/Body.mass.m.sd^2)
  Body.mass.m.mean <- 4800 # kg
  Body.mass.m.sd <- 250 # kg
  
  #body mass of Misha, used to scale parameters a, b and c
  Body.mass ~ dnorm(Body.mass.mean, 1/Body.mass.sd^2)
  Body.mass.mean <- 4800 # kg
  Body.mass.sd <- 250 # kg
  
  #supply the model with estimated parameters
  a ~ dlnorm(a.mean, 1/a.sd^2)
  
  b ~ dlnorm(b.mean, 1/b.sd^2)
  
  c ~ dlnorm(c.mean, 1/c.sd^2)

}