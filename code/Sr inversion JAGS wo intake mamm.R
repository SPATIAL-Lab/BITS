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
  
  #generate null input series
  for (i in 2:t){
    
    Rin.m[i] <- Rin.m[i - 1] + Rin.m.cps[i]
    
    Rin.m.cps[i] ~ dt(0, Rin.m.pre, 1) T(-5e-3, 5e-3) #Brownian motion, cauchy error term

  }
  # initiate the series with an reasonable prior
  Rin.m[1] ~ dnorm(Rin.int, Rin.m.pre) #allowed some variation
  
  Rin.int ~ dnorm(0.710, 1/0.01^2)  #a reasonable initial value
  
  #initial change per step
  Rin.m.cps[1] ~ dt(0, Rin.m.pre, 1) T(-5e-3, 5e-3)
  
  Rin.m.pre ~ dgamma(Rin.m.pre.shp, Rin.m.pre.rate)
  Rin.m.pre.shp = 100
  Rin.m.pre.rate = 5e-8

  ####scaling parameters a, b, c to body mass of the subject#####
  #adjusting a, b and c to the body mass of the elephant investigated
  #rate ~ scale with e3/4 body mass (basal matabolic rate)
  #pool ~ scale with 1 body mass
  #a, b and c are rate/pool, so it should scale with -1/4 body mass
  a.m <- a * ((Body.mass.m/Body.mass)^-0.25)
  b.m <- b * ((Body.mass.m/Body.mass)^-0.25)
  c.m <- c * ((Body.mass.m/Body.mass)^-0.25)
  
  #For example, male Mammuthus primigenius is estimated to be around 7500 +- 500 kg
  #for the purpose of demonstration, here we use the same parameters as in Misha
  Body.mass.m ~ dnorm(Body.mass.m.mean, 1/Body.mass.m.sd^2) T(6000, )
  Body.mass.m.mean <- 7500 # kg
  Body.mass.m.sd <- 500 # kg
  
  #body mass of Misha, used to scale parameters a, b and c
  Body.mass ~ dnorm(Body.mass.mean, 1/Body.mass.sd^2) T(2000, )
  Body.mass.mean <- 3000 # kg
  Body.mass.sd <- 150 # kg
  
  ########calibration process for parameters a, b, and c in misha######
  #Data eveluation

  for (i in 1:n.cal){
    #dist[1:t] is a descending sequence
    
    #Evaluating the mean of neighbouring data points
    #this can accommodate variable growth rate of tusk/enamel
    #the following "<" logical operations create vectors with 1s and 0s,
    #dist.mea is used as reference; if dist falls within the bracket, then 1s will be recorded
    #then inprod()/sum() is used to calculate the mean of all data points within the bracket
    
    Rs.cal.eva[i] <- inprod(Rs.cal, ((dist.cal[i] + cal.intv/2) < dist.cal.m) - ((dist.cal[i] - cal.intv/2)  < dist.cal.m))/sum(((dist.cal[i] + cal.intv/2) < dist.cal.m) - ((dist.cal[i] - cal.intv/2)  < dist.cal.m))
    #evaluate its ratio
    R.cal[i] ~ dnorm(Rs.cal.eva[i], 1/R.sd.cal[i]^2)
  }
  
  #Data model priors for ivory growth
  for (i in 2:t.cal){
    Ivo.rate.cal[i] ~ dnorm(Ivo.rate.cal.mean, 1/Ivo.rate.cal.sd^2) #ivory growth rate, micron/day
    dist.cal.m[i] <- dist.cal.m[i - 1] - Ivo.rate.cal[i] #simulate daily distance increment
  }
  
  dist.cal.m[1] <- max.dist.cal#maximum distance from the pulp cavity in microns
  
  Ivo.rate.cal[1] ~ dnorm(Ivo.rate.cal.mean, 1/Ivo.rate.cal.sd^2) #ivory growth rate, micron/day
  
  #generate time series
  for (i in 2:t.cal){
    #serum ratio
    Rs.cal[i] <- Rs.cal[i - 1] + b * (Rb.cal[i - 1] - Rs.cal[i - 1]) + a * (Rin.cal[i - 1] - Rs.cal[i - 1])

    #bone ratios
    Rb.cal[i] <- Rb.cal[i - 1] + c * (Rs.cal[i - 1] - Rb.cal[i - 1])
  }
  
  #define the three parameters with updated priors from the turnover model
  c <- b * c.coef
  c.coef ~ dunif(0.01, 1)
  b <- a * b.coef
  b.coef ~ dunif(0.01, 1)
  a ~ dunif(0, 0.5) #based on one pool model estimates

  #model initial values for bone and serum
  
  Rb.cal[1] ~ dnorm(Rs.cal[1], Sr.pre.b) #use a different error term
  Rs.cal[1] ~ dnorm(R0.mean, R0.pre)
  
  #generate time series of input values, using a correlation structure
  #Rin is the input ratio, which is modeled with a switch point in the time series
  
  for(i in 1:t.cal){
    
    Rin.cal[i] ~ dnorm(ifelse(i > switch, Re.mean, R0.mean), ifelse(i > switch, Re.pre, R0.pre))
    #When i is greater than the switch point, Rin ~ dnorm(Re.mean, Rin.m.pre)
  }
  
  switch ~ dcat(pi)
  
  #building a vector of weights for the switch point and allow some error
  #the suspected switch point is 76 out of t
  #so the values between 74 and 78 are allowed
  pi <- c(pi.int, pi.switch, pi.end)
  
  pi.end <- rep(0, t - date - err.date)
  
  pi.switch <- rep(1, 1 + 2 * err.date)
  
  pi.int <- rep(0, date - err.date - 1)
  
  #uncertainty of the date of switch = +- the number of days
  err.date <- 2 
  
  #suspected date of the switch
  date <- 76
  
  Re.mean ~ dnorm(Re, Re.pre)
  
  #allowing some uncertainty in Re and R0 values
  R0.mean ~ dnorm(R0, R0.pre)
  
  #Re has more uncertainty than Re
  Re.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.Re)
  R0.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.R0)
  Sr.pre.rate.Re <- 2e-5
  Sr.pre.rate.R0 <- 2e-6
  
  Sr.pre.b ~ dgamma(Rin.m.pre.shp, Sr.pre.b.rate)
  Sr.pre.b.rate = 5e-7
  
  Sr.pre.s ~ dgamma(Rin.m.pre.shp, Sr.pre.s.rate)
  Sr.pre.s.rate = 5e-6
  
  Sr.pre.shape <- 100
}