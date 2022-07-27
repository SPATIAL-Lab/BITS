model {
  Fb <- b * Ps 
  Fin <- a * Ps
  Pb <- Ps * b / c
  #pool size serum (mmol)
  Ps <- S.vol * S.Ca * S.Sr.Ca.ratio # mmol ~1 mmol

  #Sr/Ca ratio of serum is assumed to be the same as measured ivory
  S.Sr.Ca.ratio ~ dnorm(S.Sr.Ca.mean, 1/S.Sr.Ca.sd^2)

  S.Sr.Ca.mean <- 0.0011174 # mmol/mmol
  S.Sr.Ca.sd <- 0.000108103 # mmol/mmol

  #Serum Ca concentration
  S.Ca ~ dnorm(S.Ca.mean, 1/S.Ca.sd^2)

  S.Ca.mean <- 2.79 # mmol/L
  S.Ca.sd <- 0.11 # mmol/L

  #assuming serum volume is ~ 7% volume (L) of body mass (kg)
  S.vol <- 0.07 * Body.mass # Liter

  flux.ratio <- a/b
  pool.ratio <- c/b
  #Data eveluation
  for (i in 1:n.mea){

    #evaluate its ratio
    R.mea[i] ~ dnorm(Rs.eva[i], 1/R.sd.mea[i]^2)
  }
  #problem, 1 in 8 values are wasted
 
  #Evaluating the mean of neighbouring data points, ~ 100 micron
  
  for (i in 2:n.mea){
    for(j in 1:n.days.bef.aft){
      Rs.bef[i, j] <- Rs.m[mod.index[i] - j]
      Rs.aft[i, j] <- Rs.m[mod.index[i] + j]
    }
    Rs.eva[i] <- (Rs.m[mod.index[i]] + sum(Rs.bef[i,]) + sum(Rs.aft[i,]))/(2*n.days.bef.aft + 1)
  }
  
  Rs.eva[1] <- Rs.m [1]
  
  #matching measured distance with modeled distance
  for (i in 1:n.mea){
    min.m[i] <- min(abs(dist.mea[i] - dist))
    for(j in 1:t){
      match[i,j] <- ifelse(min.m[i] == abs(dist.mea[i] - dist[j]), j, 0)
    }
    mod.index[i] <- sum(match[i,])
    mod.dist[i] <- dist[mod.index[i]]
  }
  
  for (i in 2:t){
    Ivo.rate[i] ~ dnorm(Ivo.rate.mean, Ivo.rate.pre) #ivory growth rate, micron/day
    dist[i] <- dist[i - 1] - Ivo.rate[i] #cumulative distance
  }
  
  dist[1] <- max.dist
  max.dist <- 19200 #maximum distance from the pulp cavity in microns
  Ivo.rate[1] ~ dnorm(Ivo.rate.mean, Ivo.rate.pre) #ivory growth rate, micron/day
  
  #Data model of ivory
  #Priors for ivory sampling
  #assuming laser ablation has no averaging effects
  Ivo.rate.mean <- 14.7 #microns per day
  Ivo.rate.pre <- 1/0.6^2 # 1 sd = 0.6 according to Uno 2012
  
  #generate time series
  for (i in 2:t){
    #serum ratio
    #derivative of serum ratios is linearly correlated with the difference between serum and the two sources
    #the tow sources are bone and intake
    #Rs.m[i] <- Rs.m[i-1] + Fb/Ps (Rb.m[i-1] - Rs.m[i-1]) + Fin/Ps (Rin[i-1] - Rs.m[i-1])
    Rs.m[i] <- Rs.m[i - 1] + b * (Rb.m[i - 1] - Rs.m[i - 1]) + a * (Rin[i - 1] - Rs.m[i - 1])
    #a/b = Fin/Fb
    #c/b = Ps/Pb
    #also a > b > c
    #bone ratios
    #derivative of bone ratios is linearly correlated with the difference between serum and bone
    #Rb.m[i] <- Rb.m[i-1] + Fb/Pb (Rs.m[i-1] - Rb.m[i-1])
    Rb.m[i] <- Rb.m[i - 1] + c * (Rs.m[i - 1] - Rb.m[i - 1])
  }#can use three parameters instead
  
  #define the three parameters with uninformative priors
  
  c <- b * c.coef
  c.coef ~ dunif(0.01, 1)
  b <- a * b.coef
  b.coef ~ dunif(0.01, 1)
  a ~ dunif(0.01, 1)
  # b ~ dunif(0, 1)
  # c ~ dunif(0, 1)

  #model initial values for bone and serum
  #assume that bone value is similar to serum, but use a different error term
  Rb.m[1] ~ dnorm(Rs.m[1], Sr.pre.b) 
  Rs.m[1] ~ dnorm(R0.mean, R0.pre)
  
  #generate time series of input values
  #Rin is the input ratio, which is modeled with a switch point in the time series
  for(i in 1:t){
    
    Rin[i] ~ dnorm(ifelse(i > switch, Re.mean, R0.mean), ifelse(i > switch, Re.pre, R0.pre))
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

  #allowing some uncertainty in R0 values
  R0.mean ~ dnorm(R0, R0.pre)
  
  Re.mean ~ dnorm(Re, Re.pre)
  
  #Re has more uncertainty than Re
  Re.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.Re)
  R0.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.R0)
  Sr.pre.rate.Re <- 2e-5
  Sr.pre.rate.R0 <- 2e-6

  
  #precision for average Sr measurements in bone should be much smaller
  Sr.pre.b ~ dgamma(Sr.pre.shape, Sr.pre.rate.b)  
  Sr.pre.rate.b <- 2e-6
  
  Sr.pre.s ~ dgamma(Sr.pre.shape, Sr.pre.rate.s) 
  
  Sr.pre.shape <- 100
  Sr.pre.rate.s <- 2e-5
  

  Body.mass ~ dnorm(Body.mass.mean, 1/Body.mass.sd^2)
  Body.mass.mean <- 4800 # kg
  Body.mass.sd <- 250 # kg
  
}



