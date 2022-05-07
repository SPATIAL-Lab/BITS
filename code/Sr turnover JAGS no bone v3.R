model {

  Fin <- a * Ps

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

  Body.mass ~ dnorm(Body.mass.mean, 1/Body.mass.sd^2)
  Body.mass.mean <- 4800 # kg
  Body.mass.sd <- 250 # kg

  #Data eveluation
  for (i in 1:n.mea){

    #evaluate its ratio
    R.mea[i] ~ dnorm(Rs.eva[i], 1/R.sd.mea[i]^2)
  }
  #problem, 1 in 8 values are wasted
 
  #Evaluating the mean of neighbouring 7 data points, ~ 100 micron
  for (i in 2:(n.mea)){
    Rs.eva[i] <- (Rs.m[mod.index[i]] + Rs.m[mod.index[i]+1] + Rs.m[mod.index[i]-1]+
                    Rs.m[mod.index[i]+2] + Rs.m[mod.index[i]-2]+
                    Rs.m[mod.index[i]+3] + Rs.m[mod.index[i]-3])/7
  }
  
  Rs.eva[1] <- Rs.m [1]
  
  #also record modeled distance from pulp cavity
  mod.dist <- 19200 - dist.index * Ivo.rate

  #converting distance values to a set of indeces (interger)
  mod.index <- round(dist.index) + 1
  
  dist.index <- (19200 - dist.mea)/Ivo.rate
  
  #Data model of ivory
  #Priors for ivory sampling
  #assuming laser ablation has no averaging effects

  Ivo.rate ~ dnorm(Ivo.rate.mean, Ivo.rate.pre) #ivory growth rate, micron/day
  Ivo.rate.mean <- 14.7 #microns per day
  Ivo.rate.pre <- 1/0.6^2 # 1 sd = 0.6 according to Uno 2012
  
  #generate time series
  for (i in 2:t){
    #serum ratio
    #when bone is not a concern, the equation is the same as exponential decay
    Rs.m[i] <- Rs.m[i - 1] + a * (Rin[i - 1] - Rs.m[i - 1])
  }
  
  #define the parameter with uninformative priors
  
  a ~ dunif(0, 1)

  #model initial values for serum
  
  Rs.m[1] ~ dnorm(R0, Sr.pre.s)
  
  #generate time series of input values
  #Rin is the input ratio

  #generate time series of input values
  #Rin is the input ratio, which is modeled with a switch point in the time series
  for(i in 1:t){
    
    Rin[i] ~ dnorm(ifelse(i > switch, Re.mean, R0.mean), ifelse(i > switch, Re.pre, R0.pre))
    #When i is greater than the switch point, Rin ~ dnorm(Re.mean, Rin.m.pre)
  }
  
  switch ~ dcat(pi)
  
  #building a vector of weights for the switch point and allow some error
  #the suspected switch point is 76 out of t
  #so the values between 73 and 79 are allowed
  pi <- c(pi.int, pi.switch, pi.end)
  
  pi.end <- rep(0, t - date - err.date)
  
  pi.switch <- rep(1, 1 + 2 * err.date)
  
  pi.int <- rep(0, date - err.date - 1)
  
  #uncertainty of the date of switch = +- the number of days
  err.date <- 3 
  
  #suspected date of the switch
  date <- 77
  
  
  #allowing some uncertainty in Re and R0 values
  Re.mean ~ dnorm(Re, Re.pre)
  
  R0.mean ~ dnorm(R0, R0.pre)
  
  #Re has more uncertainty than Re
  Re.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.Re)
  R0.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.R0)
  Sr.pre.rate.Re <- 2e-5
  Sr.pre.rate.R0 <- 2e-6
  
  Sr.pre.s ~ dgamma(Sr.pre.shape, Sr.pre.rate.s)  
  
  Sr.pre.shape <- 100
  Sr.pre.rate.s <- 4e-5
}