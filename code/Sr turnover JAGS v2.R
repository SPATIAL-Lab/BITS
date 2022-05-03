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

  Body.mass ~ dnorm(Body.mass.mean, 1/Body.mass.sd^2)
  Body.mass.mean <- 4800 # kg
  Body.mass.sd <- 250 # kg

  flux.ratio <- a/b
  pool.ratio <- c/b
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
  c.coef ~ dunif(0, 1)
  b <- a * b.coef
  b.coef ~ dunif(0, 1)
  a ~ dunif(0, 1)
  # b ~ dunif(0, 1)
  # c ~ dunif(0, 1)

  #model initial values for bone and serum
  
  Rb.m[1] ~ dnorm(R0, Sr.pre.b) #use a different error term
  Rs.m[1] ~ dnorm(R0, Sr.pre.s)
  
  #generate time series of input values
  #Rin is the input ratio, which is modeled with a change point C in the time series
  for (i in 1:t){
    Rin[i] ~ dnorm(Rin.mean[i], Sr.pre.s)# not measurement error here! use a different error term
  }

  #create a vector for the means of Rins
  # Rin.mean <- Rin.mean.0 * inv.trial + Rin.mean.e * trial.sort
  # 
  # Rin.mean.e <- rep(Re, t)
  # 
  # Rin.mean.0 <- rep(R0, t)
  # 
  # #inverse of the vector, also sum of the inv.trial is the position of the change point
  # Cha.p <- sum(inv.trial)
  # inv.trial <- 1 - trial.sort
  # 
  # #sort this vector to create 1 or 0 values to calculate start values
  # trial.sort <- sort(trial)
  # 
  # #change point is modeled with a series of Bernoulli trials
  # for(i in 1:t){
  #   trial[i] ~ dbern(change.p)
  # }
  # 
  # #setting change point propability with an estimate and use a random process to generate intake values
  # change.p <- 0.1


  #precision for average Sr measurements in bone should be much smaller, here represented using rate
  Sr.pre.b ~ dgamma(Sr.pre.shape, Sr.pre.rate.b)  
  Sr.pre.rate.b <- 1e-6
  #precition for average Sr measurements in ivory and serum is estimated
  #based on 25 point average, which centers around 1e-4
  Sr.pre.s ~ dgamma(Sr.pre.shape, Sr.pre.rate.s)  

  Sr.pre.shape <- 100
  Sr.pre.rate.s <- 4e-5
}

# #physiological parameters
# #Fb: Sr flux of bone (mmol) is unknown
# #Fin: Sr flux of intake (mmol) is unknown
# 
# #pool size bone (mmol)
# 
# #Calcium mass is about 30% of bone mass, convert this to mmol, then with Sr/Ca, 20% can be mobolized?
# Pb <- 0.3 * Bone.mass /40 * 1e6 * S.Sr.Ca.ratio #mmol ~6.6k mmol
# 
# #bone mass is about 20% of body mass for an elephant
# Bone.mass <- 0.2 * Body.mass # kg
# 


