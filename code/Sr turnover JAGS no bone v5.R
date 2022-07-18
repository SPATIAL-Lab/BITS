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

  #Data evaluation
  for (i in 1:n.mea){

    #evaluate its ratio
    R.mea[i] ~ dnorm(Rs.eva[i], 1/R.sd.mea[i]^2)
  }
  #problem, 1 in 8 values are wasted
 
  #Evaluating the mean of neighboring 7 data points, ~ 100 micron
  for (i in 2:(n.mea)){
    Rs.eva[i] <- (Rs.m[mod.index[i]] + Rs.m[mod.index[i]+1] + Rs.m[mod.index[i]-1]+
                    Rs.m[mod.index[i]+2] + Rs.m[mod.index[i]-2]+
                    Rs.m[mod.index[i]+3] + Rs.m[mod.index[i]-3])/7
  }
  
  Rs.eva[1] <- Rs.m [1]
  
  #also record modeled distance from pulp cavity
  mod.dist <- 19200 - dist.index * Ivo.rate

  #converting distance values to a set of indexes (integer)
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

  #model initial values for serum assuming that it is in equilibrium with R0.mean
  
  Rs.m[1] ~ dnorm(R0.mean, R0.pre)
  
  #generate time series of input values
  #Rin is the input ratio

  #generate time series of input values
  #Rin is the input ratio, which is modeled with a switch point in the time series
  for(i in 1:t){
    
    Rin[i] ~ dnorm(ifelse(i > switch, Re.mean[i], R0.mean), ifelse(i > switch, Re.pre, R0.pre))
    #When i is greater than the switch point, Rin ~ dnorm(Re.mean, Rin.m.pre)
  }
  
  switch ~ dcat(pi)
  
  #building a vector of weights for the switch point and allow some error
  #the suspected switch point is 76 out of t
  #so the values between 73 and 79 are allowed
  pi <- c(pi.int, pi.switch, pi.end)
  
  pi.end <- rep(0, (t - date - err.date))
  
  pi.switch <- rep(1, (1 + 2 * err.date))
  
  pi.int <- rep(0, (date - err.date - 1))
  
  #uncertainty of the date of switch = +- the number of days
  err.date <- 3 
  
  #suspected date of the switch
  date <- 77
  
  Re.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.Re)
  
  Sr.pre.rate.Re <- 2e-6
  
  R0.mean ~ dnorm(R0, R0.pre) #R0.mean is modeled to be fixed

  R0.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.R0)

  Sr.pre.rate.R0 <- 2e-6
  
  # Sr.pre.s ~ dgamma(Sr.pre.shape, Sr.pre.rate.s)  #serum precision
  # 
  # Sr.pre.rate.s <- 1e-5
  
  Sr.pre.shape <- 100
  
  #generate Sr values of food using a mixing process of food and water
  for(i in 1:t){
    #weighted mean calculated by Sr concentration and ratio
    Re.mean[i] = (sum.c.Sr.pel[i] + sum.c.Sr.hay[i] + sum.c.Sr.alf[i] + c.Sr.w[i] * w.intake[i])/(w.food[i] + w.water[i])
 
    w.food[i] = w.hay[i] + w.alf[i] + w.pel[i] #total amount of Sr in food
    
    #calculate weights based on concentrations, 20kg each for hay and alfalfa
    #with the digestibility modifier for hay and alfalfa and diet ratio
    #total amount of Sr in hay, mass * concentration * digestibility * ratio
    w.hay[i] = sum(c.hay.m) * diges * f.h[i] / m.feed * f.intake[i] *(1 - f.pel[i])
    
    sum.c.Sr.hay[i] = sum(c.Sr.hay) * diges * f.h[i] / m.feed * f.intake[i] *(1 - f.pel[i])
    
    #total amount of Sr in afalfa, mass * concentration * digestibility * ratio
    w.alf[i] = sum(c.alf.m) * diges * (1 - f.h[i])/ m.feed * f.intake[i] *(1 - f.pel[i])
    
    sum.c.Sr.alf[i] = sum(c.Sr.alf) * diges * (1 - f.h[i])/ m.feed * f.intake[i] *(1 - f.pel[i])
    
    #pellet is assumed to be homogenized and completely digestible (Sr)
    #calculate weight of the pellet
    w.pel[i] = c.pel.m * f.pel[i] * f.intake[i] #total amount of Sr in pellet
    
    sum.c.Sr.pel[i] = Sr.pel.m * c.pel.m * f.pel[i] * f.intake[i]
    
    #Water Sr values may change by the day
    #calculate weight of the water
    w.water[i] = w.intake[i] * c.w.m[i] #total amount of Sr in water intake, volume * concentration
    
    c.w.m[i] ~ dlnorm(conc.w.mean, 1/conc.w.sd^2)
    
    Sr.w.m[i] ~ dnorm(Sr.w.mean, 1/Sr.w.sd^2)
    
    c.Sr.w[i] = Sr.w.m[i] * c.w.m[i] #product of concentration and Sr
    
    #pellet is fed at ~10% food by mass(Wood et al., 2020)
    
    f.pel[i] ~ dbeta(alpha.pel, beta.pel) #Pellet, Hay, Alfalfa, uninformative prior
    
    f.h[i] ~ dbeta(alpha.h, beta.h) #Hay, Alfalfa, uninformative prior
    
    f.intake[i] = f.intake.perc[i] * Body.mass
    #daily food intake ~ 1.8% +- 0.2% body mass
    f.intake.perc[i] ~ dnorm(0.018, 1/0.002^2)
    
    w.intake[i] = w.intake.perc[i] * Body.mass
    #daily water intake ~ 1.5% +- 0.2% body mass
    w.intake.perc[i] ~ dnorm(0.015, 1/0.002^2)

  }
  #Define parameters for beta distributions
  alpha.pel ~ dnorm(4, 1/0.4^2)
  
  beta.pel ~ dnorm(40, 1/4^2)
  
  alpha.h ~ dnorm(20, 1/2^2)
  
  beta.h ~ dnorm(20, 1/2^2)
  
  #digestibility of hay and alfalfa should be lower than pellet, which is estimated to be 30 +- 3%
  diges ~ dnorm(0.3, 1/0.03^2)
  
  #pellet concentration and Sr ratio
  c.pel.m ~ dlnorm(conc.pel.mean, 1/conc.pel.sd^2)
  
  Sr.pel.m ~ dnorm(Sr.pel.mean, 1/Sr.pel.sd^2)
  
  #modeling mixing of m.feed kgs of hay and alfalfa in the feeding process
  for(i in 1:m.feed){
    
    #hay
    c.hay.m[i] ~ dlnorm(conc.hay.mean, 1/conc.hay.sd^2)
    
    Sr.hay.m[i] ~ dnorm(Sr.hay.mean, 1/Sr.hay.sd^2)
    
    c.Sr.hay[i] = Sr.hay.m[i] * c.hay.m[i] #product of concentration and Sr
    
    #alfalfa
    c.alf.m[i] ~ dlnorm(conc.alf.mean, 1/conc.alf.sd^2)
    
    Sr.alf.m[i] ~ dnorm(Sr.alf.mean, 1/Sr.alf.sd^2)
    
    c.Sr.alf[i] = Sr.alf.m[i] * c.alf.m[i] #product of concentration and Sr
    
  }
  
  #time series for Sr ratio of feed and pellet
  # for(i in 2:t){
  #   
  #   pel.r[i] <- pel.r[i - 1] + pel.r.cps[i]
  #   
  #   pel.r.cps[i] ~ dnorm(pel.r.cps[i - 1] * pel.r.cps.ac, pel.r.pre)
  #   
  # }
  # pel.r.cps[1] ~ dnorm(0, pel.r.pre)
  # 
  # pel.r.cps.ac ~ dunif(0.01, 0.99)
  # 
  # pel.r.pre ~ dgamma(pel.r.pre.shape, pel.r.pre.rate)
  # 
  # pel.r.pre.shape = 100
  # 
  # pel.r.pre.rate = 1
    
  Body.mass ~ dnorm(Body.mass.mean, 1/Body.mass.sd^2)
  Body.mass.mean <- 4800 # kg
  Body.mass.sd <- 250 # kg
}