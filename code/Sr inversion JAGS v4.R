model {
  ######inversion######
  for (i in 1:n.mea){
    
    #evaluate its ratio
    R.mea[i] ~ dnorm(Rs.eva[i], 1/R.sd.mea[i]^2)
  }
  
  #Evaluating the mean of neighbouring 7 data points, ~ 100 micron
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
  mod.dist <- max.dist - dist.index * Ivo.rate
  
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
    
    #Rin.m.eps[i] ~ dnorm(Rin.m.eps[i - 1] * Rin.m.eps.ac, Rin.m.pre)
    
    #autocorrelation structure of cps (change per step)
    Rin.m.cps[i] ~ dnorm(Rin.m.cps[i - 1] * Rin.m.cps.ac[i], Rin.m.pre)
    
    #autocorrelation term is also a time series with an autocorrelation structure
    #it is centered around the previous step with some variation allowed
    #this is used to accommodate the wide range of autocorrelation values in the actual data
    #including a sharp increase in input values during the switch, and steady values before and after
    Rin.m.cps.ac[i] ~ dnorm(Rin.m.cps.ac[i - 1], Rin.m.cps.ac.pre)
  }
  
  # initiate the series with an uninformative prior
  Rin.m[1] ~ dnorm(Rin.m.int, Rin.m.pre)
  
  Rin.m.int ~ dunif(min.inv, max.inv) #any value between min and max values
  
  #setting upper and lower bounds of the uninformative prior
  min.inv = 0.701
  max.inv = 0.715
  
  #initial change per step centered around 0
  Rin.m.cps[1] ~ dnorm(0, Rin.m.pre)
  
  Rin.m.pre ~ dgamma(Rin.m.pre.shp, Rin.m.pre.rate)
  Rin.m.pre.shp = 100
  Rin.m.pre.rate = 1e-5
  
  Rin.m.cps.ac.pre ~ dgamma(Rin.m.cps.ac.pre.shp, Rin.m.cps.ac.pre.rate)
  
  Rin.m.cps.ac.pre.shp = 20
  Rin.m.cps.ac.pre.rate = 2
  
  #it can also be modeled as a sequence with independent values
  #the initial value is from an uninformative distribution
  Rin.m.cps.ac[1] ~ dunif(0.01, 0.99)

  ####scaling parameters a, b, c to body mass of the subject#####
  #adjusting a, b and c to the body mass of the elephant investigated
  #rate ~ scale with e3/4 body mass (basal matabolic rate)
  #pool ~ scale with 1 body mass
  #a, b and c are rate/pool, so it should scale with -1/4 body mass
  a.m <- a * (Body.mass.m/Body.mass)^-0.25
  b.m <- b * (Body.mass.m/Body.mass)^-0.25
  c.m <- c * (Body.mass.m/Body.mass)^-0.25
  
  #For example, Mammuthus primigenius is estimated to be around 9500 +- 500 kg
  #for the purpose of demonstration, here we use the same parameters as in Misha
  Body.mass.m ~ dnorm(Body.mass.m.mean, 1/Body.mass.m.sd^2)
  Body.mass.m.mean <- 4800 # kg
  Body.mass.m.sd <- 250 # kg
  
  ########calibration process for parameters a, b, and c in misha######
  #Data eveluation
  for (i in 1:n.cal){

    #evaluate its ratio
    #R.sd.cal is from the measurement error
    R.cal[i] ~ dnorm(Rs.cal.eva[i], 1/R.sd.cal[i]^2)
  }
 
  #Evaluating the mean of neighbouring 7 data points, ~ 100 micron
  for (i in 2:n.cal){
    Rs.cal.eva[i] <- (Rs.cal[mod.index.cal[i]] + Rs.cal[mod.index.cal[i]+1] + Rs.cal[mod.index.cal[i]-1]+
                    Rs.cal[mod.index.cal[i]+2] + Rs.cal[mod.index.cal[i]-2]+
                    Rs.cal[mod.index.cal[i]+3] + Rs.cal[mod.index.cal[i]-3])/7
  }
  
  for (i in 2:n.cal){ #it is safe to evaluate the sequence at the second measurement
    for(j in 1:n.days.bef.aft.cal){
      Rs.cal.bef[i, j] <- Rs.cal[mod.index[i] - j]
      Rs.cal.aft[i, j] <- Rs.cal[mod.index[i] + j]
    }
    Rs.cal.eva[i] <- (Rs.cal[mod.index[i]] + sum(Rs.cal.bef[i,]) + sum(Rs.cal.aft[i,]))/(2*n.days.bef.aft.cal + 1)
  }
  
  Rs.cal.eva[1] <- Rs.cal[1]
  
  #also record modeled distance from pulp cavity
  mod.dist.cal <- max.dist.cal - dist.index.cal * Ivo.rate.cal
  
  #converting distance values to a set of indexes (integer)
  #mod.index is measured in days ~1:1050
  mod.index.cal <- round(dist.index.cal) + 1
  
  dist.index.cal <- (max.dist.cal - dist.cal)/Ivo.rate.cal 
  
  max.dist.cal = 19200 #maximum distance from the pulp cavity in microns
  
  #Data model of ivory
  #Priors for ivory sampling
  #assuming laser ablation has no averaging effects

  Ivo.rate.cal ~ dnorm(Ivo.rate.cal.mean, Ivo.rate.cal.pre) #ivory growth rate, micron/day
  Ivo.rate.cal.mean <- 14.7 #microns per day
  Ivo.rate.cal.pre <- 1/0.6^2 # 1 sd = 0.6 according to Uno 2012
  
  #generate time series
  for (i in 2:t.cal){
    #serum ratio
    Rs.cal[i] <- Rs.cal[i - 1] + b * (Rb.cal[i - 1] - Rs.cal[i - 1]) + a * (Rin.cal[i - 1] - Rs.cal[i - 1])

    #bone ratios
    Rb.cal[i] <- Rb.cal[i - 1] + c * (Rs.cal[i - 1] - Rb.cal[i - 1])
  }
  
  #define the three parameters with uninformative priors
  c <- b * c.coef
  c.coef ~ dunif(0.01, 1)
  b <- a * b.coef
  b.coef ~ dunif(0.01, 1)
  a ~ dunif(0.01, 1)

  #model initial values for bone and serum
  
  Rb.cal[1] ~ dnorm(Rs.cal[1], Sr.pre.b) #use a different error term
  Rs.cal[1] ~ dnorm(R0.mean, R0.pre)
  
  #generate time series of input values, using a correlation structure
  #Rin is the input ratio, which is modeled with a switch point in the time series
  
  for(i in 1:t.cal){
    
    Rin[i] ~ dnorm(ifelse(i > switch, Re.mean[i], R0.mean), ifelse(i > switch, Re.pre, R0.pre))
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
  
  
  #allowing some uncertainty in Re and R0 values
  Re.mean ~ dnorm(Re, Re.pre)
  
  #Re has more uncertainty than Re
  Re.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.Re)
  R0.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.R0)
  Sr.pre.rate.Re <- 2e-5
  Sr.pre.rate.R0 <- 2e-6
  
  
  #precision for average Sr measurements in bone should be much smaller, here represented using rate
  Sr.pre.b ~ dgamma(Sr.pre.shape, Sr.pre.rate.b)  
  Sr.pre.rate.b <- 1e-6
  
  Sr.pre.s ~ dgamma(Sr.pre.shape, Sr.pre.rate.s) 
  
  Sr.pre.shape <- 100
  Sr.pre.rate.s <- 2e-5
  
  #generate Sr values of food using a mixing process of food and water
  for(i in 1:t.cal){
    w.contrib[i] = w.water[i]/(w.food[i] + w.water[i])
    
    h.contrib[i] = w.hay[i]/(w.food[i] + w.water[i])
    #weighted mean calculated by Sr concentration and ratio
    Re.mean[i] = (sum.c.Sr.pel[i] + sum.c.Sr.hay[i] + sum.c.Sr.sup[i] + c.Sr.w * w.intake[i])/(w.food[i] + w.water[i])
    
    w.food[i] = w.hay[i] + w.sup[i] + w.pel[i] #total amount of Sr in food
    
    #calculate weights based on concentrations, 20kg each for hay and alfalfa
    #with the digestibility modifier for hay and alfalfa and diet ratio
    #total amount of Sr in hay, mass * concentration * digestibility * ratio
    w.hay[i] = sum(c.hay.m) * diges.h * f.h[i] / m.feed * f.intake[i] 
    
    sum.c.Sr.hay[i] = sum(c.Sr.hay) * diges.h * f.h[i] / m.feed * f.intake[i]
    
    #pellet is assumed to be homogenized and 60% digestible (Sr) (half the fiber content as hay)
    #calculate weight of the pellet
    w.pel[i] = c.pel.m * f.pel[i] * diges.p * p.intake[i] #total amount of Sr in pellet
    
    sum.c.Sr.pel[i] = Sr.pel.m * c.pel.m * f.pel[i] * diges.p * p.intake[i]
    
    #calculate weight of the supplement
    w.sup[i] = c.sup.m * (1- f.pel[i]) * diges.s * p.intake[i] #total amount of Sr in pellet
    
    sum.c.Sr.sup[i] = Sr.sup.m * c.sup.m * (1- f.pel[i]) * diges.s * p.intake[i]
    
    #Water Sr values
    #calculate weight of the water
    w.water[i] = w.intake[i] * c.w.m #total amount of Sr in water intake, volume * concentration
    
    #pellet is fed at ~10% food by mass(Wood et al., 2020)
    
    p.intake[i] = f.intake[i] *(1 - f.h[i])
    
    f.pel[i] ~ dbeta(alpha.pel, beta.pel) #Pellet vs supplement
    
    f.h[i] ~ dbeta(alpha.h, beta.h) #Hay ~80%
    
    f.intake[i] = f.intake.perc[i] * Body.mass
    #daily food intake ~ 1.5% +- 0.2% body mass
    f.intake.perc[i] ~ dnorm(0.015, 1/0.002^2)
    
    w.intake[i] = w.intake.perc[i] * Body.mass
    #daily water intake ~ 1.5% +- 0.2% body mass
    w.intake.perc[i] ~ dnorm(0.015, 1/0.002^2)
    
  }
  #water concentration and Sr ratio
  c.w.m ~ dlnorm(conc.w.mean, 1/conc.w.sd^2)
  
  Sr.w.m ~ dnorm(Sr.w.mean, 1/Sr.w.sd^2)
  
  c.Sr.w = Sr.w.m * c.w.m #product of concentration and Sr
  
  #Define parameters for beta distributions
  alpha.h ~ dnorm(32, 1/0.4^2)
  
  beta.h ~ dnorm(8, 1/4^2)
  
  alpha.pel ~ dnorm(40, 1/1^2)
  
  beta.pel ~ dnorm(40, 1/1^2)
  
  #Ca absorption rate of hay is estimated to be 60 +- 2%
  diges.h ~ dnorm(0.6, 1/0.03^2)
  
  #Ca absorption rate of pellet is assumed to be slightly higher 70% (15% fiber)
  diges.p ~ dnorm(0.7, 1/0.03^2)
  
  #Ca absorption rate of supplement is assumed to be around the same as pellet at 70% (12% fiber)
  diges.s ~ dnorm(0.7, 1/0.03^2)
  
  #supplement concentration and Sr ratio
  c.sup.m ~ dlnorm(conc.sup.mean, 1/conc.sup.sd^2)
  
  Sr.sup.m ~ dnorm(Sr.sup.mean, 1/Sr.sup.sd^2)
  
  #pellet concentration and Sr ratio
  c.pel.m ~ dlnorm(conc.pel.mean, 1/conc.pel.sd^2)
  
  Sr.pel.m ~ dnorm(Sr.pel.mean, 1/Sr.pel.sd^2)
  
  #modeling mixing of m.feed kgs of hay in the feeding process
  for(i in 1:m.feed){
    
    c.hay.m[i] ~ dlnorm(conc.hay.mean, 1/conc.hay.sd^2)
    
    Sr.hay.m[i] ~ dnorm(Sr.hay.mean, 1/Sr.hay.sd^2)
    
    c.Sr.hay[i] = Sr.hay.m[i] * c.hay.m[i] #product of concentration and Sr
    
  }
  
  Body.mass ~ dnorm(Body.mass.mean, 1/Body.mass.sd^2)
  Body.mass.mean <- 4800 # kg
  Body.mass.sd <- 250 # kg
}