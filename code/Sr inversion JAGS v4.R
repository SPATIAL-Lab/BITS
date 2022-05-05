model {
  ######inversion######
  for (i in 1:n.mea){
    
    #evaluate its ratio
    R.mea[i] ~ dnorm(Rs.eva[i], 1/R.sd.mea[i]^2)
  }
  
  #Evaluating the mean of neighbouring 7 data points, ~ 100 micron
  #this should be consistent with the averaging pattern of actual data
  #The data we used here is a 100 point average, at 100 micron interval
  for (i in 2:n.mea){
    Rs.eva[i] <- (Rs.m[mod.index[i]] + 
                    Rs.m[mod.index[i]+3] + Rs.m[mod.index[i]-3]+
                    Rs.m[mod.index[i]+2] + Rs.m[mod.index[i]-2]+
                    Rs.m[mod.index[i]+1] + Rs.m[mod.index[i]-1])/7
  }

  
  Rs.eva[1] <- Rs.m [1]
  
  #also record modeled distance from pulp cavity
  mod.dist <- 19200 - dist.index * Ivo.rate
  
  #converting distance values to a set of indeces (interger)
  mod.index <- round(dist.index) + 1
  
  # new ivory distance measurement
  dist.index <- (19200 - dist.mea)/Ivo.rate
  
  #Data model of ivory
  #Priors for ivory sampling
  #assuming laser ablation has no averaging effects, which samples daily Rs values
  
  #Parameters for the ivory growth rate of the subject
  Ivo.rate ~ dnorm(Ivo.rate.mean, Ivo.rate.pre) #ivory growth rate, micron/day
  Ivo.rate.mean <- 14.7 #microns per day
  Ivo.rate.pre <- 1/0.6^2 # 1 sd = 0.6 according to Uno 2012

  for (i in 2:t){
    #serum ratio
    Rs.m[i] <- Rs.m[i - 1] + b.m * (Rb.m[i - 1] - Rs.m[i - 1]) + a.m * (Rin.m[i - 1] - Rs.m[i - 1])
    
    #bone ratios
    Rb.m[i] <- Rb.m[i - 1] + c.m * (Rs.m[i - 1] - Rb.m[i - 1])
  }
  
  # assuming the starting value of the two pools are close to the starting Rin value
  #e.g., close to equilibrium with Rin
  Rb.m[1] ~ dnorm(Rin.m[1], Sr.pre.b) #use a different error term for bone
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
  
  Rin.m.int ~ dunif(R0, Re) #any value between R0 and Re
  
  #initial change per step
  Rin.m.cps[1] ~ dnorm(0, Rin.m.pre)
  
  #parameters for auto-correlation of change per step
  #this is assumed to be a universal parameter of the entire sequence
  #Rin.m.eps.ac ~ dunif(0, 2)
  
  Rin.m.cps.ac.pre ~ dgamma(Rin.m.cps.ac.pre.shp, Rin.m.cps.ac.pre.rate)
  
  Rin.m.cps.ac.pre.shp = 20
  Rin.m.cps.ac.pre.rate = 2
  
  #try a less variable conbination
  # Rin.m.cps.ac.pre.shp = 10
  # Rin.m.cps.ac.pre.rate = 0.2
  
  #it can also be modeled as a sequence with independent values
  #the initial value is from an uninformative distribution
  Rin.m.cps.ac[1] ~ dunif(0, 2)

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
  
  Body.mass ~ dnorm(Body.mass.mean, 1/Body.mass.sd^2)
  Body.mass.mean <- 4800 # kg
  Body.mass.sd <- 250 # kg
  
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
  
  Rs.cal.eva[1] <- Rs.cal [1]
  
  #also record modeled distance from pulp cavity
  mod.dist.cal <- 19200 - dist.index.cal * Ivo.rate.cal

  #converting distance values to a set of indeces (interger)
  mod.index.cal <- round(dist.index.cal) + 1
  
  dist.index.cal <- (19200 - dist.cal)/Ivo.rate.cal
  
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
  c.coef ~ dunif(0, 1)
  b <- a * b.coef
  b.coef ~ dunif(0, 1)
  a ~ dunif(0, 1)

  #model initial values for bone and serum
  
  Rb.cal[1] ~ dnorm(R0.mean, Sr.pre.b) #use a different error term
  Rs.cal[1] ~ dnorm(R0.mean, Sr.pre.s)
  
  #generate time series of input values, using a correlation structure
  #Rin is the input ratio, which is modeled with a switch point in the time series
  
  for(i in 1:t){
    
    Rin.cal[i] ~ dnorm(ifelse(i > switch, Re.mean, R0.mean), Rin.m.pre)
    #When i is greater than the switch point, Rin ~ dnorm(Re.mean, Rin.m.pre)
  }
  
  switch ~ dcat(pi)
  
  #building a vector of weights for the switch point and allow some error
  #the suspected switch point is 75 out of t
  #so the values between 73 and 77 are allowed
  pi <- c(pi.int, pi.switch, pi.end)
  
  pi.end <- rep(0, t - date - err.date)
  
  pi.switch <- rep(1, 1 + 2 * err.date)
  
  pi.int <- rep(0, date - err.date - 1)
  
  #uncertainty of the date of switch = +- the number of days
  err.date <- 3 
  
  #suspected date of the switch
  date <- 75
  
  #allowing some uncertainty in Re and R0 values
  Re.mean ~ dnorm(Re, Sr.pre.s)
  
  R0.mean ~ dnorm(R0, Sr.pre.s)

  #precision for average Sr measurements in bone should be much smaller, here represented using rate
  Sr.pre.b ~ dgamma(Sr.pre.shape, Sr.pre.rate.b)  
  Sr.pre.rate.b <- 1e-6

  #error of the input values, has to be a reasonable value, 1sd at 1e-5 for water measurements
  Rin.m.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.Rin)
  Sr.pre.rate.Rin <- 1e-6
  
  #precition for average Sr measurements in ivory and serum is estimated
  #based on 25 point average, which centers around 1e-4
  Sr.pre.s ~ dgamma(Sr.pre.shape, Sr.pre.rate.s)
  
  Sr.pre.shape <- 100
  Sr.pre.rate.s <- 4e-5
}