model {

  flux.ratio <- a/b
  pool.ratio <- c/b
  #Data eveluation

  for (i in 1:n.mea){

    R1.eva[i] <- inprod(R1.m, ((dist.mea[i] + s.intv/2) < dist) - ((dist.mea[i] - s.intv/2)  < dist))/sum(((dist.mea[i] + s.intv/2) < dist) - ((dist.mea[i] - s.intv/2)  < dist))
    #evaluate its ratio
    R.mea[i] ~ dnorm(R1.eva[i], 1/R.sd.mea[i]^2)
  }
  
  #Data model priors for ivory growth
  for (i in 2:t){
    Ivo.rate[i] ~ dnorm(Ivo.rate.mean, 1/Ivo.rate.sd^2) #ivory growth rate, micron/day
    dist[i] <- dist[i - 1] - Ivo.rate[i] #simulate daily distance increment
  }
  
  dist[1] <- max.dist.mea#maximum distance from the pulp cavity in microns
  
  Ivo.rate[1] ~ dnorm(Ivo.rate.mean, 1/Ivo.rate.sd^2) #ivory growth rate, micron/day

  #generate time series
  for (i in 2:t){

    R1.m[i] <- R1.m[i - 1] + b * (R2.m[i - 1] - R1.m[i - 1]) + a * (Rin[i - 1] - R1.m[i - 1])

    R2.m[i] <- R2.m[i - 1] + c * (R1.m[i - 1] - R2.m[i - 1])
    
  }#can use three parameters instead

  c ~ dunif(0, 1)
  b <- a * b.coef
  b.coef ~ dunif(0, 1)
  a ~ dunif(0, 1)

  R2.m[1] ~ dnorm(R1.m[1], Sr.pre.2) 
  R1.m[1] ~ dnorm(Rpri.mean, Rpri.pre)
  
  #generate time series of input values
  #Rin is the input ratio, which is modeled with a switch point in the time series
  for(i in 1:t){
    
    Rin[i] ~ dnorm(ifelse(i > switch, Raft.mean, Rpri.mean), ifelse(i > switch, Raft.pre, Rpri.pre))
    #When i is greater than the switch point, Rin ~ dnorm(Raft.mean, Rin.m.pre)
  }
  
  switch ~ dcat(pi)
  
  pi <- c(pi.int, pi.switch, pi.end)
    
  pi.end <- rep(0, t - date - err.date)
  
  pi.switch <- rep(1, 1 + 2 * err.date)
  
  pi.int <- rep(0, date - err.date - 1)
  
  #uncertainty of the date of switch = +- the number of days
  err.date <- 2 
  
  #suspected date of the switch
  date <- 85


  Rpri.mean ~ dnorm(Rpri, Rpri.pre)
  
  Raft.mean ~ dnorm(Raft, Raft.pre)
  
  Raft.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.Raft)
  Rpri.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.Rpri)
  Sr.pre.rate.Raft <- 2e-5
  Sr.pre.rate.Rpri <- 2e-6

  Sr.pre.2 ~ dgamma(Sr.pre.shape, Sr.pre.rate.2)  
  Sr.pre.rate.2 <- 5e-7
  
  Sr.pre.1 ~ dgamma(Sr.pre.shape, Sr.pre.rate.1) 
  
  Sr.pre.shape <- 100
  Sr.pre.rate.1 <- 5e-6
  
}