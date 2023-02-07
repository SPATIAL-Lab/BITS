model {

  flux.ratio <- a/b
  pool.ratio <- c/b
  #Data eveluation

  for (i in 1:n.mea){

    R1.eva[i] <- inprod(R1.m, ((dist.mea[i] + s.intv/2) < dist) - ((dist.mea[i] - s.intv/2)  < dist))/sum(((dist.mea[i] + s.intv/2) < dist) - ((dist.mea[i] - s.intv/2)  < dist))

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

  c <- b * c.coef
  c.coef ~ dunif(0.01, 2.5) 
  b <- a * b.coef
  b.coef ~ dunif(0, 1)
  a ~ dunif(0, 1)
  
  R2.m[1] ~ dnorm(R1.m[1], Sr.pre.2) 
  R1.m[1] ~ dnorm(Rpri.mean, Rpri.pre)
  
  for(i in 1:t){
    
    Rin[i] ~ dnorm(ifelse(i > switch, Raft.mean, Rpri.mean), ifelse(i > switch, Raft.pre, Rpri.pre))

  }
  
  switch ~ dcat(pi)

  pi <- c(pi.int, pi.switch, pi.end)
    
  pi.end <- rep(0, t - date - err.date)
  
  pi.switch <- rep(1, 1 + 2 * err.date)
  
  pi.int <- rep(0, date - err.date - 1)
  
  err.date <- 2 
  
  date <- 85

  #allowing some uncertainty in R0 values
  Rpri.mean ~ dnorm(Rpri, Rpri.pre)
  
  Raft.mean ~ dnorm(Raft, Raft.pre)
  
  #Re has more uncertainty than Re
  Raft.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.Raft)
  Rpri.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.Rpri)
  Sr.pre.rate.Raft <- 2e-5
  Sr.pre.rate.Rpri <- 2e-6
  
  
  #precision for average Sr measurements in bone should be much smaller
  Sr.pre.2 ~ dgamma(Sr.pre.shape, Sr.pre.rate.2)  
  Sr.pre.rate.2 <- 5e-7
  
  Sr.pre.1 ~ dgamma(Sr.pre.shape, Sr.pre.rate.1) 
  
  Sr.pre.shape <- 100
  Sr.pre.rate.1 <- 5e-6
}