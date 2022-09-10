####helper functions#####
MCMC.CI.bound <- function (MCMC.res, CI){
  require(KernSmooth)
  require(bayestestR)
  dim.MCMC <- dim(MCMC.res)
  #typically, the first element is # of iterations
  #the second element is the time series
  map.res <- rep(0, dim.MCMC[2])
  hdi.high <- rep(0, dim.MCMC[2])
  hdi.low <- rep(0, dim.MCMC[2])
  
  for(i in 1:dim.MCMC[2]){
    map.res[i] <- map_estimate(MCMC.res[,i], method = "KernSmooth")
    hdi <- hdi(MCMC.res[,i], ci = CI)
    hdi.low[i] <- hdi$CI_low
    hdi.high[i] <- hdi$CI_high
  }
  
  return (list(map.res, hdi.low, hdi.high, CI))
}

#modeled sampling distance to compare with raw data, can be used to visualize serum values (Rs.m)
MCMC.dist.plot <- function(MCMC.res, MCMC.dist){#t is the number of days in the model
  require(scales)
  
  #typically, the first element is # of iterations
  #the second element is the time series
  #t is the number of days in the model
  
  dim.MCMC <- dim(MCMC.res)
  n <- dim.MCMC[1] #number of iterations
  
  for(i in 1:n){
    lines(MCMC.dist[i,], MCMC.res[i,], col = alpha("black", 0.02))
  }
  
}

#reconstructed time line, could be serum(Rs.m), input(Rin), modeled distance(dist).
MCMC.tl.plot <- function(MCMC.res, t){#t is the number of days in the model
  require(scales)
  
  #typically, the first element is # of iterations
  #the second element is the time series
  #t is the number of days in the model
  
  dim.MCMC <- dim(MCMC.res)
  n <- dim.MCMC[1] #number of iterations

  for(i in 1:n){
    lines(1:t, MCMC.res[i,], col = alpha("black", 0.02))
  }

}

MCMC.dist.median <- function(MCMC.res){
  dim.MCMC <- dim(MCMC.res)
  n <- dim.MCMC[2] #number of data points
  MCMC.res.med <- rep(0, n)
  MCMC.res.max <- rep(0, n)
  MCMC.res.min <- rep(0, n)
  for(i in 1:n){ #for each iteration, extract MCMC index results
    MCMC.res.med[i] <- median (MCMC.res[,i])
  }
  return(MCMC.res.med)
}

pri.multi.norm.den <- function(X.min,X.max,mu,vcov){
  require(OpenMx)
  X.range <- X.max - X.min
  X.interval <- X.range/512
  X.multinorm <- seq(from = X.min, to = X.max, by = X.interval)
  dim.vcov <- dim(vcov)
  
  #initiate density matrix
  density <- as.data.frame(matrix(0, nrow=length(X.multinorm), ncol = dim.vcov[1]))
  
  #initiate index matrix
  index <- as.data.frame(matrix(nrow=dim.vcov[1], ncol = 2))
  for (i in 1: dim.vcov[1]){
    index[1,i] <- X.min
    index[2,i] <- X.max
  }
  
  for (i in 1: dim.vcov[1]){
    for(j in 2:length(X.multinorm)){
      #assigning new values
      index[1,i] <- X.multinorm[j-1]
      index[2,i] <- X.multinorm[j]
      #integrating one interval for density[i]
      density[i,j-1] <- omxMnor(vcov,mu,index[1,],index[2,])
    }
    #revert to the default values
    index[1,i] <- X.min
    index[2,i] <- X.max
  }
  
  #scaling density values to 1
  density.sc <- density[,]/X.interval
  
  #monitoring the density scaling approximation
  #the sum should be close to 1
  
  for (i in 1: dim.vcov[1]){
    print(sum(density[i,]))
  }
  
  results<-list(x = (X.multinorm[1:512]+ X.interval/2), y = density.sc)
  return(results)
}
