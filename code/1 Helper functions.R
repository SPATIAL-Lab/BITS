####helper functions#####
MCMC.ts <- function (MCMC.res){
  require(KernSmooth)
  require(bayestestR)
  dim.MCMC <- dim(MCMC.res)
  #typically, the first element is # of interation
  #the second element is the time series
  map.res <- rep(0, dim.MCMC[2])
  hdi.high <- rep(0, dim.MCMC[2])
  hdi.low <- rep(0, dim.MCMC[2])
  
  for(i in 1:dim.MCMC[2]){
    map.res[i] <- map_estimate(MCMC.res[,i], method = "KernSmooth")
    hdi.89 <- hdi(MCMC.res[,i], ci = 0.89)
    hdi.low[i] <- hdi.89$CI_low
    hdi.high[i] <- hdi.89$CI_high
  }
  
  return (list(map.res, hdi.low, hdi.high))
}

MCMC.ts.dist <- function(MCMC.res, MCMC.dist, MCMC.index){
  require(scales)
  
  #typically, the first element is # of interation
  #the second element is the time series
  
  dim.MCMC <- dim(MCMC.res) #this is 800
  n <- dim.MCMC[1] #number of interation
  dim.MCMC.dist <- dim(MCMC.dist) #this is 100
  MCMC.res.index <- array(0, c(n, dim.MCMC.dist[2]))
  for(i in 1:n){ #for each iteration, extract MCMC index results
    MCMC.res.index[i,] <- MCMC.res[i, MCMC.index[i,]]
  }
  
  for(i in 1:n){
    lines(MCMC.dist[i,], MCMC.res.index[i,], col = alpha("black", 0.02))
  }
  
  return(MCMC.res.index)
}

MCMC.dist.median <- function(MCMC.res){
  dim.MCMC <- dim(MCMC.res)
  n <- dim.MCMC[2] #number of data points
  MCMC.res.med <- rep(0, n)
  MCMC.res.max <- rep(0, n)
  MCMC.res.min <- rep(0, n)
  for(i in 1:n){ #for each iteration, extract MCMC index results
    MCMC.res.med[i] <- median (MCMC.res[,i])
    MCMC.res.min[i] <- min (MCMC.res[,i])
    MCMC.res.max[i] <- median (MCMC.res[,i])
  }
  return(MCMC.res.med)
}