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

#modeled sampling distance to compare with raw data, can be used to visualize Pool1 values (R1.m)
MCMC.dist.plot <- function(MCMC.res, MCMC.dist){#t is the number of days in the model
  require(scales)
  
  #typically, the first element is # of iterations
  #the second element is the time series
  #t is the number of days in the model
  
  dim.MCMC <- dim(MCMC.res)
  n <- dim.MCMC[1] #number of iterations
  
  for(i in 1:n){
    lines(MCMC.dist[i,], MCMC.res[i,], col = alpha("black", 0.06))
  }
  
}

#reconstructed time line, could be serum(R1.m), input(Rin), modeled distance(dist).
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

#Function 5: plotting points with error bars
PlotPE <- function(x, y, error.y, col){
  if(length(x) != length(y)){
    warning("Error in plot.xy, x and y have different lengths")
  }
  if(length(x) != length(error.y)){
    warning("Error in plot.xy, x and error.y have different lengths")
  }
  if(length(y) != length(error.y)){
    warning("Error in plot.xy, y and error.y have different lengths")
  }
  if(is.null(col)){
    warning("Please specify plotting color")
  }
  
  n <- length(x)
  
  for(i in 1:n){
    arrows(x[i], y[i], x[i], y[i] + error.y[i], length = 0.03, angle = 90, col = col)
    arrows(x[i], y[i], x[i], y[i] - error.y[i], length = 0.03, angle = 90, col = col)
  }
  points (x, y, pch = 16, col = col)
  
}