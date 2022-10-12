##forward model with t iterations

####helper functions###
###Function 1 create input time series in number of days#####
#function to generate a number of switches between Sr 87/86 values a and b

initiate.switch <- function(t, n.switch, day.switch, a, gap, duration){

    #initial vector
    v.init <- rep(a, (day.switch -1))
    v.switch <- rep(a+gap, duration)
    v.back <- rep(a, duration)
    sw.d <- n.switch*duration
    
    v.switches <- c(v.switch, v.back)

    rep.switches <- rep(v.switches,ceiling(n.switch/2))
    
    if((day.switch -1 + sw.d) < t){#switch is completed before the total number of days
      t.end <- t-(day.switch -1 + sw.d)
      
      if((n.switch %% 2) == 0){#even number of switches
        r.end <- a
      }
      else{#odd number of switches
        r.end <- a + gap
      }
      v.end <- rep(r.end, t.end)
      all.input <- c(v.init, rep.switches[1:sw.d], v.end)
    }
    else{#switch is not completed before the total number of days
      #no need to worry about end values
      all.input <- c(v.init, rep.switches[1:sw.d])#take advantage of the vector recycling feature
    }
    Sr.input <- all.input[1:t]
    return(Sr.input)
}

###function 2: model serum and bone values####
Se.bone.forw.m <- function(t, input, a, b, c, Rs.int, Rb.int){
  if(length(input != t)){
    #warning("Length of input has to be same as t. Vector input is being recycled")
    Rin <- rep(0,t) #Initiate Rin
    Rin <- Rin + input #recycle input vector
  }
  else{
    Rin <- input
  }
  Rb <- rep(0,t) #initiate rb
  Rs <- rep(0,t) #initiate rb
  if(is.null(Rs.int)){
    Rs[1] <- input[1]
  }
  if(is.null(Rb.int)){
    Rb[1] <- input[1]
  }
  else{
    Rb[1] <- Rb.int
    Rs[1] <- Rs.int
  }
  for (i in 2:t){ #generate serum and bone series
    #serum ratio
    Rs[i] <- Rs[i - 1] + b * (Rb[i - 1] - Rs[i - 1]) + a * (Rin[i - 1] - Rs[i - 1])
    
    #bone ratios
    Rb[i] <- Rb[i - 1] + c * (Rs[i - 1] - Rb[i - 1])
  }
  return(list(Rs, Rb))
}


###function 3: ivory growth simulation using mean and sd of growth rate####

gorwth.sim <- function(t, max.dist, rate, rate.sd){
  dist.m <- rep(0,t) #initiate dist
  
  dist.m[1] <- max.dist#maximum distance from the pulp cavity in microns
  
  #generate random rates
  Ivo.rate <- rnorm(t, rate, rate.sd) #ivory growth rate, micron/day
  
  #model ivory growth
  for (i in 2:t){
    dist.m[i] <- dist.m[i - 1] - Ivo.rate[i] #simulate daily distance increment
  }
  return(dist.m)
}

###function 4: model micromill results given growth rate and sampling interval####
micromill.sim <- function(max.dist, intv, dist.m, Rs.m){
  dist.mic <- seq(max.dist, 0, by = - intv)
  n.mic <- length(dist.mic)
  #model averaging
  Rs.eva <- rep(0, n.mic)#initiate vector
  for (i in 1:n.mic){
    vect.upp <- as.integer((dist.mic[i] + intv/2) < dist.m)
    vect.low <- as.integer((dist.mic[i] - intv/2) < dist.m)
    Rs.eva[i] <- (Rs.m %*% (vect.upp - vect.low))/sum(vect.upp - vect.low)
  }
  return(list(dist.mic, Rs.eva))
}

####begin forward model####

#1 check turnover parameters#
#parameters a, b, and c from parameter estimates in file 3

#2 extract turnover parameters, use MAP in the forward model
a <- MCMC.CI.a[[1]]
b <- MCMC.CI.b[[1]]
c <- MCMC.CI.c[[1]]

#3 number of days in the simulation
t <- 700

#4 generate input series with fixed durations#
input.120 <- initiate.switch(t, n.switch=5, day.switch=100, a=0.707, gap=0.003, duration=120)

#5 generate serum and bone series based on input series and turnover parameters#
Se.bone.res <- Se.bone.forw.m(t = 700, input = input.120, a = a, b = b, c = c, Rs.int = NULL, Rb.int = NULL)

#6 ivory growth parameters, unit in microns
max.dist <- 12000 
Ivo.rate.mean <- 19.9 #microns per day
Ivo.rate.sd <- 5.4

#7 ivory growth simulation, returns a distance vector of length t
ivo.dist.m <- gorwth.sim(t = 700, max.dist = max.dist, rate = Ivo.rate.mean, rate.sd = Ivo.rate.sd)

#8 micromill parameter: milling width in microns
intv <- 1000

#9 micromill averaging simulation, returns a list of two vectors
#first vector is micromill distance, the center of the milling width
#second vector is
res.mic <- micromill.sim(max.dist = max.dist, intv = intv, dist.m = ivo.dist.m, Rs.m = Se.bone.res[[1]])
######end of forward model######

###plots###
#plot input, serum history, and simulated micromill values
plot(0,0, xlim = c(1,t), ylim = c(0.706, 0.711), xlab = "days", ylab ="Sr 87/86",main="120 day")
lines(1:t, input.120)
lines(1:t, Se.bone.res[[1]],lwd=2, col = "blue")
#converting micromill distance to ~days using rate Ivo.rate.mic
points((max(res.mic[[1]]) - res.mic[[1]])/Ivo.rate.mean + 1, res.mic[[2]], lwd = 2, col = "red")

#plot simulated ivory measurement with distance in the x axis
plot(0,0, xlim = c(max(res.mic[[1]]),-1000), ylim = c(0.706, 0.711), xlab = "dist", ylab ="Sr 87/86",
     main="")
lines(ivo.dist.m, Se.bone.res[[1]], lwd = 2, col = "blue") 
points(res.mic[[1]], res.mic[[2]], col = "red")

#########Forward model to demonstrate carry over effect, synthetic input 120 days interval#######
####begin forward model####

#1 check turnover parameters#
#parameters a, b, and c from parameter estimates in file 3

#2 extract turnover parameters, use MAP in the forward model
a <- MCMC.CI.a[[1]]
b <- MCMC.CI.b[[1]]
c <- MCMC.CI.c[[1]]

#3 number of days in the simulation


#4 generate input series with fixed durations#
#Mid -> low -> mid -> High -> Mid, 90 day interval
syn.mid <- 0.709
syn.high <- 0.713
syn.low <- 0.705

syn.input.120 <- c(rep(syn.mid,120), rep(syn.low,120), rep(syn.mid,120), rep(syn.high,120), rep(syn.mid,120))

#5 generate serum and bone series based on input series and turnover parameters#
Se.bone.res <- Se.bone.forw.m(t = length(syn.input.120), input = syn.input.120, a = a, b = b, c = c, Rs.int = NULL, Rb.int = NULL)

#6 ivory growth parameters, unit in microns
max.dist <- 12000 
Ivo.rate.mean <- 19.9 #microns per day
Ivo.rate.sd <- 5.4

#7 ivory growth simulation, returns a distance vector of length t
ivo.dist.m <- gorwth.sim(t = length(syn.input.120), max.dist = max.dist, rate = Ivo.rate.mean, rate.sd = Ivo.rate.sd)

#8 micromill parameter: milling width in microns
intv <- 1000

#9 micromill averaging simulation, returns a list of two vectors
#first vector is micromill distance, the center of the milling width
#second vector is
res.mic <- micromill.sim(max.dist = max.dist, intv = intv, dist.m = ivo.dist.m, Rs.m = Se.bone.res[[1]])
######end of forward model######

###plots###
#plot input, serum history, and simulated micromill values
plot(0,0, xlim = c(1,length(syn.input.120)), ylim = c(syn.low, syn.high), xlab = "days", ylab ="Sr 87/86",main="120 day")
lines(1:length(syn.input.120), syn.input.120)
lines(1:length(syn.input.120), Se.bone.res[[1]],lwd=2, col = "blue")
#converting micromill distance to ~days using rate Ivo.rate.mic
points((max(res.mic[[1]]) - res.mic[[1]])/Ivo.rate.mean + 1, res.mic[[2]], lwd = 2, col = "red")

#plot simulated ivory measurement with distance in the x axis
plot(0,0, xlim = c(max(res.mic[[1]]),-1000), ylim = c(syn.low, syn.high), xlab = "dist", ylab ="Sr 87/86",
     main="")
lines(ivo.dist.m, Se.bone.res[[1]], lwd = 2, col = "blue") 
points(res.mic[[1]], res.mic[[2]], col = "red")


###############synthetic input 30 day interval to demonstrate carry over effect########
####begin forward model####

#1 check turnover parameters#
#parameters a, b, and c from parameter estimates in file 3

#2 extract turnover parameters, use MAP in the forward model
a <- MCMC.CI.a[[1]]
b <- MCMC.CI.b[[1]]
c <- MCMC.CI.c[[1]]

#3 number of days in the simulation


#4 generate input series with fixed durations#
#Mid -> low -> mid -> High -> Mid, 90 day interval
syn.mid <- 0.709
syn.high <- 0.713
syn.low <- 0.705

syn.input.30 <- c(rep(syn.mid,30), rep(syn.low,30), rep(syn.mid,30), rep(syn.high,30), rep(syn.mid,30))
syn.input.30 <- rep(syn.input.30,3)

#5 generate serum and bone series based on input series and turnover parameters#
Se.bone.res <- Se.bone.forw.m(t = length(syn.input.30), input = syn.input.30, a = a, b = b, c = c, Rs.int = NULL, Rb.int = NULL)

#6 ivory growth parameters, unit in microns
max.dist <- 12000 
Ivo.rate.mean <- 19.9 #microns per day
Ivo.rate.sd <- 5.4

#7 ivory growth simulation, returns a distance vector of length t
ivo.dist.m <- gorwth.sim(t = length(syn.input.30), max.dist = max.dist, rate = Ivo.rate.mean, rate.sd = Ivo.rate.sd)

#8 micromill parameter: milling width in microns
intv <- 1000

#9 micromill averaging simulation, returns a list of two vectors
#first vector is micromill distance, the center of the milling width
#second vector is
res.mic <- micromill.sim(max.dist = max.dist, intv = intv, dist.m = ivo.dist.m, Rs.m = Se.bone.res[[1]])
######end of forward model######

###plots###
#plot input, serum history, and simulated micromill values
plot(0,0, xlim = c(1,length(syn.input.30)), ylim = c(syn.low, syn.high), 
     xlab = "days", ylab ="Sr 87/86",main="30 day")
lines(1:length(syn.input.30), syn.input.30)
lines(1:length(syn.input.30), Se.bone.res[[1]],lwd=2, col = "blue")
#converting micromill distance to ~days using rate Ivo.rate.mic
points((max(res.mic[[1]]) - res.mic[[1]])/Ivo.rate.mean + 1, res.mic[[2]], lwd = 2, col = "red")

#plot simulated ivory measurement with distance in the x axis
plot(0,0, xlim = c(max(res.mic[[1]]),-1000), ylim = c(syn.low, syn.high), xlab = "dist", ylab ="Sr 87/86",
     main="")
lines(ivo.dist.m, Se.bone.res[[1]], lwd = 2, col = "blue") 
points(res.mic[[1]], res.mic[[2]], col = "red")


######Forward model with prescribed bone value#####