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

#parameters a, b, and c from parameter estimates in file 3
params.forw.m <- mvrnorm(1, mu = turnover.params.mu, Sigma = turnover.params.vcov)

a <- exp(params.forw.m[1])# a can be substituted

b <- exp(params.forw.m[2])

c <- exp(params.forw.m[3])

#number of days in the simulation
t <- 700
input.120 <- initiate.switch(t, n.switch=5, day.switch=100, a=0.707, gap=0.003, duration=120)
plot(1:t,input,type="l")


###function 3: model micromill results given growth rate and sampling interval####
###1 declaration of fixed parameters, no stocastic relationships!!!
#try to use micromill parameters

#micromill parameters
max.dist.cal <- 12000

Ivo.rate.cal.mean <- 19.9
Ivo.rate.cal.sd <- 5.4

cal.intv <- 1000

dist.cal <- seq(12000, 0, by = - cal.intv)

n.cal <- length(dist.cal)



turnover.params.mu
turnover.params.vcov

#model components:

input.120

Se.bone.res <- Se.bone.forw.m(t = 700, input = input.120, a = a, b = b, c = c, Rs.int = NULL, Rb.int = NULL)
Se.bone.res[[1]]
Se.bone.res[[2]]


dist.cal.m <- rep(0,t) #initiate dist

dist.cal.m[1] <- max.dist.cal#maximum distance from the pulp cavity in microns

Ivo.rate.cal <- rnorm(t, Ivo.rate.cal.mean, Ivo.rate.cal.sd) #ivory growth rate, micron/day
# Ivo.rate.cal <- MCMC.ts.Ivo.rate.89.woi[[1]]

#model ivory growth
for (i in 2:t){
  dist.cal.m[i] <- dist.cal.m[i - 1] - Ivo.rate.cal[i] #simulate daily distance increment
}

#model averaging
Rs.cal.eva <- rep(0,n.cal)#initiate vector
for (i in 1:n.cal){
  vect.upp <- as.integer((dist.cal[i] + cal.intv/2) < dist.cal.m)
  vect.low <- as.integer((dist.cal[i] - cal.intv/2) < dist.cal.m)
  Rs.cal.eva[i] <- (Se.bone.res[[1]] %*% (vect.upp - vect.low))/sum(vect.upp - vect.low)
}



plot(0,0, xlim = c(1,t), ylim = c(0.706, 0.711), xlab = "days", ylab ="Sr 87/86",main="120 day")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines(1:t, input.120)
lines(1:t,Se.bone.res[[1]],lwd=2, col = "blue")
points((max(dist.cal) - dist.cal)/Ivo.rate.mic + 1, Rs.cal.eva, lwd = 2, col = "red")

plot(0,0, xlim = c(1,t), ylim = c(0.706, 0.711), xlab = "days", ylab ="Sr 87/86",main="100 day")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines(1:t, input.50)
lines(1:t,Se.bone.res[[1]],lwd=2, col = "blue")
points((max(dist.cal) - dist.cal)/Ivo.rate.mic + 1, Rs.cal.eva, lwd = 2, col = "red")

#plot distance
plot(0,0, xlim = c(8000,-1000), ylim = c(0.705, 0.716), xlab = "dist", ylab ="Sr 87/86",
     main="a = 0.02, rate.mean = 19.9")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines(dist.cal.m, Rs.cal, lwd = 2, col = "blue") 
lines(dist.mic, Rs.cal.eva, lwd = 2, col = "red") #approximate results from micromill
points(dist.mic, R.mic)

Rs.cal


#######forward model with prescribed input values of dietary switch####


######Forward model with prescribed bone value#####