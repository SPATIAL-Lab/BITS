######forward model########

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

###function 2: model P1 and P2 values####
forw.m <- function(t, input, a, b, c, R1.int, R2.int){
  if(length(input != t)){
    #warning("Length of input has to be same as t. Vector input is being recycled")
    Rin <- rep(0,t) #Initiate Rin
    Rin <- Rin + input #recycle input vector
  }
  else{
    Rin <- input
  }
  R2 <- rep(0,t) #initiate 
  R1 <- rep(0,t) #initiate 
  if(is.null(R1.int)){
    R1[1] <- input[1]
  }
  if(is.null(R2.int)){
    R2[1] <- input[1]
  }
  else{
    R2[1] <- R2.int
    R1[1] <- R1.int
  }
  for (i in 2:t){ #generate serum and bone series
    R1[i] <- R1[i - 1] + b * (R2[i - 1] - R1[i - 1]) + a * (Rin[i - 1] - R1[i - 1])
    
    R2[i] <- R2[i - 1] + c * (R1[i - 1] - R2[i - 1])
  }
  return(list(R1, R2))
}

####begin forward model####
#loading MAP estimates from posterior of the calibration
a <- MAP.a[1]
b <- MAP.b[1]
c <- MAP.c[1]
c/b

a/b
#to change pool ratio, change c; to change flux ratio, change a
# flux.ratio <- a/b
# pool.ratio <- c/b

#2 number of days in the simulation
t <- 900

#3 generate input series with fixed duration#
input.misha <- initiate.switch(t, n.switch=1, day.switch=100, a=0.706, gap=0.005, duration=360)

#4 generate serum and bone series based on input series and turnover parameters#
res <- forw.m(t = 900, input = input.misha, a = a, b = b, c = c, R1.int = NULL, R2.int = NULL)
res.misha <- res[[1]]

###############change pool ratio: make c larger##############
a <- MAP.a[1]
b <- MAP.b[1]
c <- 0.01 #MAP.c = 0.0041
c/b

#3 generate input series with fixed duration#
input.misha <- initiate.switch(t, n.switch=1, day.switch=100, a=0.706, gap=0.005, duration=360)

#4 generate serum and bone series based on input series and turnover parameters#
res <- forw.m(t = 900, input = input.misha, a = a, b = b, c = c, R1.int = NULL, R2.int = NULL)
res.prl <- res[[1]]

###############change pool ratio: make c smaller##############
a <- MAP.a[1]
b <- MAP.b[1]
c <- 0.001 #MAP.c = 0.0041
c/b

#2 number of days in the simulation
t <- 900

#3 generate input series with fixed duration#
input.misha <- initiate.switch(t, n.switch=1, day.switch=100, a=0.706, gap=0.005, duration=360)

#4 generate serum and bone series based on input series and turnover parameters#
res <- forw.m(t = 900, input = input.misha, a = a, b = b, c = c, R1.int = NULL, R2.int = NULL)
res.prs <- res[[1]]


###############change flux ratio: make a larger##############
a <- 0.04 #MAP.a = 0.0169
b <- MAP.b[1]
c <- MAP.c[1]
a/b

#2 number of days in the simulation
t <- 900

#3 generate input series with fixed duration#
input.misha <- initiate.switch(t, n.switch=1, day.switch=100, a=0.706, gap=0.005, duration=360)

#4 generate serum and bone series based on input series and turnover parameters#
res <- forw.m(t = 900, input = input.misha, a = a, b = b, c = c, R1.int = NULL, R2.int = NULL)
res.frl <- res[[1]]

###############change flux ratio: make a smaller ##############
a <- MAP.b[1] #MAP.a = 0.0169
b <- MAP.b[1] #MAP.b = 0.0141
c <- MAP.c[1]
a/b

#2 number of days in the simulation
t <- 900

#3 generate input series with fixed duration#
input.misha <- initiate.switch(t, n.switch=1, day.switch=100, a=0.706, gap=0.005, duration=360)

#4 generate serum and bone series based on input series and turnover parameters#
res <- forw.m(t = 900, input = input.misha, a = a, b = b, c = c, R1.int = NULL, R2.int = NULL)

res.frs <- res[[1]]

#########Forward model to demonstrate signal attenuation effect, synthetic inputs#######
####begin forward model####
#1 extract turnover parameters, use MAP in the forward model
a <- MAP.a[1]
b <- MAP.b[1]
c <- MAP.c[1]

#2 number of days in the simulation
#t <- 700 #omitted here due to synthetic input data 

#3 generate input series with fixed duration#
#Mid -> low -> mid -> High -> Mid
syn.mid <- 0.709
syn.high <- 0.711
syn.low <- 0.707

syn.input.150 <- c(rep(syn.mid,30), rep(syn.high,150),rep(syn.mid,30),rep(syn.low,150)  )
syn.input.150 <- rep(syn.input.150,2)
#4 generate serum and bone series based on input series and turnover parameters#
res150 <- forw.m(t = length(syn.input.150), input = syn.input.150, a = a, b = b, c = c, R1.int = NULL, R2.int = NULL)
######end of forward model######

#########Forward model to demonstrate carry over effect, synthetic input 120 days interval#######
####begin forward model####
#1 extract turnover parameters, use MAP in the forward model
a <- MAP.a[1]
b <- MAP.b[1]
c <- MAP.c[1]

#2 number of days in the simulation
#t <- 700 #omitted here due to synthetic input data 

#3 generate input series with fixed duration#
#Mid -> low -> mid -> High -> Mid
syn.mid <- 0.709
syn.high <- 0.711
syn.low <- 0.707

syn.input.120 <- c(rep(syn.mid,60), rep(syn.high,120),rep(syn.mid,60),rep(syn.low,120)  )
syn.input.120 <- rep(syn.input.120,2)
#4 generate serum and bone series based on input series and turnover parameters#
res120 <- forw.m(t = length(syn.input.120), input = syn.input.120, a = a, b = b, c = c, R1.int = NULL, R2.int = NULL)

###############synthetic input 60 day interval to demonstrate carry over effect########

#1 extract turnover parameters, use MAP in the forward model
#parameters a, b, and c from parameter estimates in file 3
a <- MAP.a[1]
b <- MAP.b[1]
c <- MAP.c[1]

#2 number of days in the simulation
#t <- 700 #omitted here due to synthetic input data 

#3 generate input series with fixed duration#
#Mid -> low -> mid -> High -> Mid
syn.mid <- 0.709
syn.high <- 0.711
syn.low <- 0.707

syn.input.60 <- c(rep(syn.mid,30), rep(syn.high,60),rep(syn.mid,30),rep(syn.low,60) )
syn.input.60 <- rep(syn.input.60,4)

#4 generate serum and bone series based on input series and turnover parameters#
res60 <- forw.m(t = length(syn.input.60), input = syn.input.60, a = a, b = b, c = c, R1.int = NULL, R2.int = NULL)
######end of forward model######

#######compared to the reaction progress approach######
####begin forward model####
#loading MAP estimates from posterior of the calibration
a <- MAP.a[1]
b <- MAP.b[1]
c <- MAP.c[1]
c/b

a/b
#to change pool ratio, change c; to change flux ratio, change a
# flux.ratio <- a/b
# pool.ratio <- c/b

#2 number of days in the simulation
t <- 900

#3 generate input series with fixed duration#
input.misha <- initiate.switch(t, n.switch=1, day.switch=100, a=0.706, gap=0.005, duration=360)

#4 generate serum and bone series based on input series and turnover parameters#
res <- forw.m(t = 900, input = input.misha, a = a, b = b, c = c, R1.int = NULL, R2.int = NULL)
res.misha <- res[[1]]

#convert to reaction progress

input.misha #from 0.706 to 0.711
res.misha #dxt

reac.prog <- (0.711-res.misha)/(0.711-0.706)

reac.prog.tk <- reac.prog[101:length(reac.prog)]

#after the inflection point: [200:800]

lm.res.misha.2 <- lm(log(reac.prog.tk[200:800]) ~ c(200:800))
summary(lm.res.misha.2)
#intercept = -0.6539
#slope = -0.002102
exp(-0.6539) #f2 = 0.52
-1/-0.002102*log(2) #t1/2 = 330 days

#estimate of parameter c:
#c ~= -0.002102 /(1-exp(-0.654))
0.002102/(1-exp(-0.6539))

#fit for pool 1

reac.prog.res1 <- log(reac.prog.tk[1:150]-exp(1:150 * -0.002102 -0.6539)) #first residual
plot(c(1:150),reac.prog.res1,ylab="residual: ln(1-F)",xlab="Days")
lm.res.misha.1 <- lm(reac.prog.res1 ~ c(1:150))
summary(lm.res.misha.1)
#intercept = -0.7252
#slope = -0.03378
exp(-0.7252) #f1 = 0.48
-1/-0.03378*log(2) #t1/2 = 20.5 days

#estimate of parameter a:
#a ~= 0.03378*exp(-0.7252)
0.03378*exp(-0.7252)

0.52/0.48 # f2/f1 ~= flux ratio

#################use a different ratio###################
####begin forward model####
#loading MAP estimates from posterior of the calibration
a <- 0.04
b <- MAP.b[1]
c <- MAP.c[1]
a/b
#to change pool ratio, change c; to change flux ratio, change a
# flux.ratio <- a/b
# pool.ratio <- c/b

#2 number of days in the simulation
t <- 900

#3 generate input series with fixed duration#
input.misha <- initiate.switch(t, n.switch=1, day.switch=100, a=0.706, gap=0.005, duration=360)

#4 generate serum and bone series based on input series and turnover parameters#
res <- forw.m(t = 900, input = input.misha, a = a, b = b, c = c, R1.int = NULL, R2.int = NULL)
res.frl <- res[[1]]

#convert to reaction progress

input.misha #from 0.706 to 0.711
res.frl #dxt

reac.prog.frl <- (0.711-res.frl)/(0.711-0.706)

reac.prog.frl.tk <- reac.prog.frl[101:length(reac.prog.frl)]

#after the inflection point: [200:800]

lm.res.frl.2 <- lm(log(reac.prog.frl.tk[200:800]) ~ c(200:800))
plot(1:length(reac.prog.frl.tk), log(reac.prog.frl.tk),ylab="ln(1-F)",xlab="Days",main="Slow turnover pool")
summary(lm.res.frl.2)
#intercept = -1.235
#slope = -0.00297
exp(-1.235) #f2 = 0.291
-1/-0.00297*log(2) #t1/2 = 233 days

#estimate of parameter c:
#for pool 2, c ~= -1 *  (1-exp(-1.235)) * -0.00297
0.00297/(1-exp(-1.235))

#fit for pool 1
reac.prog.frl.res1 <- log(reac.prog.frl.tk[1:150]-exp(1:150 * -0.00297 -1.235)) #first residual
plot(c(1:150),reac.prog.frl.res1,ylab="residual: ln(1-F)",xlab="Days")
lm.res.frl.1 <- lm(reac.prog.frl.res1 ~ c(1:150))
summary(lm.res.frl.1)
#intercept = -0.3111
#slope = -0.05775
exp(-0.3111) #f1 = 0.73
-1/-0.05775*log(2) #t1/2 = 12 days

#estimate of parameter a:
0.05775*exp(-0.3111) 

0.73/0.291 # f2/f1 ~= flux ratio
