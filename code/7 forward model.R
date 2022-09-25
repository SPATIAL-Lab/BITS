##forward model with t iterations

###1 declaration of fixed parameters
#try to use micromill parameters
#t <- 450
t <- 500

Sr.pre.rate.Re <- 2e-5
Sr.pre.rate.R0 <- 2e-6

Sr.pre.shape <- 100
Sr.pre.s.rate <- 5e-6

Sr.pre.b.rate = 5e-7

Re.pre <- rgamma(1, shape=Sr.pre.shape, rate=Sr.pre.rate.Re)
R0.pre <- rgamma(1, shape=Sr.pre.shape, rate=Sr.pre.rate.R0)

Sr.pre.s <- rgamma(1, shape=Sr.pre.shape, rate=Sr.pre.s.rate)

Re.mean <- rnorm(1, Re, 1/sqrt(Re.pre))

R0.mean <- rnorm(1, R0, 1/sqrt(R0.pre))

Sr.pre.b <- rgamma(1, shape=Sr.pre.shape, rate=Sr.pre.b.rate)

# date <- 76
# 
# err.date <- 2
# 
# switch <- runif(1, date - err.date - 1, t - date - err.date)
# 
# switch <- round(switch, digits = 0)

# for(i in 1:(switch-1)){ #before switch
#   
#   Rin.cal[i] <- rnorm(1, R0.mean, 1/sqrt(R0.pre))
#   #When i is greater than the switch point, Rin ~ dnorm(Re.mean, Rin.m.pre)
# }
# 
# for(i in switch:t){ #after switch
#   
#   Rin.cal[i] <- rnorm(1, Re.mean, 1/sqrt(Re.pre))
#   #When i is greater than the switch point, Rin ~ dnorm(Re.mean, Rin.m.pre)
# }

#Rin.cal <- MCMC.inv.bmca.Rin.m.89[[1]] #use inversion generated median value as input
Rin.cal <- MCMC.ts.Rin.m.89.woi[[1]]

#micromill parameters
max.dist.cal <- 8000

Ivo.rate.cal.mean <- 19.9
Ivo.rate.cal.sd <- 5.4

cal.intv <- 400

dist.cal <- dist.mic

n.cal <- 20

Rb.cal <- rep(0,t) #initiate vector for input
Rs.cal <- rep(0,t) #initiate vector for input

Rs.cal[1] <- rnorm(1, Rin.cal[1], 1/sqrt(R0.pre))

Rb.cal[1] <- rnorm(1, Rs.cal[1], 1/sqrt(Sr.pre.b)) #use a different error term

# params.forw.m <- mvrnorm(1, mu = turnover.params.mu, Sigma = turnover.params.vcov)
# 
# a <- exp(params.forw.m[1])# a can be substituted
# 
# b <- exp(params.forw.m[2])
# 
# c <- exp(params.forw.m[3])
a <- 0.02

b <- 0.004

c <- 0.003


for (i in 2:t){ #generate serum and bone series
  #serum ratio
  Rs.cal[i] <- Rs.cal[i - 1] + b * (Rb.cal[i - 1] - Rs.cal[i - 1]) + a * (Rin.cal[i - 1] - Rs.cal[i - 1])

  #bone ratios
  Rb.cal[i] <- Rb.cal[i - 1] + c * (Rs.cal[i - 1] - Rb.cal[i - 1])
}

# Rs.cal <- MCMC.ts.Rs.m.89.woi[[1]]
# Rb.cal <- MCMC.ts.Rb.m.89.woi[[1]]



dist.cal.m <- rep(0,t) #initiate dist

dist.cal.m[1] <- max.dist.cal#maximum distance from the pulp cavity in microns

Ivo.rate.cal <- rnorm(t, Ivo.rate.cal.mean, Ivo.rate.cal.sd) #ivory growth rate, micron/day
# Ivo.rate.cal <- MCMC.ts.Ivo.rate.89.woi[[1]]

#model ivory growth
for (i in 2:t){
  dist.cal.m[i] <- dist.cal.m[i - 1] - Ivo.rate.cal[i] #simulate daily distance increment
}

Rs.cal.eva <- rep(0,n.cal)

#model averaging
for (i in 1:n.cal){
  vect.upp <- as.integer((dist.cal[i] + cal.intv/2) < dist.cal.m)
  vect.low <- as.integer((dist.cal[i] - cal.intv/2) < dist.cal.m)
  Rs.cal.eva[i] <- (Rs.cal %*% (vect.upp - vect.low))/sum(vect.upp - vect.low)

}

plot(0,0, xlim = c(1,500), ylim = c(0.705, 0.716), xlab = "days", ylab ="Sr 87/86")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines((max(dist.mic+200) - dist.mic)/Ivo.rate.mic + 1, R.mic, lwd = 2, col = "blue") #approximate results from micromill
lines((max(dist.mic+200) - dist.mic)/Ivo.rate.mic + 1, Rs.cal.eva, lwd = 2, col = "red")

#plot distance
plot(0,0, xlim = c(8000,-1000), ylim = c(0.705, 0.716), xlab = "dist", ylab ="Sr 87/86",
     main="a = 0.02, rate.mean = 19.9")
#converting micromill distance to ~days using rate Ivo.rate.mic
lines(dist.cal.m, Rs.cal, lwd = 2, col = "blue") 
lines(dist.mic, Rs.cal.eva, lwd = 2, col = "red") #approximate results from micromill
points(dist.mic, R.mic)

Rs.cal
