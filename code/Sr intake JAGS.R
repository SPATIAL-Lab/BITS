model {
  #evaluation
  R.intake.mean ~ dnorm(Rin, 1/R.intake.sd^2)
  
  #generate Sr values of food using a mixing process of food and water
  w.contrib = w.water/(w.food + w.water)
  
  h.contrib = w.hay/(w.food + w.water)
  #weighted mean calculated by Sr concentration and ratio
  Rin <- (sum.c.Sr.pel + sum.c.Sr.hay + sum.c.Sr.sup + c.Sr.w * w.intake)/(w.food + w.water)
  
  w.food = w.hay + w.sup + w.pel #total amount of Sr in food in mg
  
  #with the digestibility modifier for hay 
  #total amount of Sr in hay (mg), mass * concentration * digestibility * ratio
  w.hay = sum(c.hay.m) * diges.h * f.h / m.feed * f.intake 
  
  sum.c.Sr.hay = sum(c.Sr.hay) * diges.h * f.h / m.feed * f.intake
  
  #pellet is assumed to be homogenized and 60% digestible (Sr) (half the fiber content as hay)
  #calculate weight of the pellet
  w.pel = c.pel.m * f.pel * diges.p * p.intake #total amount of Sr in pellet
  
  sum.c.Sr.pel = Sr.pel.m * c.pel.m * f.pel * diges.p * p.intake
  
  #calculate weight of the supplement in kg
  w.sup = c.sup.m * (1- f.pel) * diges.s * p.intake #total amount of Sr in pellet
  
  sum.c.Sr.sup = Sr.sup.m * c.sup.m * (1- f.pel) * diges.s * p.intake
  
  #Water Sr values
  #calculate weight of the water
  w.water = w.intake * c.w.m #total amount of Sr in water intake, volume * concentration
  
  #pellet and supplement is fed at ~10% food intake each by mass(Wood et al., 2020)
  
  p.intake = f.intake *(1 - f.h)
  
  f.pel ~ dbeta(alpha.pel, beta.pel) #Pellet vs supplement: about half-and-half
  
  f.h ~ dbeta(alpha.h, beta.h) #Hay ~80% by mass
  
  f.intake = f.intake.perc * Body.mass
  #daily food intake ~ 1.5% +- 0.2% body mass
  f.intake.perc ~ dnorm(0.015, 1/0.002^2)
  
  w.intake = w.intake.perc * Body.mass
  #daily water intake ~ 2.5% +- 0.5% body mass
  w.intake.perc ~ dnorm(0.025, 1/0.005^2)
  
  #water concentration and Sr ratio
  c.Sr.w = Sr.w.m * c.w.m #product of concentration and Sr
  
  c.w.m ~ dlnorm(conc.w.mean, 1/conc.w.sd^2)
  
  Sr.w.m ~ dnorm(Sr.w.mean, 1/Sr.w.sd^2)
  
  #Define parameters for beta distributions
  #~80% hay in diet
  alpha.h ~ dnorm(32, 1/0.4^2)
  
  beta.h ~ dnorm(8, 1/4^2)
  
  #~20% supplement plus pellets
  alpha.pel ~ dnorm(40, 1/1^2)
  
  beta.pel ~ dnorm(40, 1/1^2)
  
  #Ca/Sr absorption rate of hay is estimated to be 60 +- 5%
  diges.h ~ dnorm(0.6, 1/0.05^2)
  
  #Ca/Sr absorption rate of pellet is assumed to be slightly higher 70% (15% fiber)
  diges.p ~ dnorm(0.7, 1/0.05^2)
  
  #Ca/Sr absorption rate of supplement is assumed to be around the same as pellet at 70% (12% fiber)
  diges.s ~ dnorm(0.7, 1/0.05^2)
  
  #supplement concentration and Sr ratio
  c.sup.m ~ dlnorm(conc.sup.mean, 1/conc.sup.sd^2)
  
  Sr.sup.m ~ dnorm(Sr.sup.mean, 1/Sr.sup.sd^2)
  
  #pellet concentration and Sr ratio
  c.pel.m ~ dlnorm(conc.pel.mean, 1/conc.pel.sd^2)
  
  Sr.pel.m ~ dnorm(Sr.pel.mean, 1/Sr.pel.sd^2)
  
  #modeling mixing of m.feed kgs of hay in the feeding process
  for(i in 1:m.feed){
    
    c.hay.m[i] ~ dlnorm(conc.hay.mean, 1/conc.hay.sd^2)
    
    Sr.hay.m[i] ~ dnorm(Sr.hay.mean, 1/Sr.hay.sd^2)
    
    c.Sr.hay[i] = Sr.hay.m[i] * c.hay.m[i] #product of concentration and Sr
    
  }
  
  #body mass of misha, used to estimate daily food and water intake in kg
  Body.mass ~ dnorm(Body.mass.mean, 1/Body.mass.sd^2)
  Body.mass.mean <- 3000 # kg
  Body.mass.sd <- 250 # kg
}



