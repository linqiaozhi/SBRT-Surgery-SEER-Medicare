library(dplyr)
library(survival)

set.seed(3)
calc.err  <-  function (cox.out, beta_A)  100*(coef(cox.out)['A'] - beta_A)  / beta_A
logistic  <-  function(x) 1/(1+exp( x))

# Generate data
N  <-  5E4
U  <-   rnorm( N) # Unmeasured confounder
W  <-  rnorm (N) + U # proxy 1
Z  <-  rnorm (N) + U # proxy 2
X  <-   rnorm( N)
A  <-   rbinom( N, size =1, prob = logistic (0.3 + 0.4*U + X))
table( A, useNA="ifany")

beta_A  <-  0.5
x  <- A*beta_A + U + X
#x  <- A*beta_A
unif  <- runif( N ) 
lambda  <-  0.01
rateC  <-  0.01
#https://stats.stackexchange.com/questions/135124/how-to-create-a-toy-survival-time-to-event-data-with-right-censoring/135129#135129
T_  <-  - log (unif) / ( lambda * exp ( x ) ) 
summary(T_)
C_ <- rexp(n=N, rate=rateC)
outcome.time  <-  pmin( T_, C_)
outcome.observed  <-  T_ < C_
table( outcome.observed, useNA="ifany")
plot(survfit(Surv(outcome.time, outcome.observed)~A))

data.sim  <- data.frame ( outcome.time, outcome.observed, A, U, W, Z, X) 

## With knowledge of U
#cox.out  <-  coxph ( Surv(outcome.time, outcome.observed) ~ A + U + X , data = data.sim)
#print(sprintf('Adjusting for U, presuming it were measured: %.2f%% error', calc.err( cox.out, beta_A)))
#
#
## Without knowledge of U
#cox.out  <-  coxph ( Surv(outcome.time, outcome.observed) ~ A + X, data = data.sim)
#print(sprintf('No adjustment: %.2f%% error', calc.err( cox.out, beta_A)))
#
## Without knowledge of U, but knowing W and Z
#cox.out  <-  coxph ( Surv(outcome.time, outcome.observed) ~ A + W + X + Z, data = data.sim)
#print(sprintf('Adjusting directly for W and Z: %.2f%% error', calc.err( cox.out, beta_A)))
#
## Without knowledge of U
##stage1  <-  lm( W ~ A + Z + X, data = subset ( data.sim, outcome.observed == T))
#stage1  <-  lm( W ~ A + Z + X, data =  data.sim )
#W_hat  <-  predict(stage1, data.sim)
#cox.out  <-  coxph ( Surv(outcome.time, outcome.observed) ~ A+ X + W_hat, data = cbind(data.sim, W_hat))
#coef(cox.out)['A']
#print(sprintf('Two stage approach to adjustment: %.2f%% error', calc.err( cox.out, beta_A)))

################################
# Poisson 
################################

kfit1  <-  glm( outcome.observed ~ A + U + X  + offset(log(outcome.time)),
               family = poisson, 
               data = data.sim)
print(sprintf('Poisson with knowledge of U: %.2f%% error', calc.err( kfit1, beta_A)))

kfit  <-  glm( outcome.observed ~ A +  X  + offset(log(outcome.time)),
               family = poisson, 
               data = data.sim)
print(sprintf('Poisson without knowledge of U: %.2f%% error', calc.err( kfit, beta_A)))

kfit2  <-  glm( outcome.observed ~ A + X + W + Z  + offset(log(outcome.time)),
               family = poisson, 
               data = data.sim)
print(sprintf('Poisson without knowledge of U, but adusting for W and Z: %.3f%% error', calc.err( kfit2, beta_A)))

stage1  <-  lm( W ~ A + Z + X , data =  data.sim )
W_hat  <-  predict(stage1, data.sim)
data.sim.2  <-  data.frame( data.sim, W_hat = W_hat)
kfit4  <-  glm( outcome.observed ~ A +  X + W_hat +  offset(log(outcome.time)),
               family = poisson, 
               data = data.sim.2)
print(sprintf('Poisson with two stage: %.2f%% error', calc.err( kfit4, beta_A)))


################################
# Poisson with splits
################################

data.sim.split  <-  survSplit ( Surv(outcome.time, outcome.observed)~., cut = c(5,10,50), episode = 'interval', data = data.sim)
kfit1  <-  glm( outcome.observed ~ A + U + X + factor(interval) -1 + offset(log(outcome.time-tstart)),
               family = poisson, 
               data = data.sim.split)
summary(kfit1)
print(sprintf('Poisson with knowledge of U: %.2f%% error', calc.err( kfit1, beta_A)))

kfit2  <-  glm( outcome.observed ~ A +  X + factor(interval) -1 + offset(log(outcome.time-tstart)),
               family = poisson, 
               data = data.sim.split)
print(sprintf('Poisson without knowledge of U: %.2f%% error', calc.err( kfit2, beta_A)))


kfit3  <-  glm( outcome.observed ~ A +  X + W + Z+ factor(interval) -1 + offset(log(outcome.time-tstart)),
               family = poisson, 
               data = data.sim.split)
print(sprintf('Poisson without knowledge of U, but adusting for W and Z: %.3f%% error', calc.err( kfit3, beta_A)))



stage1  <-  lm( W ~ A + Z + X , data =  data.sim )
W_hat  <-  predict(stage1, data.sim)
data.sim.2  <-  data.frame( data.sim, W_hat = W_hat)
data.sim.2.split  <-  survSplit ( Surv(outcome.time, outcome.observed)~., cut = c(5,10,50), episode = 'interval', data = data.sim.2)
kfit4  <-  glm( outcome.observed ~ A +  X + W_hat +  factor(interval) -1 + offset(log(outcome.time-tstart)),
               family = poisson, 
               data = data.sim.2.split)
print(sprintf('Poisson with two stage: %.2f%% error', calc.err( kfit4, beta_A)))

