library(dplyr)
library(survival)

calc.err  <-  function (cox.out, beta_A)  100*(coef(cox.out)['A'] - beta_A)  / beta_A
logistic  <-  function(x) 1/(1+exp( x))


# Generate data
N  <-  1E4
U  <-   rnorm( N) # Unmeasured confounder
W  <-  rnorm (N) + U # proxy 1
Z  <-  rnorm (N) + U # proxy 2
X  <-  rep( 0, N)
X  <-   rnorm( N)
A  <-   rbinom( N, size =1, prob = logistic (0.3 + 0.4*U + X))
table( A, useNA="ifany")
#U <-   rbinom( N, size =1, prob = 0.5) # Unmeasured confounder
# W  <-   rbinom( N, size =1, prob = 0.2*U+0.2)
# Z  <-   rbinom( N, size =1, prob = 0.2*U+0.2)

beta_A  <-  1.2
x  <- A*beta_A + U*2 + X
unif  <- runif( N ) 
lambda  <-  0.01
rho  <-  1
rateC  <-  0.1
#https://stats.stackexchange.com/questions/135124/how-to-create-a-toy-survival-time-to-event-data-with-right-censoring/135129#135129


piecewise.T  <-  function ( unif, x, cut.point = 5, lambda1 = 0.1, lambda2=0.5) {
    out  <-  -log (unif)/ ( lambda1 * exp(x) ) 
    out [ unif > cut.point ] =  out [ unif > cut.point ] -cut.point * ( lambda1 - lambda2)
    return(out)
}

#T_  <-  - log (unif) / ( lambda * exp ( x ) ) ^ (1 / rho ) 
T_  <-  piecewise.T (unif, x, cut.point = 0.3) 
summary(T_)
C_ <- rexp(n=N, rate=rateC)
outcome.time  <-  pmin( T_, C_)
outcome.observed  <-  T_ < C_
table( outcome.observed, useNA="ifany")
summary(outcome.time)


data.sim  <- data.frame ( outcome.time, outcome.observed, A, U, W, Z, X) 

# With knowledge of U
cox.out  <-  coxph ( Surv(outcome.time, outcome.observed) ~ A + U + X , data = data.sim)
print(sprintf('Adjusting for U, presuming it were measured: %.2f%% error', calc.err( cox.out, beta_A)))


# Without knowledge of U
cox.out  <-  coxph ( Surv(outcome.time, outcome.observed) ~ A + X, data = data.sim)
print(sprintf('No adjustment: %.2f%% error', calc.err( cox.out, beta_A)))

# Without knowledge of U, but knowing W and Z
cox.out  <-  coxph ( Surv(outcome.time, outcome.observed) ~ A + W + X + Z, data = data.sim)
print(sprintf('Adjusting directly for W and Z: %.2f%% error', calc.err( cox.out, beta_A)))

# Without knowledge of U
#stage1  <-  lm( W ~ A + Z + X, data = subset ( data.sim, outcome.observed == T))
stage1  <-  lm( W ~ A + Z + X, data =  data.sim )
W_hat  <-  predict(stage1, data.sim)
cox.out  <-  coxph ( Surv(outcome.time, outcome.observed) ~ A+ X + W_hat, data = cbind(data.sim, W_hat))
coef(cox.out)['A']
print(sprintf('Two stage approach to adjustment: %.2f%% error', calc.err( cox.out, beta_A)))


################################
# Poisson 
################################

data.sim.split  <-  survSplit ( Surv(outcome.time, outcome.observed)~., cut = c(5,10,50), episode = 'interval', data = data.sim)
#kfit1  <-  glm( outcome.observed ~ A + U + X + factor(interval) -1 + offset(log(outcome.time-tstart)),
kfit1  <-  glm( outcome.observed ~ A + U + X +  offset(log(outcome.time-tstart)),
               family = poisson, 
               data = data.sim.split)
summary(kfit1)
coef(cox.out)['A']
coef(kfit1)['A']


table( data.sim.split$outcome.observed, useNA="ifany")

utime

