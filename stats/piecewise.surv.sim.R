################################
# Piecewise constant
################################
library(survival)
library(dplyr)
# https://stats.stackexchange.com/questions/105881/how-to-simulate-survival-times-using-true-base-line-hazard-function
# https://stats.stackexchange.com/questions/135124/how-to-create-a-toy-survival-time-to-event-data-with-right-censoring/135129#135129
# https://grodri.github.io/glms/notes/c7s4
# https://cran.r-project.org/web/packages/survival/vignettes/approximate.pdf
# https://stats.stackexchange.com/questions/359970/checking-the-proportional-hazard-assumption/400981#400981
# https://jamanetwork-com.ezp-prod1.hul.harvard.edu/journals/jama/fullarticle/2763185
# https://web.stanford.edu/~lutian/coursepdf/STAT331unit3.pdf
tdom <- seq(0, 10, by=0.001)
haz <- rep(0, length(tdom))
haz[tdom <= 1] <- 0.5
haz[tdom > 1 & tdom <= 2.5] <- 0.1
haz[tdom > 2.5] <- 0.3
cumhaz <- cumsum(haz*0.001) # This is integral of lambda from 0 to t
S <- exp(-cumhaz)
# par(mfrow=c(3,1))
# plot(tdom, haz, type='l', xlab='Time domain', ylab='Hazard')
# plot(tdom, cumhaz, type='l', xlab='Time domain', ylab='Cumulative hazard')
# plot(tdom, S, type='l', xlab='Time domain', ylab='Survival')

N  <-  10000
u <- runif(N) *0.99 # Avoid T = 0 
fofo  <-  colSums(outer(S, u, `>`))
fofo[fofo==0] = 1
T_ <- tdom[fofo]
C_ <- rexp(n=N, rate=0.0001)
outcome.time  <-  pmin( T_, C_)
outcome.observed  <-  T_ < C_
table( outcome.observed, useNA="ifany") # 

#sum(outcome.time[outcome.observed] > 3)
if (min(outcome.time ) ==0 ) stop('T must be > 0')

data.sim  <-  data.frame ( outcome.time=outcome.time+0.0001, outcome.observed)
summary(data.sim$outcome.time)

# dev.off()
 plot(survfit(Surv(outcome.time, outcome.observed)~1))


data.sim.split  <-  survSplit ( Surv(outcome.time, outcome.observed)~., cut = c(1,2.5), episode = 'interval', data = data.sim)
kfit1  <-  glm( outcome.observed ~ factor(interval) -1 + offset(log(outcome.time-tstart)),
               family = poisson, 
               data = data.sim.split)
summary(kfit1)
exp(coef(kfit1))


################################
# Piecewise constant with covariates 
################################
rm( list = ls())
library(survival)
library(dplyr)
# https://stats.stackexchange.com/questions/105881/how-to-simulate-survival-times-using-true-base-line-hazard-function
delta  <- 0.001
tdom <- seq(0, 20, by=delta)
haz0 <- rep(0, length(tdom))
haz0[tdom <= 1] <- 0.3
haz0[tdom > 1 & tdom <= 2.5] <- 0.1
haz0[tdom > 2.5] <- 0.2

N  <-  5000
logistic  <-  function(x) 1/(1+exp( x))
U  <-   rnorm( N) # Unmeasured confounder
W  <-  rnorm (N) + U # proxy 1
Z  <-  rnorm (N) + U # proxy 2
X  <-   rnorm( N)
A  <-   rbinom( N, size =1, prob = logistic (0.3 + 0.4*U + X))

beta_ = 0.5
ee  <-  exp( beta_*A + 0.3*U)
#ee  <-  exp( beta_*A )
haz_i  <- (ee) %*% t(haz0) # Each row is an individual S_i
cumhaz_i  <-  t(apply( haz_i*delta, 1, FUN=cumsum))
S_i <- exp(-cumhaz_i)
# par(mfrow=c(3,1))
# plot(tdom, haz, type='l', xlab='Time domain', ylab='Hazard')
# plot(tdom, cumhaz, type='l', xlab='Time domain', ylab='Cumulative hazard')
# plot(tdom, S, type='l', xlab='Time domain', ylab='Survival')

u <- runif(N) *0.99 # Avoid T = 0 
T_  <- rep(-1, N)
for (i in 1:N) {
    fofo  <-  sum( S_i[i,] > u[i] )
    T_[i]  <-  tdom[fofo]
}

C_ <- rexp(n=N, rate=0.01)
outcome.time  <-  pmin( T_, C_)
outcome.observed  <-  T_ < C_
table( outcome.observed, useNA="ifany") # 

#sum(outcome.time[outcome.observed] > 3)
if (min(outcome.time ) ==0 ) stop('T must be > 0')

data.sim  <-  data.frame ( outcome.time=outcome.time+0.0001, outcome.observed, A, X, U, W,Z)
summary(data.sim$outcome.time)

# dev.off()
 plot(survfit(Surv(outcome.time, outcome.observed)~A))


data.sim.split  <-  survSplit ( Surv(outcome.time, outcome.observed)~., cut = c(1,2.5), episode = 'interval', data = data.sim)
kfit1  <-  glm( outcome.observed ~ A+ U + factor(interval) -1 + offset(log(outcome.time-tstart)),
               family = poisson, 
               data = data.sim.split)
summary(kfit1)
exp(coef(kfit1))


kfit2  <-  coxph( Surv(outcome.time, outcome.observed) ~ A+ U ,  data = data.sim.split)
summary(kfit2)
(coef(kfit2))

