library(dplyr)
library(survival)

calc.err  <-  function (cox.out, beta_A)  100*(coef(cox.out)['A'] - beta_A)  / beta_A
logistic  <-  function(x) 1/(1+exp( -x))

sim  <-  function( N = 1E4, eps = 1,lambda = 0.01, rateC = 0.01, beta_A=0.3,  verbose = F , plotit=F)  { 
    # Generate data
    U  <-   rnorm( N)           # Unmeasured confounder
    W  <-  U + eps*rnorm (N)    # proxy 1
    Z  <-  U + eps*rnorm (N)    # proxy 2
    X  <-   rnorm( N)           # Measured confounder
    A  <-   rbinom( N, size =1, prob = logistic (0.3 + U + X )) # treatment
    if (verbose) print(table( A, useNA="ifany"))
    x  <- A*beta_A + U + X
    unif  <- runif( N )         
    T_  <-  - log (unif) / ( lambda * exp ( x ) )       # Generating survival times with constant baseline hazard, e.g. as https://stats.stackexchange.com/a/135129/103007
    if (rateC == -1 ) {
        C_  <- rep(Inf,N) # No censoring
    }else{
        C_ <- rexp(n=N, rate=rateC)     # Censoring times drawn from Exponential distribution
    }
    outcome.time  <-  pmin( T_, C_)
    outcome.observed  <-  T_ < C_
    if (verbose) print(table( outcome.observed, useNA="ifany"))
    if (plotit) plot(survfit(Surv(outcome.time, outcome.observed)~A))
    data.sim  <- data.frame ( outcome.time, outcome.observed, A, U, W, Z ) 
    out  <-  list()
    kfit  <-  glm( outcome.observed ~ A + X + U +  offset(log(outcome.time)),
                   family = poisson, 
                   data = data.sim)
    out[['with.u']]  <-  calc.err( kfit, beta_A)
    if (verbose) print(summary(kfit))
    kfit  <-  glm( outcome.observed ~ A + X +  offset(log(outcome.time)),
                   family = poisson, 
                   data = data.sim)
    out[['without.u']]  <-  calc.err( kfit, beta_A)
    kfit  <-  glm( outcome.observed ~ A + X +  W + Z  + offset(log(outcome.time)),
                   family = poisson, 
                   data = data.sim)
    out[['without.u.with.w.z']]  <-  calc.err( kfit, beta_A)
    stage1  <-  lm( W ~ A + X+ Z  , data =  data.sim )
    W_hat  <-  predict(stage1, data.sim)
    data.sim.2  <-  data.frame( data.sim, W_hat = W_hat)
    kfit  <-  glm( outcome.observed ~ A + X+  W_hat +  offset(log(outcome.time)),
                   family = poisson, 
                   data = data.sim.2)
    if (verbose) print(summary(kfit))
    out[['without.u.two.stage']]  <-  calc.err( kfit, beta_A)
    kfit  <-  coxph( Surv(outcome.time, outcome.observed) ~ A + X + U ,  data = data.sim)
    if (verbose) print(summary(kfit))
    out[['with.u.cox']]  <-  calc.err( kfit, beta_A)
    kfit  <-  coxph( Surv(outcome.time, outcome.observed) ~ A + X + W_hat ,  data = data.sim.2)
    if (verbose) print(summary(kfit))
    out[['without.u.two.stage.cox']]  <-  calc.err( kfit, beta_A)
    return(out)
}

print.results  <-  function ( outs) {
    out.means  <-  apply( outs , 2, mean)
    out.se  <-  apply( outs , 2, FUN = function(x) 1.96*sd(x)/sqrt(length(x)))
    print(sprintf('Poisson with knowledge of U: %.1f%% +/- %.1f%% error', out.means['with.u'], out.se['with.u']))
    print(sprintf('Poisson without knowledge of U: %.1f%% +/- %.1f%% error', out.means['without.u'], out.se['without.u']))
    print(sprintf('Poisson without knowledge of U, but adjusting for W and Z: %.1f%%+/- %.1f%%  error', out.means['without.u.with.w.z'], out.se['without.u.with.w.z']))
    print(sprintf('Poisson with two stage: %.1f%% +/- %.1f%% error', out.means['without.u.two.stage'], out.se['without.u.two.stage']))
    print(sprintf('Cox with knowledge of U: %.1f%% +/- %.1f%% error', out.means['with.u.cox'], out.se['with.u.cox']))
    print(sprintf('Cox with two stage: %.1f%% +/- %.1f%% error', out.means['without.u.two.stage.cox'], out.se['without.u.two.stage.cox']))
}


################################
# Simulations 
################################
pdf ('curves.pdf')
sim(N= 1E4, eps = 1,  rateC = 0.01,beta_A=0.3,  verbose=T, plotit=F)
dev.off ()

Nits  <- 200
set.seed(5)
outs  <-  data.frame () 
for (i in 1:Nits) {
    outs  <-  rbind ( outs, as.data.frame(sim(N= 1E4, eps = 1,  rateC = -1,  verbose=F))) 
}
print('Medium noise regime (eps=1), no censoring')
print.results(outs)


set.seed(5)
outs  <-  data.frame () 
for (i in 1:Nits) {
    outs  <-  rbind ( outs, as.data.frame(sim(N= 1E4, eps = 3,  rateC = -1,  verbose=F))) 
}
print('High noise regime (eps=3), no censoring')
print.results(outs)


set.seed(5)
outs  <-  data.frame () 
for (i in 1:Nits) {
    outs  <-  rbind ( outs, as.data.frame(sim(N= 1E4, eps = 1))) 
}
print('Medium noise regime, eps=1, with censoring')
print.results(outs)

ggplot2::qplot (without.u.two.stage.cox, data = outs, geom = 'histogram')


set.seed(5)
outs  <-  data.frame () 
for (i in 1:Nits) {
    outs  <-  rbind ( outs, as.data.frame(sim(N= 1E4, eps = 3,  verbose =F))) 
}
print('High noise regime, eps=3, with censoring')
print.results(outs)

