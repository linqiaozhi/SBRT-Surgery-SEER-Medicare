library(dplyr)
library(survival)
library(timereg)

calc.err  <-  function (cox.out, beta_A, var.name = 'const(A)')  100*(coef(cox.out)[var.name,'Coef.'] - beta_A)  / beta_A
logistic  <-  function(x) 1/(1+exp( -x))

sim  <-  function( N = 1E5, eps = 1,rateC = 0.01, lambda = 0.015, beta_A = 0.001, beta_U = 0.001,  verbose = F , plotit=F)  { 
 # N = 1E5
# eps = 1
# lambda = 0.015 
# rateC = 0.01
# beta_A=0.001
# beta_U = 0.001
# verbose = T 
# plotit=T
    # Generate data
    U  <-   rnorm( N)+3           # Unmeasured confounder
    W  <-  U + eps*rnorm (N)    # proxy 1
    Z  <-  U + eps*rnorm (N)    # proxy 2
    X  <-   runif(N)           # Measured confounder
    A  <-   rbinom( N, size =1, prob = logistic (-3 + U + X )) # treatment
    if (verbose) print(table( A, useNA="ifany"))
    x  <- lambda + A*beta_A +beta_U*U + 0.001*X
    if (min(x) <0 ) stop('Negative hazard')
    unif  <- runif( N )         
    T_  <-  (- log (unif) / ( x ))       # Generating survival times with constant baseline hazard
    if (verbose) print(summary(T_))
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
    kfit  <-  aalen( Surv(outcome.time, outcome.observed) ~ 1+ const(A) + const(X) + const(U) ,  data = data.sim, robust = 0)
    if (verbose) print(summary(kfit))
    out[['with.u']]  <-  calc.err( kfit, beta_A)
    kfit  <-  aalen( Surv(outcome.time, outcome.observed) ~ 1+ const(A) + const(X)  ,  data = data.sim, robust = 0)
    if (verbose) print(summary(kfit))
    out[['without.u']]  <-  calc.err( kfit, beta_A)
    kfit  <-  aalen( Surv(outcome.time, outcome.observed) ~ 1+ const(A) + const(X) + const(W) + const(Z) ,  data = data.sim, robust = 0)
    if (verbose) print(summary(kfit))
    out[['without.u.with.w.z']]  <-  calc.err( kfit, beta_A)
    stage1  <-  lm( W ~ A + X+ Z  , data =  data.sim )
    W_hat  <-  predict(stage1, data.sim)
    data.sim.2  <-  data.frame( data.sim, W_hat = W_hat)
    kfit  <-  aalen( Surv(outcome.time, outcome.observed) ~ 1+ const(A) + const(X) + const(W_hat)  ,  data = data.sim.2, robust = 0)
    if (verbose) print(summary(kfit))
    out[['without.u.two.stage']]  <-  calc.err( kfit, beta_A)
    return(out)
}

print.results  <-  function ( outs) {
    out.means  <-  apply( outs , 2, mean)
    out.se  <-  apply( outs , 2, FUN = function(x) 1.96*sd(x)/sqrt(length(x)))
    print(sprintf('Aalen\'s with knowledge of U: %.1f%% +/- %.1f%% error', out.means['with.u'], out.se['with.u']))
    print(sprintf('Aalen\'s without knowledge of U: %.1f%% +/- %.1f%% error', out.means['without.u'], out.se['without.u']))
    print(sprintf('Aalen\'s without knowledge of U, but adjusting for W and Z: %.1f%%+/- %.1f%%  error', out.means['without.u.with.w.z'], out.se['without.u.with.w.z']))
    print(sprintf('Aalen\'s with two stage: %.1f%% +/- %.1f%% error', out.means['without.u.two.stage'], out.se['without.u.two.stage']))
}


################################
# Simulations 
################################


pdf ('curves.pdf')
sim(N= 1E4, eps = 1,  rateC = -1, lambda = 0.015, beta_A = 0.001, beta_U = 0.001,verbose=T, plotit=T)
dev.off()

Nits  <- 1000
set.seed(5)
outs  <-  data.frame () 
for (i in 1:Nits) {
    outs  <-  rbind ( outs, as.data.frame(sim(N= 1E4, eps = 1,  rateC = 0.00001,  verbose=F))) 
}
print('Medium noise regime (eps=1), no censoring')
print.results(outs)


#dev.off ()
set.seed(5)
outs  <-  data.frame () 
for (i in 1:Nits) {
    outs  <-  rbind ( outs, as.data.frame(sim(N= 1E4, eps = 1, verbose =F))) 
}
print('Medium noise regime, eps=1, with censoring')
print.results(outs)



