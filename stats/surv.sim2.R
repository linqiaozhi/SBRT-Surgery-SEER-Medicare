library(dplyr)
library(survival)

calc.err  <-  function (cox.out, beta_A)  100*(coef(cox.out)['A'] - beta_A)  / beta_A
logistic  <-  function(x) 1/(1+exp( x))

sim  <-  function( N = 1E4, eps = 1,lambda = 0.01, rateC = 0.01, beta_A=1,  verbose = F )  { 
    # Generate data
    U  <-   rnorm( N)           # Unmeasured confounder
    W  <-  U + eps*rnorm (N)    # proxy 1
    Z  <-  U + eps*rnorm (N)    # proxy 2
    A  <-   rbinom( N, size =1, prob = logistic (0.3 + U ))
    x  <- A*beta_A + U 
    unif  <- runif( N ) 
    T_  <-  - log (unif) / ( lambda * exp ( x ) ) 
    C_ <- rexp(n=N, rate=rateC)
    outcome.time  <-  pmin( T_, C_)
    outcome.observed  <-  T_ < C_
    if (verbose) print(table( outcome.observed, useNA="ifany"))
    if (verbose ) print(table( outcome.observed, useNA="ifany"))
    if (verbose) plot(survfit(Surv(outcome.time, outcome.observed)~A))
    data.sim  <- data.frame ( outcome.time, outcome.observed, A, U, W, Z ) 
    ################################
    # Poisson 
    ################################
    out  <-  list()
    kfit  <-  glm( outcome.observed ~ A + U +  offset(log(outcome.time)),
                   family = poisson, 
                   data = data.sim)
    out[['with.u']]  <-  calc.err( kfit, beta_A)
    kfit  <-  glm( outcome.observed ~ A +   offset(log(outcome.time)),
                   family = poisson, 
                   data = data.sim)
    out[['without.u']]  <-  calc.err( kfit, beta_A)
    kfit  <-  glm( outcome.observed ~ A +  W + Z  + offset(log(outcome.time)),
                   family = poisson, 
                   data = data.sim)
    out[['without.u.with.w.z']]  <-  calc.err( kfit, beta_A)
    stage1  <-  lm( W ~ A + Z  , data =  data.sim )
    W_hat  <-  predict(stage1, data.sim)
    data.sim.2  <-  data.frame( data.sim, W_hat = W_hat)
    kfit  <-  glm( outcome.observed ~ A +   W_hat +  offset(log(outcome.time)),
                   family = poisson, 
                   data = data.sim.2)
    out[['without.u.two.stage']]  <-  calc.err( kfit, beta_A)
    return(out)
}

print.results  <-  function ( outs) {
    out.means  <-  apply( outs , 2, mean)
    print(sprintf('Poisson with knowledge of U: %.2f%% error', out.means['with.u']))
    print(sprintf('Poisson without knowledge of U: %.2f%% error', out.means['without.u']))
    print(sprintf('Poisson without knowledge of U, but adjusting for W and Z: %.3f%% error', out.means['without.u.with.w.z']))
    print(sprintf('Poisson with two stage: %.2f%% error', out.means['without.u.two.stage']))
}


set.seed(5)
outs  <-  data.frame () 
for (i in 1:10) {
    outs  <-  rbind ( outs, as.data.frame(sim(N= 1E4, eps = 1))) 
}
print('Medium noise regime, eps=1')
print.results(outs)


set.seed(5)
outs  <-  data.frame () 
for (i in 1:10) {
    outs  <-  rbind ( outs, as.data.frame(sim(N= 1E4, eps = 3,  verbose =F))) 
}
print('High noise regime, eps=3')
print.results(outs)
