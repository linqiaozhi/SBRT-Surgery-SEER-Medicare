library(pci2s)
library(glmnet)
library(ahaz)
library(arsenal)
library(tidyr)
library(dplyr)
library(patchwork)
library(lubridate)
library(timereg)
library('ggtext')
library(ggplot2)
source('utilities.R')
source('two.step.variable.selection.R')
source('load.data.R')

# subset.names  <- list(      'all.gte.65', 'sens1')
nc_time_days = 90
vs  <-  F
nboot = 1000
library(foreach) 
library(doParallel)
subset.name  <-  'all.gte.65'
sm  <- load.data(1002, subset.name, nc_time_days = nc_time_days, W1s.global = T)
analysis.name  <- sprintf('s1.nctime%d_%s', nc_time_days, subset.name)
X1_all  <-  model.matrix(as.formula(sprintf('~ %s', paste(sm$X1s, collapse = '+'))),  sm$A.final)[,-1]
Z_all  <- sm$A.final[,sm$Zs]
A_ = (sm$A.final$tx == 'sbrt')*1.0

Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
cl  <- makeCluster(5, type = "FORK")
registerDoParallel(cl)
comb <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}
# sensitivity.outcomes  <- foreach(exclude.variable.i=1:length(colnames(Z_all)), .combine ='rbind', .multicombine = T,.packages=c('dplyr','timereg', 'pci2s'))  %dopar% {
system.time(sensitivity.outcomes  <- foreach(exclude.variable.i=1:length(colnames(Z_all)), .combine ='comb', .multicombine = T,.packages=c('dplyr','timereg', 'pci2s'))  %dopar% {
    exclude.variable  <-  colnames(Z_all)[exclude.variable.i]
    selected.columns <- list()
    for (W1i in 1:length(sm$W1s)) {
            outcome.name  <- sm$W1s[[W1i]]
            varnames  <-  c(colnames(X1_all), colnames(Z_all))
            selected.columns[[outcome.name]]  <- varnames[varnames != exclude.variable]
    }
        outcome.i  <- 2
        outcome.name  <-  sm$outcome.names[outcome.i]
        W1  <- sm$W1s[[outcome.i]]
        print(outcome.name)
        A.temp  <-  sm$A.final %>% mutate( 
                                          outcome.time = if_else (tt == 0, 0.5, tt)/ 365+ runif(dim(sm$A.final)[1] ,0,1)*1e-8,
                                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                                          W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
                                          cause = case_when (
                                                             nna(!!rlang::sym(outcome.name))~ 0, # primary event
                                                             nna(!!rlang::sym(W1))~ 1,
                                                             T ~ -1
                                          )
        )
        # Proximal
        set.seed(3)
            A_ = (A.temp$tx == 'sbrt')*1.0
            X1_  <-  model.matrix(as.formula(sprintf('~ %s', paste(sm$X1s, collapse = '+'))),  A.temp)[,-1]
            X1_selected  <-  X1_[,colnames(X1_)  %in% selected.columns[[W1]] ]
            X2_  <-  model.matrix(as.formula(sprintf('~ %s', paste(sm$X2s, collapse = '+'))),  A.temp)[,-1]
            Z_  <- A.temp[,sm$Zs]
            Z_selected  <-  Z_[,colnames(Z_)  %in% selected.columns[[W1]] ]
            Y_ = A.temp$outcome.time
            print(system.time(mout  <- p2sls.cprisk.nc  (times = Y_, cause = A.temp$cause, A = A_,  X1 = X1_selected,X2 = X2_,Z = Z_selected, nc_time = nc_time_days/365,  bootstrap = T, nboot = nboot, conf.level = 0.95) ))
            mout
            est <- mout$beta_a
            hd.outcome.proximal  <-  c(outcome=exclude.variable.i, estimate=est, low_ci=mout$beta_a_ci[1], high_ci=mout$beta_a_ci[2])
            print(hd.outcome.proximal)
           list(hd.outcome.proximal)
})
stopCluster(cl)
str(sensitivity.outcomes)
rownames(sensitivity.outcomes[[1]])  <- colnames(Z_all)
saveRDS(sensitivity.outcomes[[1]], 'data/sensitivity.outcomes2.RDS')

# 16 hours

