# devtools::load_all('/Users/george/Research_Local/pci2s_gcl/pci2s')
library(pci2s)
# devtools::install_github('KenLi93/pci2s')
# devtools::load_all('/Users/george/Research_Local/pci2s/')
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

subset.names  <- list(      'all.gte.65', 'sens1')
nc_time_days = 90
vs  <-  F
nboot = 1000
library(foreach) 
library(doParallel)
 # subset.name  <-  subset.names[[1]]
for (subset.name in subset.names ){
    # Each analysis is separate, so set the seed for reproducibility. If
    # it's set outside the loop then the results can only be reproduced by
    # rerunning the entire script.
    set.seed(3)
    sm  <- load.data(1003, subset.name, nc_time_days = nc_time_days, W1s.global = T)
    analysis.name  <- sprintf('data49.nctime%d_%s', nc_time_days, subset.name)
    print('====================')
    print(analysis.name)
    print('====================')
    #################################
    ## Table 1 
    #################################
    if (T) {
        sm$A.final$mediastinal.staging  <- sm$A.final$med_pre_count > 0 |  sm$A.final$ebus_pre_count > 0
        sm$A.final$histology3  <- sm$A.final$histology
        sm$A.final$histology3[sm$A.final$histology3 == 'Adenocarcinoma, papillary']  <- 'Adenocarcinoma, NOS/other'

        Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")
        tblcontrol <- tableby.control(numeric.stats = c('Nmiss', 'meansd'), numeric.simplify = T, cat.simplify =T, digits = 1,total = T,test = F)
        f  <-  sprintf( 'tx ~ %s', paste( c(sm$X.numeric, sm$X.factor,'race', 'histology3','mediastinal.staging', gsub('_count', '_count_bool', c(sm$Zs))), collapse = "+"))
        labels(sm$A.final)  <-  label_list
        tt <- tableby(as.formula(f), data=sm$A.final, control = tblcontrol)
        summary(tt) %>% write2html(sprintf('/Users/george/Research_Local/SEER-Medicare/tbls/table1_2_%s.htm', analysis.name), quiet=T)
    }


    #################################
    # # variable selection for W1
    ################################
    X1_all  <-  model.matrix(as.formula(sprintf('~ %s', paste(sm$X1s, collapse = '+'))),  sm$A.final)[,-1]
    Z_all  <- sm$A.final[,sm$Zs]
    A_ = (sm$A.final$tx == 'sbrt')*1.0
    selected.columns <- list()
    for (W1i in 1:length(sm$W1s)) {
        outcome.name  <- sm$W1s[[W1i]]
        if (vs ) {
            # design.mat  <-  as.matrix(cbind(A_, X1_all, Z_all))
            # A.temp  <-  sm$A.final %>% mutate( 
            #                                   outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
            #                                   outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
            #                                   outcome.time  = if_else (outcome.time == 0, 0.5, outcome.time)/ 365+ runif(dim(sm$A.final)[1] ,0,1)*1e-8
            # ) 
            # mm  <-   tune.ahazpen( Surv(A.temp$outcome.time, A.temp$outcome.bool), design.mat )
            # selected.columns[[outcome.name]]  <- get.selected.columns.ahaz(mm, s = 'lambda.min', colnames(design.mat), verbose=T, min.vars = 1)
            # print(sprintf( 'Included variables in stage 1 for %s: %s', outcome.name, paste( selected.columns[[outcome.name]], collapse = ', ')))
        }else {
            selected.columns[[outcome.name]]  <- c(colnames(X1_all), colnames(Z_all))
        }
    }


    ################################
    # Obtain estimates using the selected variables 
    ################################
    # Time-to-event outcomes 
    hazard.differences.outcomes  <-  make.odds.ratio.df ( sm$outcome.names) 
    hazard.differences.outcomes.adj  <-  make.odds.ratio.df ( sm$outcome.names) 
    hazard.differences.outcomes.proximal  <-  make.odds.ratio.df ( sm$outcome.names) 
    mouts  <-  list()
    
    cl  <- makeCluster(10)
    registerDoParallel(cl)

    comb <- function(...) {
      mapply('rbind', ..., SIMPLIFY=FALSE)
    }
    # outcome.i = 2
    hd.outcomes  <- foreach(outcome.i=1:length(sm$outcome.names), .combine ='comb', .multicombine = T,.packages=c('dplyr','timereg', 'pci2s') ) %dopar% {

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
        # Raw
        f  <- as.formula('Surv(outcome.time, outcome.bool) ~ const(tx)')
        m  <-  aalen( f ,  data = A.temp, robust = 0)
        hd.outcome  <-  c(outcome=outcome.i, coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
        # Adjust for X
        f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + %s',  paste(sprintf('const(%s)', c(sm$X2s, sm$Zs)), collapse="+") )
        m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
        hd.outcome.adj  <-  c(outcome=outcome.i,  coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
        print.martingales  <-  T
        # Proximal
        set.seed(3)
        if (!outcome.name %in% c('death', 'death.other.cause')) {
            A_ = (A.temp$tx == 'sbrt')*1.0
            X1_  <-  model.matrix(as.formula(sprintf('~ %s', paste(sm$X1s, collapse = '+'))),  A.temp)[,-1]
            X1_selected  <-  X1_[,colnames(X1_)  %in% selected.columns[[W1]] ]
            X2_  <-  model.matrix(as.formula(sprintf('~ %s', paste(sm$X2s, collapse = '+'))),  A.temp)[,-1]
            Z_  <- A.temp[,sm$Zs]
            Z_selected  <-  Z_[,colnames(Z_)  %in% selected.columns[[W1]] ]
            Y_ = A.temp$outcome.time
            table( A.temp$cause, useNA="ifany")
            print(system.time(mout  <- p2sls.cprisk.nc  (times = Y_, cause = A.temp$cause, A = A_,  X1 = X1_selected,X2 = X2_,Z = Z_selected, nc_time = nc_time_days/365,  bootstrap = T, nboot = nboot, conf.level = 0.95) ))
            if (print.martingales && outcome.name == 'death.lc.specific') {
                res  <- p2sls.cprisk.nc.martingale.residuals(times = Y_, cause = A.temp$cause, A = A_,  a = 1, X=  X1_selected, Z = as.matrix(Z_selected), nc_time = nc_time_days/365  )
                saveRDS(res, sprintf('data/%s.%s.martingale.residuals.rds', analysis.name, outcome.name))
                # A.temp %>% filter (PVD_pre_12months_unique_count > 12) %>% t
                # head(Z_selected)
                # max of each column of Z_selected
                # apply(Z_selected, 2, max)
            }
            mouts[[outcome.name]]  <-  mout
            est <- mout$beta_a
            hd.outcome.proximal  <-  c(outcome=outcome.i, estimate=est, low_ci=mout$beta_a_ci[1], high_ci=mout$beta_a_ci[2])
        }else {
            print('Proximal causal estimates not defined for this outcome')
            hd.outcome.proximal  <-  c(outcome=outcome.i, NA, NA, NA)
        }
        list(hd.outcome, hd.outcome.adj, hd.outcome.proximal)
    }
    # order each element of list by outcome
    # hd.outcomes2  <- laply(hd.outcomes, function(x) x[order(x[,1]),])
    # str(hd.outcomes)
    # join this
    hazard.differences.outcomes[,1:3]  <- hd.outcomes[[1]][,-1]
    hazard.differences.outcomes.adj[,1:3]  <- hd.outcomes[[2]][,-1]
    hazard.differences.outcomes.proximal[,1:3]  <- hd.outcomes[[3]][,-1]
    # print(selected.columns)
    print(hazard.differences.outcomes[,1:3])
    print(hazard.differences.outcomes.adj[,1:3])
    print(hazard.differences.outcomes.proximal[,1:3])
    stopCluster(cl)
    saveRDS(hazard.differences.outcomes, sprintf('data/%s.hazard.differences.outcomes.raw.rds', analysis.name))
    saveRDS(hazard.differences.outcomes.adj, sprintf('data/%s.hazard.differences.outcomes.adj.rds', analysis.name))
    saveRDS(hazard.differences.outcomes.proximal, sprintf('data/%s.hazard.differences.outcomes.proximal.rds', analysis.name))
}


