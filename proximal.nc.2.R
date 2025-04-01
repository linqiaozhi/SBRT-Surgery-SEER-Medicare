# library(pci2s)
devtools::load_all('/Users/george/Research_Local/pci2s_gcl/pci2s')
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
for (subset.name in subset.names ){
    # Each analysis is separate, so set the seed for reproducibility. If
    # it's set outside the loop then the results can only be reproduced by
    # rerunning the entire script.
    set.seed(3)
    sm  <- load.data(33, subset.name, nc_time_days = nc_time_days)
    analysis.name  <- sprintf('v44.tte.nc_time_ss_no_vs_%d_%s', nc_time_days, subset.name)
    print('====================')
    print(analysis.name)
    print('====================')
    #################################
    ## Table 1 
    #################################
    if (T) {
        Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")
        tblcontrol <- tableby.control(numeric.stats = c('Nmiss', 'meansd'), numeric.simplify = T, cat.simplify =T, digits = 1,total = T,test = F)
        f  <-  sprintf( 'tx ~ %s', paste( c(sm$X.numeric, sm$X.factor,'histology', gsub('_count', '_count_bool', c(sm$Zs))), collapse = "+"))
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
            design.mat  <-  as.matrix(cbind(A_, X1_all, Z_all))
            A.temp  <-  sm$A.final %>% mutate( 
                                              outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                                              outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                                              outcome.time  = if_else (outcome.time == 0, 0.5, outcome.time)/ 365+ runif(dim(sm$A.final)[1] ,0,1)*1e-8
            ) 
            mm  <-   tune.ahazpen( Surv(A.temp$outcome.time, A.temp$outcome.bool), design.mat )
            selected.columns[[outcome.name]]  <- get.selected.columns.ahaz(mm, s = 'lambda.min', colnames(design.mat), verbose=T, min.vars = 1)
            print(sprintf( 'Included variables in stage 1 for %s: %s', outcome.name, paste( selected.columns[[outcome.name]], collapse = ', ')))
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
    for (outcome.i in 1:length(sm$outcome.names)){ 
        outcome.name  <-  sm$outcome.names[outcome.i]
        W1  <- sm$W1s[[outcome.i]]
        print(outcome.name)
        A.temp  <-  sm$A.final %>% mutate( 
                                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                                          outcome.time  = if_else (outcome.time == 0, 0.5, outcome.time)/ 365+ runif(dim(sm$A.final)[1] ,0,1)*1e-8,
                                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                                          W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt),
                                          W1.time  = if_else (W1.time == 0, 0.5, W1.time)/365,
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
        hazard.differences.outcomes[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
        # Adjust for X
        f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + %s',  paste(sprintf('const(%s)', c(sm$X2s, sm$Zs)), collapse="+") )
        m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
        hazard.differences.outcomes.adj[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
        # Proximal
        if (!outcome.name %in% c('death', 'death.other.cause')) {
            A_ = (A.temp$tx == 'sbrt')*1.0
            X1_  <-  model.matrix(as.formula(sprintf('~ %s', paste(sm$X1s, collapse = '+'))),  A.temp)[,-1]
            X1_selected  <-  X1_[,colnames(X1_)  %in% selected.columns[[W1]] ]
            X2_  <-  model.matrix(as.formula(sprintf('~ %s', paste(sm$X2s, collapse = '+'))),  A.temp)[,-1]
            Z_  <- A.temp[,sm$Zs]
            Z_selected  <-  Z_[,colnames(Z_)  %in% selected.columns[[W1]] ]
            Y_ = A.temp$outcome.time
            print(system.time(mout  <- p2sls.cprisk.nc  (times = Y_, cause = A.temp$cause, A = A_,  X1 = X1_selected,X2 = X2_,Z = Z_selected, nc_time = nc_time_days/365,  bootstrap = T, nboot = nboot, conf.level = 0.95) ))
            mouts[[outcome.name]]  <-  mout
            est <- mout$beta_a
            hazard.differences.outcomes.proximal[outcome.i,1:3]  <-  c(est, mout$beta_a_ci[1], mout$beta_a_ci[2])
        }else {
            print('Proximal causal estimates not defined for this outcome')
        }
        print(hazard.differences.outcomes[outcome.i,1:3])
        print(hazard.differences.outcomes.adj[outcome.i,1:3])
        print(hazard.differences.outcomes.proximal[outcome.i,1:3])
    } 
    saveRDS(hazard.differences.outcomes, sprintf('data/%s.hazard.differences.outcomes.raw.rds', analysis.name))
    saveRDS(hazard.differences.outcomes.adj, sprintf('data/%s.hazard.differences.outcomes.adj.rds', analysis.name))
    saveRDS(hazard.differences.outcomes.proximal, sprintf('data/%s.hazard.differences.outcomes.proximal.rds', analysis.name))
}

A.temp %>% count(outcome.bool, W1.bool,cause)
A.temp %>% filter( cause == 1) %>% select(outcome.time, W1.time, tt, outcome.bool, W1.bool, COD_TO_SITE_RECODE) %>% left_join(cod.df, by = c('COD_TO_SITE_RECODE'='COD_TO_SITE_RECODE')) %>% count(COD_TO_SITE_RECODE) %>% arrange( -n)

A.temp %>% filter( cause == 1) %>% select(outcome.time, W1.time, tt, outcome.bool, W1.bool, COD_TO_SITE_RECODE) %>% left_join(cod.df, by = c('COD_TO_SITE_RECODE'='COD_TO_SITE_RECODE')) %>% count(COD_TO_SITE_RECODE) %>% arrange( -n) %>% print(n=100)


# library(ggsurvfit)
# library(survival)
# library(tidycmprsk)
# dframe  <-  data.frame( time = A.temp$tt/365, status = A.temp$cause, A=A.temp$tx )
# dframe$A  <-  case_match(dframe$A,
#                          'sbrt' ~ 'SBRT',
#                          'sublobar' ~ 'Surgery')
# dframe$status  <-  case_match(
#                               dframe$status,
#                               0 ~ 'Cause-specific death',
#                               1 ~ 'Other-cause death',
#                               -1 ~ 'Censored'
#                               ) %>% as.factor
# dframe$status  <- factor(dframe$status, levels=c('Censored', 'Cause-specific death', 'Other-cause death'))
# cuminc.obj <- cuminc(Surv(time, status) ~ A, data = dframe)
# g2   <- ggcuminc(cuminc.obj, outcome = "Cause-specific death", linewidth = 1)  + ggtitle ("Cause-specific death")
# g2
# g3  <-  ggcuminc(cuminc.obj, outcome = "Other-cause death", linewidth = 1)+ ggtitle ("Other-cause death") + ylim(0, 0.6)
# g3
# g2 / g3
