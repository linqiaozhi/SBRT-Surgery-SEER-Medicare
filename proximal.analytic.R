library(pci2s)
library(glmnet)
library(ahaz)
library(arsenal)
library(dplyr)
library(tidyr)
library(MASS)
library(patchwork)
library(lubridate)
library(timereg)
library('ggtext')
library(ggplot2)
source('utilities.R')
source('two.step.variable.selection.R')
alpha  <- 1

variable.selection  <- 'automatic'
subset.names  <- list(  'sens1', 'all.gte.65')
lambda.s  <- 'lambda.min'

for (subset.name in subset.names ){
    # Cache the fits for when testing  multiple lambdas. Not done in this code.
    stage1.fits  <- list()
    stage2.fits  <- list()
        # Each analysis is separate, so set the seed for reproducibility. If
        # it's set outside the loop then the results can only be reproduced by
        # rerunning the entire script.
        set.seed(3)
        analysis.name  <- sprintf('v16.%s.%s', subset.name, lambda.s)
        print('====================')
        print(analysis.name)
        print('====================')

        ################################
        # Load data 
        ################################
        filename.in  <-  sprintf('data/A.final16.%s.RDS', subset.name)
        A.final  <-  readRDS(filename.in)  %>% 
            mutate(treatment.year = year(tx.date),
                   death.90.day = if_else ( ninety.day.mortality, death, as.Date(NA_character_)),
                   death.cause.specific = if_else ( cause.specific.mortality == 'Death', death, as.Date(NA_character_)),
                   death.other.cause = if_else ( other.cause.mortality == 'Death', death, as.Date(NA_character_)),
                   death.other.cause.gt90day = if_else ( other.cause.mortality == 'Death' & tt > 90, death, as.Date(NA_character_)),
                   pre.tx.days = pre.tx.months * 30.5,
            )

        A.final$tnm.n[is.na(A.final$tnm.n)]  <- 'X'
         print(sprintf('%.3f%% of the sublobar patients are N+', 100*sum(A.final$tnm.n != 0 & A.final$tx == 'sublobar')/ sum(A.final$tx == 'sublobar')))
        A.final$tx  <-  droplevels(A.final$tx)

        # preprocessing
        A.final <- A.final %>% mutate( race2 = ifelse (race == 'White' , 'White', 'Other'),
                                      treatment.year2 = as.character(treatment.year),
                                      treatment.year2 = (ifelse (treatment.year2 %in% c('2019', '2020'), '2019-2020', treatment.year2)),
                                      histology2 = ifelse (grepl('Adeno', histology), 'Adenocarcinoma', 'Squamous cell'))

        # Define X
        X.factor  <-  c('sex', 'race',  'histology2', 'treatment.year2' ) 
        X.numeric  <-  c('age', 'size') 
        Xs  <-  c(sprintf('%s_z', X.numeric),  X.factor)

        Z.count  <- c('O2accessories', 'walking_aids' , 'hospital_beds_and_supplies' , 'wheelchairs_accessories' , 'transportation_services', 'other_supplies', 'diabetic_footwear' ,  'pressure_ulcer', 'ischemic_heart_disease', 'CHF', 'PVD', 'CVD', 'dementia',   'MILDLD','MSLD', 'DIAB_UC', 'DIAB_C',  'RD', 'mental_disorders', 'nervous_system',   'dialysis', 'echo', 'Insulin', 'Anticoags',  'smoking', 'o2',  'pneumonia_and_influenza','asthma', 'COPD','interstitial_lung')
        Z.count.unscaled = sprintf( '%s_pre_month_count', Z.count )
        Zs  <-   sprintf( '%s_pre_month_count_s', Z.count )

        W.count  <-  c( 'fall',  'other_injury', 'diverticular_disease', 'hernia',  'arthropathy','GU_sx')
        Ws = sprintf( '%s_pre_month_count', W.count )

        outcome.names  <-  c( 'death', 'death.cause.specific', 'death.90.day',  'death.other.cause.gt90day')

        A.final  <- A.final %>% mutate( time.offset = pre.tx.months, 
                                       across( all_of(c(Z.count.unscaled, Ws)), function(x) (x >0), .names = "{.col}_unbinned" ),
                                       across( all_of(c(Z.count.unscaled)), function(x) quartile(x/time.offset), .names = "{.col}_s" ),
                                       across( all_of(c(X.numeric)), scale_, .names = "{.col}_z" ))

        if (F) {
            print('X')
            for (i in 1:length(Xs)) cat(i, Xs[i], '\n')
            print('Z')
            for (i in 1:length(Zs)) cat(i, Zs[i], '\n')
            print('W')
            for (i in 1:length(Ws)) cat(i, Ws[i], '\n')
            cat(i+1,'other.cause.mortality\n')
        }
        #################################
        ## Table 1 
        #################################
        if (T) {
            Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")
            tblcontrol <- tableby.control(numeric.stats = c('Nmiss', 'meansd'), numeric.simplify = T, cat.simplify =T, digits = 1,total = T,test = F)
            f  <-  sprintf( 'tx ~ %s', paste( c(X.numeric, X.factor,'histology', gsub('_count_s', '_count_unbinned', Zs), gsub('_count', '_count_unbinned', c(Ws))), collapse = "+"))
            labels(A.final)  <-  label_list
            tt <- tableby(as.formula(f), data=A.final, control = tblcontrol)
            summary(tt) %>% write2html(sprintf('/Users/george/Research_Local/SEER-Medicare/tbls/table1_2_%s.htm', analysis.name), quiet=T)
        }

        A_ = (A.final$tx == 'sbrt')*1.0
        X_  <-  model.matrix(as.formula(sprintf('~ %s', paste(Xs, collapse = '+'))),  A.final)[,-1]
        Z_  <- A.final[,Zs]
        W1  <-  'death.other.cause.gt90day'
            print('Automatic variable selection: stage 1 models')
            ################################
            # Variable selection: Stage 1. For each W,  first fit a lasso to determine
            # which X and Z will be included, which will make up the Xw matrix. Then, refit
            # the model without penalization to construct W_hat.
            ################################
            nfolds  <- 5
            # Prepare the lists and matrices.
            Xw =  replicate (length(Ws)+1,   list())
            names(Xw)  <- c(W1, Ws)
            W_hat =  matrix(NA, nrow = dim(A.final)[1], ncol = length(W1) + length(Ws))
            colnames(W_hat)  <- c(W1, Ws)

            # First, for W1, using the penalized additive hazard model implemented in ahaz
            A.temp  <-  A.final %>% mutate( 
                                           W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt),
                                           # Scale so units are in years and a small jitter to the time to avoid ties
                                           W1.time  = if_else (W1.time == 0, 0.5, W1.time)/365+runif(dim(A.final)[1] ,0,1)*1e-8,
                                           W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
            ) 

            design.mat  <-  as.matrix(cbind(A_, X_, Z_))
            if (!W1 %in% names(stage1.fits)) {
                stage1.fits[[W1]]  <-  tune.ahazpen( Surv(A.temp$W1.time, A.temp$W1.bool), design.mat )
                # pdf(sprintf('figs/%s.%s.lasso.pdf', analysis.name, W1))
                # plot(stage1.fits[[W1]])
                # dev.off()
            }else{
                print('Using cached fit')
            }
            # print(W1)
            selected.columns  <- get.selected.columns.ahaz(stage1.fits[[W1]], s = lambda.s, colnames(design.mat), verbose=F, min.vars = 10, denom =5)
            Xw[[1]]  <-  as.matrix(cbind( A_, cbind(X_, Z_ )[,selected.columns]))
            ## Refit model without penalization
            W_hat[,W1]  <-  lin_ah( time = A.temp$W1.time, 
                                   event = A.temp$W1.bool, 
                                   covariates = Xw[[1]]
                                   )$PREDICT

            # Next, for the rest of the Ws (the counts)
            for (i in 1:length(Ws)) {
                ## Select with a poisson L1 model
                W_i  <- A.final[,Ws[i]]
                # print(Ws[i])
                # use a glmnet to regress W_i onto A, X, Zin a lasso
                if (!W_i %in% names(stage1.fits)) {
                    stage1.fits[[Ws[i]]]  <- cv.glmnet(as.matrix(cbind(A_, X_, Z_)), as.matrix(W_i), offset= log(A.final$time.offset), family = 'poisson', alpha = 1, nfolds = nfolds)
                    # pdf(sprintf('figs/%s.%s.lasso.pdf', analysis.name, Ws[i]))
                    # plot(stage1.fits[[Ws[i]]])
                    # dev.off()
                }
                selected.columns  <- get.selected.columns(stage1.fits[[Ws[i]]],s=lambda.s,  verbose=F, min.vars = 1)
                Xw[[Ws[i]]]  <-  as.matrix(cbind(1, A_,cbind( X_, Z_ )[,selected.columns]))
                ## Refit model without penalization
                W_hat[,Ws[i]]  <- negbin_fit(y = as.matrix(W_i), 
                                             x = Xw[[Ws[i]]], 
                                             offset = log(A.final$time.offset))$PREDICT
            }

            ################################
            # Variable selection: Stage 2. For each outcome, fit a lasso to determine which
            # of the Ws and Xs to include.
            ################################
            print('Automatic variable selection: stage 2 models')
            # Next, find the variables to include in each stage 2 model. 
            # Next, for the count Ys
            selected.Ws  <- list()
            selected.Xs  <- list()
            for (outcome.i in 1:length(Ws)){ 
                outcome.name  <-  Ws[outcome.i]
                # print(outcome.name)
                A.temp  <-  A.final %>% mutate( outcome.count  = !!rlang::sym(outcome.name))
                noc.names.temp  <- setdiff( Ws, c(outcome.name))
                if (!W_i %in% names(stage2.fits)) {
                    stage2.fits[[Ws[outcome.i]]]  <-  cv.glmnet(cbind (A_, X_, W_hat[,c('death.other.cause.gt90day', noc.names.temp)] ), as.matrix(A.temp$outcome.count), offset= log(A.temp$time.offset), family = 'poisson', alpha = 1, nfolds = nfolds)
                }
                selected.columns  <- get.selected.columns(stage2.fits[[Ws[outcome.i]]], s=lambda.s, verbose=F, min.vars = 1)
                selected.Ws[[outcome.name]] <-  Ws[check.if.present(Ws,selected.columns)]
                selected.Xs[[outcome.name]] <-  Xs[check.if.present(Xs,selected.columns)]
            }

            for (outcome.i in 1:length(outcome.names)){ 
                outcome.name  <-  outcome.names[outcome.i]
                # print(outcome.name)
                A.temp  <- A.temp %>% mutate(
                                             Y.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                                             Y.time  = if_else (Y.time == 0, 0.5, Y.time)/365 +runif(dim(A.temp)[1] ,0,1)*1e-8,
                                             Y.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                )
                if (outcome.name =='death.other.cause.gt90day' ) {
                    W.hat.temp  <- W_hat[,c(noc.names.temp)]
                }else{
                    W.hat.temp  <- W_hat[,c('death.other.cause.gt90day', noc.names.temp)]
                }
                design.mat  <-  as.matrix(cbind(A_, X_, W.hat.temp))
                if (!W_i %in% names(stage2.fits)) {
                   stage2.fits[[outcome.name]] <-  tune.ahazpen( Surv(A.temp$Y.time, A.temp$Y.bool)  , design.mat )
                }
                #TODO: that shouldn't be i...
                # print(outcome.name)
                selected.columns  <- get.selected.columns.ahaz(stage2.fits[[outcome.name]], s=lambda.s, colnames(design.mat), verbose=F, min.vars = 1)
                selected.Ws[[outcome.name]] <-  Ws[check.if.present(Ws,selected.columns)]
                selected.Xs[[outcome.name]] <-  Xs[check.if.present(Xs,selected.columns)]
            }

        # Print out the variables selected at each stage
        sink(sprintf('data/variable.selection.%s.txt', analysis.name))
        # First, print the Xs and Zs for each W
        print('Stage 1')
        for (W_i in c(W1, Ws)) {
            cat( sprintf('\n\nW: %s, with the following Xs and Zs:', W_i))
            # past the colnames
            cat(paste(colnames(Xw[[W_i]]), collapse = '\n'))
        }
        print('====================')
        print('Stage 2')
        for (outcome.i in 1:length(outcome.names)){ 
            outcome.name  <-  outcome.names[outcome.i]
            cat( sprintf('\n\nOutcome: %s, with the following Ws and Xs:', outcome.name))
            for (W_i in selected.Ws[[outcome.name]]) {
                print( sprintf('-->W: %s', W_i))
            }
            for (X_i in selected.Xs[[outcome.name]]) {
                print( sprintf('-->X: %s', X_i))
            }
        }
        for (outcome.i in 1:length(Ws)){ 
            outcome.name  <-  Ws[outcome.i]
            cat( sprintf('\n\nOutcome: %s, with the following Ws and Xs:', outcome.name))
            for (W_i in selected.Ws[[outcome.name]]) {
                print( sprintf('-->W: %s', W_i))
            }
            for (X_i in selected.Xs[[outcome.name]]) {
                print( sprintf('-->X: %s', X_i))
            }
        }
        sink()



        ################################
        # Obtain estimates using the selected variables 
        ################################
        print('Estimate the hazard differences for time-to-event outcomes')

        # Time-to-event outcomes 
        hazard.differences.outcomes  <-  make.odds.ratio.df ( outcome.names) 
        hazard.differences.outcomes.adj  <-  make.odds.ratio.df ( outcome.names) 
        hazard.differences.outcomes.proximal  <-  make.odds.ratio.df ( outcome.names) 
        mouts  <-  list()
        # for (outcome.i in 2:2){ 
        for (outcome.i in 1:length(outcome.names)){ 
            outcome.name  <-  outcome.names[outcome.i]
            print(outcome.name)
            A.temp  <-  A.final %>% mutate( 
                                           outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                                           outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                                           outcome.time  = if_else (outcome.time == 0, 0.5, outcome.time)/ 365+ runif(dim(A.final)[1] ,0,1)*1e-8
            ) 
            # Raw
            f  <- as.formula('Surv(outcome.time, outcome.bool) ~ const(tx)')
            m  <-  aalen( f ,  data = A.temp, robust = 0)
            hazard.differences.outcomes[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
            # Adjust for X
            f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + %s',  paste(sprintf('const(%s)', c(Xs, Zs)), collapse="+") )
            m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
            hazard.differences.outcomes.adj[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
            # Proximal
            mout  <-  two.step.variable.selection(
                                                  data.mat = A.final, 
                                                  outcome.name = outcome.name,
                                                  X.selected.Y  = selected.Xs[[outcome.name]],
                                                  W.selected.Y  = selected.Ws[[outcome.name]],
                                                  Xw = Xw,
                                                  Y.count.bool =F, 
                                                  verbose=T,
                                                  skip.W1 = outcome.name != 'death.cause.specific'
            )
            mouts[[outcome.name]]  <-  mout
            est <- mout$ESTIMATE[1]
            se  <- mout$SE[1]
            hazard.differences.outcomes.proximal[outcome.i,1:3]  <-  c(est, est - 1.96*se, est + 1.96*se)
            print(hazard.differences.outcomes.proximal[outcome.i,1:3])
        } 


        # print('Estimate the risk ratios for count outcomes')
         # Counts
         odds.ratios.nocs  <-  make.odds.ratio.df ( Ws) 
         odds.ratios.nocs.adj  <-  make.odds.ratio.df ( Ws) 
         odds.ratios.nocs.proximal  <-  make.odds.ratio.df ( Ws) 
         for (outcome.i in 1:length(Ws)){ 
          # for (outcome.i in 1:0){ 
             outcome.name  <-  Ws[outcome.i]
             print(outcome.name)
             A.temp  <-  A.final %>% mutate( outcome.count  = !!rlang::sym(outcome.name))
             # Raw
             f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.offset))' )
             mm  <- model.matrix(as.formula(f), data=A.temp)
             m  <-  glm.nb(f, data = A.temp, control = glm.control(maxit=100), init.theta = 0.1)
             odds.ratios.nocs[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'] , confint(m)['txsbrt',]))
             # Adjust for X
             f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.offset))  + %s',  paste(sprintf('%s', c(Xs, Zs)), collapse="+") )
             mm  <- model.matrix(as.formula(f), data=A.temp)
             nb_fit.out  <- negbin_fit(A.temp$outcome.count, mm, offset = log(A.temp$time.offset))
             est  <-  nb_fit.out$ESTIMATE[3]
             se  <- nb_fit.out$SE[3]
             odds.ratios.nocs.adj[outcome.i,1:3]  <-  exp(c( est, est-1.96*se, est+1.96*se))
             # Proximal
             mout  <-  two.step.variable.selection(
                                                   data.mat = A.final, 
                                                   outcome.name = outcome.name,
                                                   X.selected.Y  = selected.Xs[[outcome.name]],
                                                   W.selected.Y  = selected.Ws[[outcome.name]],
                                                   Xw = Xw,
                                                   Y.count.bool =T, 
                                                   verbose=F)
             est <- mout$ESTIMATE['A']
             se  <- mout$SE['A']
             odds.ratios.nocs.proximal[outcome.i,1:3]  <-  exp(c( est, est- 1.96*se, est + 1.96*se ))
             # print(odds.ratios.nocs[outcome.i,1:3])
             # print(odds.ratios.nocs.adj[outcome.i,1:3])
             # print(odds.ratios.nocs.proximal[outcome.i,1:3])
             cat('==========\n\n\n')
         }


        #################################
        ## Individual plots
        #################################
         saveRDS(hazard.differences.outcomes, sprintf('data/%s.hazard.differences.outcomes.raw.rds', analysis.name))
         saveRDS(hazard.differences.outcomes.adj, sprintf('data/%s.hazard.differences.outcomes.adj.rds', analysis.name))
        saveRDS(hazard.differences.outcomes.proximal, sprintf('data/%s.hazard.differences.outcomes.proximal.rds', analysis.name))
         saveRDS(odds.ratios.nocs, sprintf('data/%s.odds.ratios.nocs.raw.rds', analysis.name))
         saveRDS(odds.ratios.nocs.adj, sprintf('data/%s.odds.ratios.nocs.adj.rds', analysis.name))
         saveRDS(odds.ratios.nocs.proximal, sprintf('data/%s.odds.ratios.nocs.proximal.rds', analysis.name))
}

sum(nna(A.final$death.cause.specific[ which(A.final$tnm.n != '0')]))/sum(A.final$tnm.n != '0')
sum(nna(A.final$death.cause.specific[ which(A.final$tnm.n == '0')]))/sum(A.final$tnm.n == '0')
