# library(pci2s)
# devtools::load_all('/Users/george/Research_Local/pci2s_gcl/pci2s')
devtools::load_all('/Users/george/Research_Local/pci2s/')
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

subset.names  <- list(      'all.gte.65', 'sens1')
nc_time_days = 90
nboot = 1000
for (subset.name in subset.names ){
        # Each analysis is separate, so set the seed for reproducibility. If
        # it's set outside the loop then the results can only be reproduced by
        # rerunning the entire script.
        set.seed(3)
        # outcome.names  <-  c( 'death', 'death.cause.specific', 'death.other.cause.gt90day',  'death.copd','death.heart',   'death.noncopd.nonheart')
        W1s  <-  list(
                    # 'death' = 'death.other.cause.gt90day', 
                    'death.cause.specific' = 'death.other.cause.gt90day', 
                     'death.other.cause.gt90day' = 'death.other.cause.gt90day', 
                    'death.copd' = 'death.noncopd', 
                    'death.heart' = 'death.nonheart',
                    'death.stroke' = 'death.nonstroke',
                     'death.noncopd.nonheart.nonstroke' = 'death.copd.heart.stroke' 
                    )
        outcome.names  <- names(W1s)
        # After tt.min (days), W is a valid negative outcome
        tt.min <- -1
        analysis.name  <- sprintf('v12.tte.nc_time_ss_%d_%s', nc_time_days, subset.name)
        print('====================')
        print(analysis.name)
        print('====================')

        ################################
        # Load data 
        ################################
        filename.in  <-  sprintf('data/A.final12.%s.RDS', subset.name)
        A.final  <-  readRDS(filename.in)  %>% 
            mutate(treatment.year = year(tx.date),
                   death.90.day = if_else ( ninety.day.mortality, death, as.Date(NA_character_)),
                   death.cause.specific = if_else ( cause.specific.mortality == 'Death', death, as.Date(NA_character_)),
                   # death.other.cause = if_else ( other.cause.mortality == 'Death', death, as.Date(NA_character_)),
                   death.other.cause.gt90day = if_else ( other.cause.mortality == 'Death' & tt > tt.min, death, as.Date(NA_character_)),
  		     death.copd = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE == '50130' & tt > tt.min ,  death, as.Date(NA_character_)),
		     death.heart = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE == '50060'& tt > tt.min , death, as.Date(NA_character_)),
		     death.other = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE == '50300' & tt > tt.min , death, as.Date(NA_character_)),
		     death.nonother = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50300'& tt > tt.min , death, as.Date(NA_character_)),
		     death.stroke = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE == '50080'& tt > tt.min , death, as.Date(NA_character_)),
		     death.nonstroke = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50080'& tt > tt.min , death, as.Date(NA_character_)),
  		     death.noncopd = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50130' & tt > tt.min ,  death, as.Date(NA_character_)),
		     death.nonheart = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50060'& tt > tt.min , death, as.Date(NA_character_)),
		     death.noncopd.nonheart = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50130' & COD_TO_SITE_RECODE != '50060'& tt > tt.min  , death, as.Date(NA_character_)),
		     death.noncopd.nonheart.nonstroke = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50130' & COD_TO_SITE_RECODE != '50060' & COD_TO_SITE_RECODE != '50080'& tt > tt.min  , death, as.Date(NA_character_)),
		     death.copd.heart.stroke = if_else ( other.cause.mortality == 'Death' & ( COD_TO_SITE_RECODE == '50130' | COD_TO_SITE_RECODE == '50060' | COD_TO_SITE_RECODE == '50080' ) & tt > tt.min  , death, as.Date(NA_character_)),
		    # death.noncopd = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50130'  , death, as.Date(NA_character_)),
		    # death.nonheart = if_else ( other.cause.mortality == 'Death' &  COD_TO_SITE_RECODE != '50060' & tt > 90, death, as.Date(NA_character_)),
                   pre.tx.days = pre.tx.months * 30.5,
                   post.tx.months = tt/30.5
            )
        A.final %>%   
            mutate(COD = ifelse( COD_TO_SITE_RECODE %in% cod.df$COD_TO_SITE_RECODE, COD_TO_SITE_RECODE, '50300')) %>%
            group_by(COD) %>% summarise(  n = n(), p = n()/nrow(.) ) %>% arrange(desc(n)) %>%   left_join(cod.df, by = c('COD'='COD_TO_SITE_RECODE'))  %>% mutate (np = sprintf('%d (%.1f%%)', n ,p*100 ))%>% select( COD, Name,  np) %>% print(n =10)
        A.final %>% summarise(sum(other.cause.mortality == 'Death'),  mean(other.cause.mortality == 'Death') *100) 

        table( A.final$tx, useNA="ifany")
        # For the sensitivity analysis, some node positive patients are included
        A.final$tnm.n[is.na(A.final$tnm.n)]  <- 'X'
        A.final  <- A.final  %>% filter (tnm.n %in% c('0', '1', '2')) 
        A.final  <- A.final %>% filter (tx == 'sbrt' |
                                        ( tx == 'sublobar' & tnm.n == '0' ) | 
                                        (tx == 'sublobar' & tnm.n != '0' & REGIONAL_NODES_EXAMINED_1988 != '00' ) )
        print(sprintf('%.3f%% of the sublobar patients are N+', 100*sum(A.final$tnm.n != 0 & A.final$tx == 'sublobar')/ sum(A.final$tx == 'sublobar')))
        # print(table( A.final$REGIONAL_NODES_EXAMINED_1988, A.final$tnm.n, useNA="ifany"))
        A.final$tx  <-  droplevels(A.final$tx)

        # preprocessing
        A.final <- A.final %>% mutate( race2 = ifelse (race == 'White' , 'White', 'Other'),
                                      treatment.year2 = as.character(treatment.year),
                                      treatment.year2 = (ifelse (treatment.year2 %in% c('2019', '2020'), '2019_2020', treatment.year2)),
                                      treatment.year2 = (ifelse (treatment.year2 %in% c('2010', '2011'), '2010_2011', treatment.year2)),
                                      histology2 = ifelse (grepl('Adeno', histology), 'Adenocarcinoma', 'Squamous cell'))

        # Define X
        X.factor  <-  c('sex', 'race2',  'histology2', 'treatment.year2' ) 
        X.numeric  <-  c('age', 'size') 
        Xs  <-  c(sprintf('%s_z', X.numeric),  X.factor)

        Z.count  <- c('O2accessories', 'walking_aids' ,  'wheelchairs_accessories' , 'transportation_services', 'other_supplies',   'pressure_ulcer', 'ischemic_heart_disease', 'CHF', 'PVD', 'CVD',    'MILDLD','MSLD', 'DIAB_UC', 'DIAB_C',  'RD', 'mental_disorders', 'nervous_system',    'echo',  'Anticoags',  'smoking', 'o2',  'pneumonia_and_influenza','asthma', 'COPD','interstitial_lung')
        Z.count.unscaled = sprintf( '%s_pre_12months_count', Z.count )
         Zs  <-   sprintf( '%s_pre_12months_count', Z.count )
        # Q.count  <-  c( 'fall',  'other_injury', 'diverticular_disease', 'hernia',  'arthropathy','GU_sx')
        # Qs.unscaled = sprintf( '%s_pre_12months_count', Q.count )

        A.final  <- A.final %>% mutate( 
			     time.offset = pre.tx.months, 
                                       across( all_of(c(Z.count.unscaled)), function(x) (x >0), .names = "{.col}_bool" ),
                                       # across( all_of(c(Z.count.unscaled)), function(x) quartile(x), .names = "{.col}_s" ),
                                       across( all_of(c(X.numeric)), scale_, .names = "{.col}_z" ))

        if (F) {
            print('X')
            for (i in 1:length(Xs)) cat(i, Xs[i], '\n')
            print('Z')
             # for (i in 1:length(Zs)) cat(i, label_list[[gsub('_count_s', '_count_bool', Zs)[i]]], '\n')
             # for (i in 1:length(Zs)) cat(sprintf('%s, ',label_list[ Zs[i]]))
              for (i in 1:length(Zs)) cat(sprintf('%s\n ', Zs[i]))
        }

        #################################
        ## Table 1 
        #################################
        if (T) {
            Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")
            tblcontrol <- tableby.control(numeric.stats = c('Nmiss', 'meansd'), numeric.simplify = T, cat.simplify =T, digits = 1,total = T,test = F)
             f  <-  sprintf( 'tx ~ %s', paste( c(X.numeric, X.factor,'histology', gsub('_count', '_count_bool', c(Zs))), collapse = "+"))
            # f  <-  sprintf( 'tx ~ %s', paste( c(X.numeric, X.factor, Zs, Qs.unscaled), collapse = "+"))
            labels(A.final)  <-  label_list
            tt <- tableby(as.formula(f), data=A.final, control = tblcontrol)
            # summary(tt) %>% write2html(sprintf('/Users/george/Research_Local/SEER-Medicare/tbls/table1_2_%s.bool.htm', analysis.name), quiet=T)
            summary(tt) %>% write2html(sprintf('/Users/george/Research_Local/SEER-Medicare/tbls/table1_2_%s.htm', analysis.name), quiet=T)
        }


        ################################
        # Obtain estimates using the selected variables 
        ################################
        # Time-to-event outcomes 
        hazard.differences.outcomes  <-  make.odds.ratio.df ( outcome.names) 
        hazard.differences.outcomes.adj  <-  make.odds.ratio.df ( outcome.names) 
        hazard.differences.outcomes.proximal  <-  make.odds.ratio.df ( outcome.names) 
        mouts  <-  list()
        outcome.i  <-  2
        for (outcome.i in 1:length(outcome.names)){ 
            outcome.name  <-  outcome.names[outcome.i]
            W1  <- W1s[[outcome.i]]
            print(outcome.name)
            A.temp  <-  A.final %>% mutate( 
                                           outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                                           outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                                           outcome.time  = if_else (outcome.time == 0, 0.5, outcome.time)/ 365+ runif(dim(A.final)[1] ,0,1)*1e-8
            ) 
            A.temp  <-  A.temp %>% mutate( 
                                          W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt),
                                          W1.time  = if_else (W1.time == 0, 0.5, W1.time)/365,
                                          W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
                                          Y.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                                          Y.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
            )
            A_ = (A.temp$tx == 'sbrt')*1.0
            X_  <-  model.matrix(as.formula(sprintf('~ %s', paste(Xs, collapse = '+'))),  A.temp)[,-1]
            W_ <- as.matrix( A.temp$W1.time)
            D2_  <- A.temp$W1.bool*1.0
            Z_  <- A.temp[,Zs]
            N_  <- dim(A.temp)[1]
            Y_ = A.temp$Y.time
            D_ = A.temp$Y.bool*1.0
            A.temp <- A.temp %>% mutate(
                                        cause = case_when (
                                                           nna(!!rlang::sym(outcome.name))~ 0, # primary event
                                                           nna(!!rlang::sym(W1))~ 1,
                                                           T ~ -1
                                                           )
                                        )


            A.temp %>% count(cause, cause.specific.mortality, other.cause.mortality, )
            A.temp %>% count(cause, cause.specific.mortality, nna(death.copd), nna(death.noncopd) )
            # Raw
            f  <- as.formula('Surv(outcome.time, outcome.bool) ~ const(tx)')
            m  <-  aalen( f ,  data = A.temp, robust = 0)
            hazard.differences.outcomes[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
            # Adjust for X
            f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + %s',  paste(sprintf('const(%s)', c(Xs, Zs)), collapse="+") )
            m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
            hazard.differences.outcomes.adj[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
            # Proximal
            if (outcome.name != W1) {
                print(system.time(mout  <- p2sls.cprisk.nc  (times = Y_, cause = A.temp$cause, A = A_,  X = X_, Z = Z_, nc_time = nc_time_days/365,  bootstrap = T, nboot = nboot, conf.level = 0.95) ))
                mouts[[outcome.name]]  <-  mout
                est <- mout$beta_a
                hazard.differences.outcomes.proximal[outcome.i,1:3]  <-  c(est, mout$beta_a_ci[1], mout$beta_a_ci[2])
            }else {
                hazard.differences.outcomes.proximal[outcome.i,1:3]  <-  c(NULL, NULL, NULL)

            }
            print(hazard.differences.outcomes[outcome.i,1:3])
            print(hazard.differences.outcomes.adj[outcome.i,1:3])
            print(hazard.differences.outcomes.proximal[outcome.i,1:3])
        } 

         saveRDS(hazard.differences.outcomes, sprintf('data/%s.hazard.differences.outcomes.raw.rds', analysis.name))
         saveRDS(hazard.differences.outcomes.adj, sprintf('data/%s.hazard.differences.outcomes.adj.rds', analysis.name))
        saveRDS(hazard.differences.outcomes.proximal, sprintf('data/%s.hazard.differences.outcomes.proximal.rds', analysis.name))
}

