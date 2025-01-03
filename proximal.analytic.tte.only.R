# library(pci2s)
devtools::load_all('/Users/george/Research_Local/pci2s_gcl/pci2s')
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
alpha  <- 1

variable.selection  <- 'automatic'
subset.names  <- list(    'all.gte.65')
 subset.name  <- 'sens1'
 subset.name  <- 'all.gte.65'
 subset.names  <- list(     'sens1', 'all.gte.65')
  # subset.names  <- list(     'sens1')
 for (subset.name in subset.names ){
        # Each analysis is separate, so set the seed for reproducibility. If
        # it's set outside the loop then the results can only be reproduced by
        # rerunning the entire script.
        set.seed(3)
        outcome.names  <-  c( 'death', 'death.cause.specific', 'death.90.day',  'death.other.cause.gt90day')
        outcome.names  <-  c( 'death', 'death.cause.specific', 'death.other.cause.gt90day',  'death.copd','death.heart',   'death.noncopd.nonheart')
          # W1  <-  'death.nonheart'
           # W1  <-  'death.noncopd'
           # W1  <-  'death.other.cause.gt90day'
           # W1  <-  'death.nonstroke'
           # W1  <-  'death.nonstroke'
        # W1s  <- c('death.other.cause.gt90day', 'death.copd', 'death.heart', 'death.noncopd.nonheart', 'death.noncopd.nonheart.nonstroke')
        W1s  <-  list(
                    'death' = 'death.other.cause.gt90day', 
                    'death.cause.specific' = 'death.other.cause.gt90day', 
                    'death.other.cause.gt90day' = 'death.other.cause.gt90day', 
                    'death.copd' = 'death.noncopd', 
                    'death.heart' = 'death.nonheart',
                    # 'death.other' = 'death.nonother',
                    'death.stroke' = 'death.nonstroke',
                     'death.noncopd.nonheart.nonstroke' = 'death.copd.heart.stroke' 
                    )
        outcome.names  <- names(W1s)
        analysis.name  <- sprintf('v3.tte.%s', subset.name)
        print('====================')
        print(analysis.name)
        print('====================')

tt.min <- 90
        ################################
        # Load data 
        ################################
        filename.in  <-  sprintf('data/A.final05.%s.RDS', subset.name)
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

        A.final %>% filter (other.cause.mortality == 'Death') %>% count(COD_TO_SITE_RECODE) %>% arrange(desc(n)) %>% print(n =Inf)
        # For the sensitivity analysis, some node positive patients are included
        A.final$tnm.n[is.na(A.final$tnm.n)]  <- 'X'
        A.final  <- A.final  %>% filter (tnm.n %in% c('0', '1', '2')) 
         print(sprintf('%.3f%% of the sublobar patients are N+', 100*sum(A.final$tnm.n != 0 & A.final$tx == 'sublobar')/ sum(A.final$tx == 'sublobar')))
        A.final$tx  <-  droplevels(A.final$tx)

        # preprocessing
        A.final <- A.final %>% mutate( race2 = ifelse (race == 'White' , 'White', 'Other'),
                                      treatment.year2 = as.character(treatment.year),
                                      treatment.year2 = (ifelse (treatment.year2 %in% c('2019', '2020'), '2019_2020', treatment.year2)),
                                      treatment.year2 = (ifelse (treatment.year2 %in% c('2010', '2011'), '2010_2011', treatment.year2)),
                                      histology2 = ifelse (grepl('Adeno', histology), 'Adenocarcinoma', 'Squamous cell'))

        # Define X
        X.factor  <-  c('sex', 'race2',  'histology2', 'treatment.year2' ) # X.factor  <-  c('sex', 'race',  'histology2') 
        X.numeric  <-  c('age', 'size') 
        Xs  <-  c(sprintf('%s_z', X.numeric),  X.factor)

           Z.count  <- c('O2accessories', 'walking_aids' ,  'wheelchairs_accessories' , 'transportation_services', 'other_supplies',   'pressure_ulcer', 'ischemic_heart_disease', 'CHF', 'PVD', 'CVD',    'MILDLD', 'DIAB_UC', 'DIAB_C',  'RD', 'mental_disorders', 'nervous_system',    'echo',  'Anticoags',  'smoking', 'o2',  'pneumonia_and_influenza','asthma', 'COPD','interstitial_lung')
        Z.count.unscaled = sprintf( '%s_pre_12months_count', Z.count )
        # Zs  <-   sprintf( '%s_pre_12months_count_bool', Z.count )
         Zs  <-   sprintf( '%s_pre_12months_count', Z.count )
        Q.count  <-  c( 'fall',  'other_injury', 'diverticular_disease', 'hernia',  'arthropathy','GU_sx')
        Qs.unscaled = sprintf( '%s_pre_12months_count', Q.count )

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
            for (i in 1:length(Zs)) cat(label_list[[gsub('_count_s', '_count_bool', Zs)[i]]], ',')
            print('Q')
        }

        #################################
        ## Table 1 
        #################################
        if (T) {
            Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")
            tblcontrol <- tableby.control(numeric.stats = c('Nmiss', 'meansd'), numeric.simplify = T, cat.simplify =T, digits = 1,total = T,test = F)
            # f  <-  sprintf( 'tx ~ %s', paste( c(X.numeric, X.factor,'histology', gsub('_count_s', '_count_bool', Zs)), collapse = "+"))
            f  <-  sprintf( 'tx ~ %s', paste( c(X.numeric, X.factor,'histology', Zs, Qs.unscaled), collapse = "+"))
            labels(A.final)  <-  label_list
            tt <- tableby(as.formula(f), data=A.final, control = tblcontrol)
            summary(tt) %>% write2html(sprintf('/Users/george/Research_Local/SEER-Medicare/tbls/table1_2_%s.htm', analysis.name), quiet=T)
        }


        #################################
        # # variable selection for W1
        ################################
        exclude.from.stage1  <-  c('size_z', 'histology2')
        X_  <-  model.matrix(as.formula(sprintf('~ %s', paste(Xs[!(Xs) %in% exclude.from.stage1], collapse = '+'))),  A.final)[,-1]
        Z_  <- A.final[,Zs]
        A_ = (A.final$tx == 'sbrt')*1.0
        design.mat  <-  as.matrix(cbind(A_, X_, Z_))
        selected.columns <- list()
        for (W1i in 1:length(W1s)) {
            outcome.name  <- W1s[[W1i]]
            A.temp  <-  A.final %>% mutate( 
                                           outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                                           outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                                           outcome.time  = if_else (outcome.time == 0, 0.5, outcome.time)/ 365+ runif(dim(A.final)[1] ,0,1)*1e-8
            ) 
            mm  <-   tune.ahazpen( Surv(A.temp$outcome.time, A.temp$outcome.bool), design.mat )
            selected.columns[[outcome.name]]  <- get.selected.columns.ahaz(mm, s = 'lambda.min', colnames(design.mat), verbose=T, min.vars = 1)
            # lambda.min.idx  <-  which.min(mm$tunem)
            # mm$fit$beta[,lambda.min.idx]
            # colnames(design.mat)
            # plot(mm)
            print(sprintf( 'Included variables in stage 1 for %s: %s', outcome.name, paste( selected.columns[[outcome.name]], collapse = ', ')))
        }

        ################################
        # Obtain estimates using the selected variables 
        ################################
        # Time-to-event outcomes 
        hazard.differences.outcomes  <-  make.odds.ratio.df ( outcome.names) 
        hazard.differences.outcomes.adj  <-  make.odds.ratio.df ( outcome.names) 
        hazard.differences.outcomes.proximal  <-  make.odds.ratio.df ( outcome.names) 
        mouts  <-  list()
        for (outcome.i in 1:length(outcome.names)){ 
            outcome.name  <-  outcome.names[outcome.i]
            W1  <- W1s[[outcome.i]]
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
	     mout  <-  two.step.loglinear.variable.selection(
				data.mat = A.final, 
				outcome.name = outcome.name,
				Xs = Xs,
				Zs = Zs,
                                W1.selected.variables = selected.columns[[W1]],
				Y.count.bool =	F,
				verbose=T,
				W1 = W1
)
            mouts[[outcome.name]]  <-  mout
            est <- mout$ESTIMATE[1]
            se  <- mout$SE[1]
            hazard.differences.outcomes.proximal[outcome.i,1:3]  <-  c(est, est - 1.96*se, est + 1.96*se)
            print(hazard.differences.outcomes[outcome.i,1:3])
            print(hazard.differences.outcomes.adj[outcome.i,1:3])
            print(hazard.differences.outcomes.proximal[outcome.i,1:3])
        } 

         saveRDS(hazard.differences.outcomes, sprintf('data/%s.hazard.differences.outcomes.raw.rds', analysis.name))
         saveRDS(hazard.differences.outcomes.adj, sprintf('data/%s.hazard.differences.outcomes.adj.rds', analysis.name))
        saveRDS(hazard.differences.outcomes.proximal, sprintf('data/%s.hazard.differences.outcomes.proximal.rds', analysis.name))
}

