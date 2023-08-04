library(dplyr)
source("stats/nc_ph.R")
library(patchwork)
library(lubridate)
library(foreach)
library(doParallel)
library(timereg)
library('ggtext')
library(arsenal)
library(ggplot2)
source('utilities.R')
set.seed(3)
#subset.name <- 'age.gte.80'
subset.name <- 'all.gte.65'
filename.in  <-  sprintf('data/A.final4.%s.RDS', subset.name)
A.final  <-  readRDS(filename.in)  %>% filter (time.enrolled > 0 ) %>% 
    mutate(treatment.year = year(tx.date),
            death.90.day = if_else ( ninety.day.mortality, death, as.Date(NA_character_)),
            death.cause.specific = if_else ( cause.specific.mortality == 'Death', death, as.Date(NA_character_)),
            death.cause.specific = if_else ( cause.specific.mortality == 'Death', death, as.Date(NA_character_)),
            death.other.cause = if_else ( other.cause.mortality == 'Death', death, as.Date(NA_character_)),
    
    )

label_list  <-  readRDS('data/label.list.RDS')

comorbidities  <-  c( 'DM','DMcx', 'LiverMild', 'Pulmonary', 'PVD', 'CHF', 'MI', 'Renal',   'PUD', 'Rheumatic', 'Dementia', 'LiverSevere', 'Paralysis', 't_stage_8')
W.prefix  <- 'any_'


# [ ]  Use other outcomes, especially non-cause mortality

# Not using: non_diabetes_endocrine
X.factor  <-  c('t_stage_8','sex', 'race', 'marital.status', 'histology' ) 
X.numeric  <-  c('age', 'treatment.year', sprintf('%s_pre_count', c(
                    # DMEs
                      'hospital_beds_and_supplies', 'wheelchairs_accessories', 'walking_aids', 'O2accessories', 'other_supplies', 'diabetic_footwear', 'transportation_services', 
                    # Diagnoses
                       'smoking', 'o2', 'other_bacterial_diseases', 'pneumonia_and_influenza', 'pressure_ulcer', 'ischemic_heart_disease', 'CHF', 'PVD', 'CVD', 'dementia', 'COPD', 'PUD', 'MILDLD', 'DIAB_UC', 'DIAB_C', 'PARA', 'RD', 'cancer_nonlung', 'MSLD', 'METS',  'mental_disorders', 'nervous_system', 'other_heart_disease', 'veins_lymphatics_other_circulatory', 'rheum',
                    # Drugs
                    'Insulin', 'Anticoags')
                        )) 

# Get number of non-zero values for each X.numeric variable
A.final %>% select(all_of(X.numeric)) %>% summarise_all(list(~sum(. > 0, na.rm = T))) %>% t() %>% as.data.frame() %>% arrange(desc(V1))



adjust.for  <-  setdiff( c(X.numeric, X.factor) , c()) # For now, removing some of the rarer comorbiditis
# Use dplyr to scale all variables in the numeric.X list
scale_  <-  function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T)

# Make each X.numeric variable into a z-score, but name them with suffix _z
A.final  <- A.final %>% mutate( 
    across( all_of(X.numeric), scale_ , .names = "{.col}_z" )
)
adjust.for.scaled  <-  setdiff( c(sprintf('%s_z', X.numeric), X.factor) , c()) # For now, removing some of the rarer comorbiditis

################################
# File structure 
################################
# Each row is a patient, identified by PATIENT_ID. The date of treatment is tx.date, and that is considered time zero. tt is the time from treatment to either death or end of trial.  The various comorbidities
# as obtained by the Quan et al. scores (e.g. DM, or Diabetes Mellitus) are
# boolean, and TRUE if there is a corresponding code BEFORE treatment. The negative
# outcomes, such as fall, have multiple corresponding variables, as follows:
# fall: date of first code after treatment
# fall_pre: date of first code before treatment
# fall_pre_count: count of codes after treatment
# fall_pre_date_count: count of dates with codes after treatment
# fall_any: date of first code 
# fall_post: date of first code after treatment
# fall_post_count: count of codes after treatment
# fall_post_date_count: count of dates with codes after treatment
# fall_any_count: count of codes at any time
# fall_any_date_count: count of dates with codes at any time


negative.outcomes.oi  <-  c( 'fall',  'other_injury', 'diverticular_disease', 'hernia',  'arthropathy','GU_sx',  'pancreatic',   'oral', 'optho' )
outcome.names  <-  c( 'death', 'death.cause.specific', 'death.90.day', 'death.other.cause' )
noc.any.count.names = sprintf( '%s_any_count', negative.outcomes.oi )
noc.post.count.names = sprintf( '%s_post_count', negative.outcomes.oi )
noc.pre.count.names = sprintf( '%s_pre_count', negative.outcomes.oi )
noc.tte.names = sprintf( '%s', negative.outcomes.oi )
label_list2  <-  c( label_list,
                   death = '**Overall mortality**', 
                   death.cause.specific = '**Cancer-specific mortality**', 
                   death.other.cause = '**Other mortality**', 
                   death.90.day = '**90-day mortality**', 
                   fall = '*Fall*',
                   other_injury = '*Injury*',
                   GU_sx = '*GU-related*',
                   arthropathy = '*Arthropathy*',
                   cholelithiasis = '*Cholelithiasis-related*',
                   gout = '*Gout*',
                   obstruction = '*Intestinal obstruction*',
                   hernia = '*Abdominal hernia*',
                   diverticular_disease = '*Diverticular disease*',
                   hemorrhoids = '*Hemorrhoids*',
                   pancreatic = '*Pancreatic*',
                   optho = '*Ophthalmic*',
                   oral = '*Oral*'
)

################################
#   Section I: Cox models for outcomes and negative controls
################################
# Raw
print('Cox model with and without no adjustment')
outcome.names.temp  <-  c(outcome.names, noc.tte.names)
hazard.ratios.outcomes  <-  make.odds.ratio.df ( outcome.names.temp) 
outcome.i  <-  1
for (outcome.i in 1:length(outcome.names.temp)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    # print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365,
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        )
    m  <-  coxph( Surv(outcome.time, outcome.bool) ~ tx,  data = A.temp)
    # print(summary(m))
    hazard.ratios.outcomes[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'], confint(m,'txsbrt'))) 
} 
hazard.ratios.outcomes
g1.a  <-  make.OR.plot(hazard.ratios.outcomes, label_list2) + ggtitle('Raw')

outcome.i  <-  1
hazard.ratios.outcomes.adj  <-  make.odds.ratio.df ( outcome.names.temp) 
for (outcome.i in 1:length(outcome.names.temp)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    # print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365,
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        )
    f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ tx + %s + treatment.year',  paste(sprintf('%s', adjust.for), collapse="+") )
    m  <-  coxph( as.formula(f) ,  data = A.temp)
     # print(summary(m))
    hazard.ratios.outcomes.adj[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'], confint(m,'txsbrt'))) 
} 
hazard.ratios.outcomes.adj
g2.a  <-  make.OR.plot(hazard.ratios.outcomes.adj, label_list2) + ggtitle('Adjusting for X')
g  <- g1.a / g2.a #+ plot_annotation( title = 'Cox model, W as TTE')
 ggsave(g, width=7, height=5, filename = sprintf('figs/cox.adjX.%s.pdf', subset.name))

################################
# Section II: Aalen's additive hazards model for outcomes and negative controls  
################################
print('Aalen model with and without no adjustment')

hazard.differences.outcomes  <-  make.odds.ratio.df ( outcome.names.temp) 
outcome.i  <-  1
for (outcome.i in 1:length(outcome.names.temp)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    # print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365,
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        )
    m  <-  aalen( Surv(outcome.time, outcome.bool) ~ const(tx) ,  data = A.temp, robust = 0)
    # print(summary(m))
    hazard.differences.outcomes[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
} 
hazard.differences.outcomes
g1.a  <-  make.HD.plot(hazard.differences.outcomes, label_list2)



# adjusting
hazard.differences.outcomes.adj  <-  make.odds.ratio.df ( outcome.names.temp) 
for (outcome.i in 1:length(outcome.names.temp)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365,
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        )
    f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + %s',  paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
    m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
    # print(summary(m))
    hazard.differences.outcomes.adj[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
} 
hazard.differences.outcomes.adj
g1.b  <-  make.HD.plot(hazard.differences.outcomes.adj, label_list2)
g  <- g1.a / g1.b 
ggsave(g, width=7, height=5, filename = sprintf('figs/aalen.adjX.%s.pdf', subset.name))



# # Adjusting,b ut using time-varying covariates
# hazard.differences.outcomes.adj  <-  make.odds.ratio.df ( outcome.names.temp) 
# for (outcome.i in 1:length(outcome.names.temp)){ 
#     outcome.name  <-  outcome.names.temp[outcome.i]
#     print(outcome.name)
#     A.temp  <-  A.final %>% mutate( 
#                           outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365,
#                           outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
#         )
#     f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + %s',  paste(sprintf('(%s)', adjust.for.scaled), collapse="+") )
#     m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
#      print(summary(m))
#     hazard.differences.outcomes.adj[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
# } 
# hazard.differences.outcomes.adj
# g1.b  <-  make.HD.plot(hazard.differences.outcomes.adj, label_list2) 
# ggsave(g1.b, width=7, height=2.5, filename = sprintf('figs/aalen.adjX.time.varying.%s.pdf', subset.name))

# # Diagnosing the time-varying covariates

# outcome.name  <-  'gout'
# A.temp  <-  A.final %>% mutate( 
#                       outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365,
#                       outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
#     )
# m  <-  aalen( as.formula('Surv(outcome.time, outcome.bool) ~ const(tx) + age' ) ,  data = A.temp, robust = 0, silent = F)
# print(summary(m))
# m  <-  aalen( as.formula('Surv(outcome.time, outcome.bool) ~ const(tx) + const(age)' ) ,  data = A.temp, robust = 0, silent = F)
# print(summary(m))

# data(sTRACE)
# sTRACE %>% filter (no == '3763')
# out<-aalen(Surv(time,status==9)~age+sex+diabetes+chf+vf, sTRACE,max.time=7,n.sim=100, silent =F, id = sTRACE$no)
# sum(duplicated(out$no))
# summary(out)

################################
# Section IVa: Two-step proximal adjustment with Aalen's additive hazards model for outcomes and negative controls , using other mortality as an egative control
################################
Zs  <-  c('O2accessories_pre_count_z', 'walking_aids_pre_count_z' , 'hospital_beds_and_supplies_pre_count_z' , 'wheelchairs_accessories_pre_count_z' , 'transportation_services_pre_count_z', 'other_supplies_pre_count_z', 'diabetic_footwear_pre_count_z' )
adjust.for  <-  setdiff( c(X.numeric, X.factor) , Zs) # For now, removing some of the rarer comorbiditis
adjust.for.scaled  <-  setdiff( c(sprintf('%s_z', X.numeric), X.factor) , Zs) # For now, removing some of the rarer comorbiditis
#TODO: Figureo ut why some people died but are not of cause specific or other cause
print("Two-step proximal adjustment with Aalen's additive hazards model for outcomes and negative controls, using other mortality as an egative control")
# A.final2 %>% select( death, death.cause.specific, death.other.cause, tx.date, tt, diverticular_disease_pre_count) %>% print ( n = 100)

 f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ (tx) +%s+ %s',  paste(sprintf('(%s)', Zs), collapse="+"), paste(sprintf('(%s)', adjust.for.scaled), collapse="+") )
 summary(coxph(as.formula(f), data = A.temp))


two.step.aalen.othermortality  <-  function(A.final2, Zs,  outcome.name, adjust.for.scaled, W= NA){
    A.temp  <-  A.final2 %>% mutate( 
                                          outcome.time  = if_else (nna(!!rlang::sym(W)), as.numeric( !!rlang::sym(W) - tx.date, units = "days" ), tt)/365,
                                          outcome.bool = ifelse( nna(!!rlang::sym(W)), T, F),
                                    )
        # f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + const(treatment.year_z)+  const(Z) ' )
        # f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ (tx) +  (arthropathy_pre_count)+ (treatment.year_z) ' )
        # f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ tx + ' )
         f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) +%s+ %s',  paste(sprintf('const(%s)', Zs), collapse="+"), paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
        # f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ (tx) +%s+ %s',  paste(sprintf('(%s)', Zs), collapse="+"), paste(sprintf('(%s)', adjust.for.scaled), collapse="+") )
          # summary(coxph(as.formula(f), data = A.temp))
        m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0, silent = 0 )
         # summary(m)
        mm  <-  model.matrix(as.formula(f), A.temp)[,-1] 
        coefs  <-  as.matrix(coef( m)[,'Coef.'])
        negative_outcome_pred  <-  mm %*% coefs
        A.temp  <-  A.final2 %>% mutate( 
                              outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                              outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                                negative_outcome_pred = scale(negative_outcome_pred[,1]))
        f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + const(negative_outcome_pred) +  %s',  paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
        m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0, silent = 0)
        #TODO: If you set silent = 0, you'll see some of these are singular, are those biasing the results?
        boot.res  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')])
        return(boot.res)
}
set.seed(3)
B  <-  1000
hazard.differences.outcomes.two.step  <-  make.odds.ratio.df ( outcome.names.temp) 
for (outcome.i in 1:( length(outcome.names.temp))){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    mout  <-  two.step.aalen.othermortality(A.final, Zs,  outcome.name, adjust.for.scaled, 'death.other.cause')
    est_boot <- parallel::mclapply(1:B, function(bb){
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        mout  <-  two.step.aalen.othermortality(A.final2, Zs,  outcome.name, adjust.for.scaled, 'death.other.cause')
        return(mout)
    }, mc.cores =8)
    est_boot
    se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
   hazard.differences.outcomes.two.step[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
   print(hazard.differences.outcomes.two.step[outcome.i,1:3])
}

g3.a  <-  make.HD.plot(hazard.differences.outcomes.two.step, label_list2, xlims = c(-0.08, 0.25))
g  <- (g1.a + ggtitle('No adjustment')) / (g1.b + ggtitle('Adjusting for X'))  / (g3.a  + ggtitle( 'Two stage, DME as Zs, other mortality as W'))
ggsave(g, width=7, height=7.5, filename = sprintf('figs/aalen.adjX.othermortality.%s.pdf', subset.name))
ggsave(g3.a, width=7, height=2.5, filename = sprintf('figs/aalen.othermortality.%s.pdf', subset.name))



################################
# Section IVaa: Two-step proximal adjustment with Aalen's additive hazards model for outcomes and negative controls , using other mortality as an egative control and sumW
################################
Zs  <-  c('O2accessories_pre_count_z', 'walking_aids_pre_count_z' , 'hospital_beds_and_supplies_pre_count_z' , 'wheelchairs_accessories_pre_count_z' , 'transportation_services_pre_count_z', 'other_supplies_pre_count_z', 'diabetic_footwear_pre_count_z' )
adjust.for  <-  setdiff( c(X.numeric, X.factor) , Zs) # For now, removing some of the rarer comorbiditis
adjust.for.scaled  <-  setdiff( c(sprintf('%s_z', X.numeric), X.factor) , Zs) # For now, removing some of the rarer comorbiditis
print("Two-step proximal adjustment with Aalen's additive hazards model for outcomes and negative controls, using other mortality as an egative control")

 # f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ (tx) +%s+ %s',  paste(sprintf('(%s)', Zs), collapse="+"), paste(sprintf('(%s)', adjust.for.scaled), collapse="+") )
 # summary(coxph(as.formula(f), data = A.temp))


two.step.aalen.othermortality  <-  function(A.final2, Zs,  outcome.name,noc.names.temp, adjust.for.scaled, W= NA){
    A.temp  <-  A.final2 %>% mutate( 
                                          outcome.time  = if_else (nna(!!rlang::sym(W)), as.numeric( !!rlang::sym(W) - tx.date, units = "days" ), tt)/365,
                                          outcome.bool = ifelse( nna(!!rlang::sym(W)), T, F),
                                    )
         f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) +%s+ %s',  paste(sprintf('const(%s)', Zs), collapse="+"), paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
        m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0, silent = 0 )
        mm  <-  model.matrix(as.formula(f), A.temp)[,-1] 
        coefs  <-  as.matrix(coef( m)[,'Coef.'])
        negative_outcome_pred  <-  mm %*% coefs
        Ws  <-  noc.names.temp
        A.temp  <-  A.final2 %>% mutate( 
                             outcome.count  = !!rlang::sym(outcome.name),
                              W = rowSums( ( across( all_of(Ws))))
        )
        f  <-  sprintf( 'W ~ tx + %s + %s', paste(sprintf('const(%s)', Zs), collapse="+"),  paste(sprintf('%s', adjust.for), collapse="+") )
        stage1  <- glm( as.formula(f) , family = poisson(link='log'), data = A.temp)
        negative_outcome_pred2  <-  predict(stage1, type = 'link')
        A.temp  <-  A.final2 %>% mutate( 
                              outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                              outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                                negative_outcome_pred = scale(negative_outcome_pred[,1]))
        f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + const(negative_outcome_pred) + const(negative_outcome_pred2) +  %s',  paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
        m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0, silent = 0)
        boot.res  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')])
        return(boot.res)
}

set.seed(3)
B  <-  500
hazard.differences.outcomes.two.step  <-  make.odds.ratio.df ( outcome.names.temp) 
for (outcome.i in 1:( length(outcome.names.temp))){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    noc.names.temp  <- setdiff( noc.any.count.names, c(outcome.name))
    mout  <-  two.step.aalen.othermortality(A.final, Zs,  outcome.name, noc.names.temp, adjust.for.scaled, 'death.other.cause')
    est_boot <- parallel::mclapply(1:B, function(bb){
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        mout  <-  two.step.aalen.othermortality(A.final2, Zs,  outcome.name, noc.names.temp, adjust.for.scaled, 'death.other.cause')
        return(mout)
    }, mc.cores =8)
    est_boot
    se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
   hazard.differences.outcomes.two.step[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
   print(hazard.differences.outcomes.two.step[outcome.i,1:3])
}

g3.a  <-  make.HD.plot(hazard.differences.outcomes.two.step, label_list2, xlims = c(-0.08, 0.25))
g  <- (g1.a + ggtitle('No adjustment')) / (g1.b + ggtitle('Adjusting for X'))  / (g3.a  + ggtitle( 'Two stage, DME as Zs, other mortality as W'))
ggsave(g, width=7, height=7.5, filename = sprintf('figs/aalen.adjX.othermortality.sumW.%s.pdf', subset.name))











################################
# Section IVb: Two-step proximal adjustment with Aalen's additive hazards model for outcomes and negative controls 
################################
print("Two-step proximal adjustment with Aalen's additive hazards model for outcomes and negative controls")

two.step.aalen.sumW  <-  function(A.final2, Z, noc.names.temp, outcome.name, adjust.for.scaled, W= NA){
        A.temp  <-  A.final2 %>% mutate( 
                              outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                              outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        )
        if (is.na(W)) {
            A.temp <- A.temp %>% mutate( W = rowSums( ( across( all_of(noc.names.temp)))) )
        }else{
            A.temp <- A.temp %>% mutate( W = !!rlang::sym(W) )
        }
        f  <-  sprintf( 'W ~ tx + %s + %s', Z,  paste(sprintf('%s', adjust.for.scaled), collapse="+") )
        stage1  <- glm( as.formula(f) , family = poisson(link='log'), data = A.temp)
        summary(stage1)
        negative_outcome_pred  <-  predict(stage1, type = 'link')
        A.temp  <-  A.final2 %>% mutate( 
                              outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                              outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                                negative_outcome_pred = negative_outcome_pred)
        f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + const(negative_outcome_pred) +  %s',  paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
        m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
        boot.res  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')])
}

set.seed(3)
Z  <-  'optho_pre_count'
B  <-  1000
hazard.differences.outcomes.two.step  <-  make.odds.ratio.df ( outcome.names.temp) 
for (outcome.i in 1:( length(outcome.names.temp) - 1)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    noc.names.temp  <- setdiff( noc.any.count.names, c(Z, sprintf('%s_any_count', outcome.name)))
    print(outcome.name)
    mout  <-  two.step.aalen.sumW(A.final, Z, noc.names.temp, outcome.name, adjust.for.scaled)
    est_boot <- parallel::mclapply(1:B, function(bb){
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        boot.res  <-  two.step.aalen.sumW(A.final2, Z, noc.names.temp, outcome.name, adjust.for.scaled)
        return(boot.res)
    }, mc.cores =8)
    est_boot
    se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
   hazard.differences.outcomes.two.step[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
   print(hazard.differences.outcomes.two.step[outcome.i,1:3])
}
g3.a  <-  make.HD.plot(hazard.differences.outcomes.two.step[-nrow(hazard.differences.outcomes.two.step),], label_list2, xlims = c(-0.08, 0.25))
ggsave(g3.a, width=7, height=2.5, filename = sprintf('figs/aalen.sumanyW.%s.pdf', subset.name))



################################
# Section IVc: Two-step proximal adjustment with Aalen's additive hazards model for outcomes and negative controls , multiple Zs
################################
print("Two-step proximal adjustment with Aalen's additive hazards model for outcomes and negative controls, multiple Zs")

two.step.aalen.sumW  <-  function(A.final2, Z, noc.names.temp, outcome.name, adjust.for.scaled, W= NA){
        A.temp  <-  A.final2 %>% mutate( 
                              outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                              outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        )
        if (is.na(W)) {
            A.temp <- A.temp %>% mutate( W = rowSums( ( across( all_of(noc.names.temp)))) )
        }else{
            A.temp <- A.temp %>% mutate( W = !!rlang::sym(W) )
        }
        f  <-  sprintf( 'W ~ tx + %s', paste(sprintf('%s', c( Z, adjust.for.scaled)), collapse="+") )
        stage1  <- glm( as.formula(f) , family = poisson(link='log'), data = A.temp)
        negative_outcome_pred  <-  predict(stage1, type = 'link')
        A.temp  <-  A.final2 %>% mutate( 
                              outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                              outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                                negative_outcome_pred = negative_outcome_pred)
        f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + const(negative_outcome_pred) +  %s',  paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
        m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
        boot.res  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')])
}

set.seed(3)
Z  <-  c( 'optho_pre_count', 'diverticular_disease_pre_count', 'arthropathy_pre_count' )
B  <-  1000
outcome.names.temp.2  <-  setdiff( outcome.names.temp, Z)
hazard.differences.outcomes.two.step  <-  make.odds.ratio.df ( outcome.names.temp.2) 
for (outcome.i in 1:( length(outcome.names.temp.2))){ 
    outcome.name  <-  outcome.names.temp.2[outcome.i]
    noc.names.temp  <- setdiff( noc.pre.count.names, c(Z, sprintf('%s_any_count', outcome.name)))
    print(outcome.name)
    mout  <-  two.step.aalen.sumW(A.final, Z, noc.names.temp, outcome.name, adjust.for.scaled)
    est_boot <- parallel::mclapply(1:B, function(bb){
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        boot.res  <-  two.step.aalen.sumW(A.final2, Z, noc.names.temp, outcome.name, adjust.for.scaled)
        return(boot.res)
    }, mc.cores =8)
    est_boot
    se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
   hazard.differences.outcomes.two.step[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
   print(hazard.differences.outcomes.two.step[outcome.i,1:3])
}
g3.a  <-  make.HD.plot(hazard.differences.outcomes.two.step[-nrow(hazard.differences.outcomes.two.step),], label_list2, xlims = c(-0.08, 0.25))
ggsave(g3.a, width=7, height=2.5, filename = sprintf('figs/aalen.sumanyW.Zoptho.divertic.arthropathy.%s.pdf', subset.name))



# Choose the W to be the one with most correlation with A and Z, as it is thus
# closest related to U, since any correlation with A and Z is necessarily
# mediated by U


# rbind(cor (A.final[[Z]], A.final[noc.pre.count.names] ), cor (A.final[['tx']] == 'sbrt', A.final[noc.pre.count.names] )) %>% t
# set.seed(3)
# Z  <-  'optho_pre_count'
# W  <-  'other_injury_pre_count'
# B  <-  1000
# outcome.names.temp  <-  c(outcome.names, noc.tte.names)
# hazard.differences.singleW.outcomes.two.step  <-  make.odds.ratio.df ( outcome.names.temp) 
# for (outcome.i in 1:( length(outcome.names.temp) - 1)){ 
#     outcome.name  <-  outcome.names.temp[outcome.i]
#     noc.names.temp  <- setdiff( noc.pre.count.names, c(Z, sprintf('%s_pre_count', outcome.name)))
#     print(outcome.name)
#     mout  <-  two.step.aalen.sumW(A.final, Z, noc.names.temp, outcome.name, adjust.for.scaled, W=W)
#     est_boot <- parallel::mclapply(1:B, function(bb){
#         A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
#         boot.res  <-  two.step.aalen.sumW(A.final2, Z, noc.names.temp, outcome.name, adjust.for.scaled, W=W)
#         return(boot.res)
#     }, mc.cores =8)
#     est_boot
#     se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
#    hazard.differences.singleW.outcomes.two.step[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
#    print(hazard.differences.singleW.outcomes.two.step[outcome.i,1:3])
# }
# g3.a  <-  make.HD.plot(hazard.differences.singleW.outcomes.two.step[-nrow(hazard.differences.singleW.outcomes.two.step),], label_list2, xlims = c(-0.08, 0.25))
# ggsave(g3.a, width=7, height=2.5, filename = sprintf('figs/aalen.W%s.%s.pdf', W,subset.name))

# # Try to project a random variable onto the covariates, it will be high variance. That is what's happening.

# fofo.temp  <-  setdiff(noc.pre.count.names,c(Z,W))
# A.temp <- A.final %>% mutate( W.all = rowSums( ( across( all_of(noc.names.temp)))) )
# rbind(cor (A.temp[[Z]], A.temp$W.all ), cor (A.temp[['tx']] == 'sbrt', A.temp$W.all )) %>% t
# fofo.temp  <-  setdiff(noc.pre.count.names,c(Z,W, 'gout_pre_count', 'pancreatic_pre_count'))
# A.temp <- A.final %>% mutate( W.all = rowSums( ( across( all_of(fofo.temp)))) )
# rbind(cor (A.temp[[Z]], A.temp$W.all ), cor (A.temp[['tx']] == 'sbrt', A.temp$W.all )) %>% t

# fofo.temp  <-  setdiff(noc.pre.count.names,c(Z,W))
# f.temp  <-  as.formula(sprintf('(tx=="sbrt") ~ %s', paste(sprintf('%s', fofo.temp), collapse="+") ))
# A.temp$predictZ  <-  predict(lm(f.temp, data = A.final))
# rbind(cor (A.temp[[Z]], A.temp$predictZ ), cor (A.temp[['tx']] == 'sbrt', A.temp$predictZ )) %>% t
# sqrt(0.06506)


## Instead of summing all the W, use individual Ws
## Make an empty data frame with the following columns: outcome, W_i, A.Z.cor, Z.outcome.cor, W.A.cor, se
#se.df  <- data.frame(outcome = character(), W_i = numeric(), A.Z.cor = numeric(), A.W.cor = numeric(), se = numeric(), stringsAsFactors = F)
## add a new row to se.df
#set.seed(3)
#Z  <-  'optho_pre_count'
#B  <-  500
#for (outcome.i in 1:( length(outcome.names.temp) - 1)){ 
#    outcome.name  <-  outcome.names.temp[outcome.i]
#    hazard.differences.outcomes.two.step.Wi  <-  make.odds.ratio.df ( outcome.names.temp) 
#    #Get first part of string Z, before the character _
#    W_is  <- setdiff( negative.outcomes.oi, c(strsplit(Z, '_')[[1]][1], outcome.name))
#    for (noc.j  in 1:length(W_is) ) {
#        W_i  <- sprintf('%s_pre_count', W_is[noc.j])
#        noc.names.temp  <- setdiff( noc.pre.count.names, c(Z, sprintf('%s_pre_count', outcome.name)))
#        print(sprintf('%s - %s', outcome.name, W_i))
#        mout  <-  two.step.aalen.sumW(A.final, Z, noc.names.temp, outcome.name, adjust.for.scaled, W= W_i)
#        est_boot <- parallel::mclapply(1:B, function(bb){
#            A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
#            boot.res  <-  two.step.aalen.sumW(A.final2, Z, noc.names.temp, outcome.name, adjust.for.scaled, W= W_i)
#            return(boot.res)
#        }, mc.cores =8)
#        est_boot
#        se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
#       hazard.differences.outcomes.two.step.Wi[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
#       print(hazard.differences.outcomes.two.step.Wi[outcome.i,1:3])
#       A.Z.cor  <-  cor(A.final[[Z]], A.final[['tx']] == 'sbrt')
#       A.W.cor  <-  cor(A.final[[W_i]], A.final[['tx']] == 'sbrt')
#       se.df[nrow(se.df) + 1,]  <-  c(outcome.name, W_i, A.Z.cor, A.W.cor, se)
#    }
#    g3.a  <-  make.HD.plot(hazard.differences.outcomes.two.step.Wi[-nrow(hazard.differences.outcomes.two.step.Wi),], label_list2, xlims = c(-0.08, 0.25))
#    ggsave(g3.a, width=7, height=2.5, filename = sprintf('figs/aalen.W%s.%s.pdf', W_i, subset.name))
#}
#se.df %>% filter(outcome =='fall') %>% select( A.W.cor, se)
#cor(se.df$A.W.cor %>% as.numeric, se.df$se %>% as.numeric, use = 'complete.obs')


# # only with X = treatment.year
# set.seed(3)
# Z  <-  'optho_pre_count'
# B  <-  1000
# adjust.for.2  <-  c('treatment.year', 'histology')
# hazard.differences.outcomes.two.step  <-  make.odds.ratio.df ( outcome.names.temp) 
# for (outcome.i in 1:( length(outcome.names.temp) - 1)){ 
#     outcome.name  <-  outcome.names.temp[outcome.i]
#     noc.names.temp  <- setdiff( noc.pre.count.names, c(Z, sprintf('%s_pre_count', outcome.name)))
#     print(outcome.name)
#     mout  <-  two.step.aalen.sumW(A.final, Z, noc.names.temp, outcome.name, adjust.for.2)
#     est_boot <- parallel::mclapply(1:B, function(bb){
#         A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
#         boot.res  <-  two.step.aalen.sumW(A.final2, Z, noc.names.temp, outcome.name, adjust.for.2)
#         return(boot.res)
#     }, mc.cores =8)
#     est_boot
#     se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
#    hazard.differences.outcomes.two.step[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
#    print(hazard.differences.outcomes.two.step[outcome.i,1:3])
# }
# g3.a  <-  make.HD.plot(hazard.differences.outcomes.two.step[-nrow(hazard.differences.outcomes.two.step),], label_list2, xlims = c(-0.08, 0.25))
# ggsave(g3.a, width=7, height=2.5, filename = sprintf('figs/aalen.sumW.Xonlytxyear.%s.pdf', subset.name))

################################
# Section IV: Poisson Models for negative control outcomes
################################

print('Poisson models without and without adjustment')
# Raw
odds.ratios.nocs  <-  make.odds.ratio.df ( noc.any.count.names) 
for (outcome.i in 1:length(noc.any.count.names)){ 
    outcome.name  <-  noc.any.count.names[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( outcome.count  = !!rlang::sym(outcome.name))
    f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.enrolled))' )
    m  <-  glm( f  ,  data = A.temp, family = poisson(link=log))
    odds.ratios.nocs[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'] , confint(m)['txsbrt',]))
}
odds.ratios.nocs 
rownames(odds.ratios.nocs)  <-  gsub('_any_count', '', rownames(odds.ratios.nocs))
g2.a  <-  make.OR.plot(odds.ratios.nocs, label_list2) + ggtitle('Raw')

# Adjusting for X
odds.ratios.nocs.adj  <-  make.odds.ratio.df ( noc.any.count.names) 
for (outcome.i in 1:length(noc.any.count.names)){ 
    outcome.name  <-  noc.any.count.names[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( outcome.count  = !!rlang::sym(outcome.name))
    f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.enrolled))  + %s',  paste(sprintf('%s', adjust.for), collapse="+") )
    m  <-  glm( f  ,  data = A.temp, family = poisson(link=log))
    # print(summary(m))
    odds.ratios.nocs.adj[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'] , confint(m)['txsbrt',]))
}
odds.ratios.nocs.adj 
rownames(odds.ratios.nocs.adj)  <-  gsub('_any.count', '', rownames(odds.ratios.nocs.adj))
g2.b  <-  make.OR.plot(odds.ratios.nocs.adj, label_list2) + ggtitle('Adjusting for X')
g2  <- g2.a / g2.b 
ggsave(g2, width=7, height=5, filename = sprintf('figs/poisson.adjX.%s.pdf', subset.name))


################################
# Section Va: Two-step proximal adjustment with poisson models for negative controls,
# using other mortality outcomes as negative controls
################################

two.step.poisson.sumW  <-  function(A.final2, Zs,  outcome.name, adjust.for.scaled){
    W  <- 'death.other.cause'
    A.temp  <-  A.final2 %>% mutate( 
                                          outcome.time  = if_else (nna(!!rlang::sym(W)), as.numeric( !!rlang::sym(W) - tx.date, units = "days" ), tt)/365,
                                          outcome.bool = ifelse( nna(!!rlang::sym(W)), T, F),
                                    )
         f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) +%s+ %s',  paste(sprintf('const(%s)', Zs), collapse="+"), paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
        m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0, silent = 0 )
        mm  <-  model.matrix(as.formula(f), A.temp)[,-1] 
        coefs  <-  as.matrix(coef( m)[,'Coef.'])
        negative_outcome_pred  <-  mm %*% coefs
        A.temp  <-  A.final2 %>% mutate( 
                               outcome.count  = !!rlang::sym(outcome.name),
                                negative_outcome_pred = negative_outcome_pred)
        f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.enrolled))  + %s + negative_outcome_pred',  paste(sprintf('%s', adjust.for.scaled), collapse="+") )
        m  <-  glm( f  ,  data = A.temp, family = poisson(link=log))
        summary(m)
       boot.res  <-  exp(c( coef(m)['txsbrt'] ))
}

set.seed(3)
B  <-  1000
poisson.nocs.two.step  <-  make.odds.ratio.df ( noc.any.count.names) 
ses  <- rep(NA, length(noc.any.count.names))
for (outcome.i in 1:( length(noc.any.count.names) - 1)){ 
    outcome.name  <-  noc.any.count.names[outcome.i]
    print(outcome.name)
    mout  <-  two.step.poisson.sumW(A.final, Zs,  outcome.name, adjust.for.scaled)
    est_boot <- parallel::mclapply(1:B, function(bb){
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        boot.res  <-  two.step.poisson.sumW(A.final2, Zs,  outcome.name, adjust.for.scaled)
        return(boot.res)
    }, mc.cores =8)
    est_boot
    se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
    ses[outcome.i]  <-  se
   poisson.nocs.two.step[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
   print(poisson.nocs.two.step[outcome.i,1:3])
}
rownames(poisson.nocs.two.step)  <-  gsub('_any_count', '', rownames(poisson.nocs.two.step))
g3.a  <-  make.OR.plot(poisson.nocs.two.step[-nrow(poisson.nocs.two.step),], label_list2)
ggsave(g3.a, width=7, height=2.5, filename = sprintf('figs/poisson.sumW.%s.noXinstage1.pdf', subset.name))
names(ses)  <- noc.any.count.names
################################
# Section Vb: Two-step proximal adjustment with poisson models for negative controls,
# using other mortality outcomes as negative controls and also an extra W
################################

two.step.poisson.sumW  <-  function(A.final2, Zs,  outcome.name, noc.names.temp, adjust.for.scaled){
    W  <- 'death.other.cause'
    A.temp  <-  A.final2 %>% mutate( 
                                          outcome.time  = if_else (nna(!!rlang::sym(W)), as.numeric( !!rlang::sym(W) - tx.date, units = "days" ), tt)/365,
                                          outcome.bool = ifelse( nna(!!rlang::sym(W)), T, F),
                                    )
         f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) +%s+ %s',  paste(sprintf('const(%s)', Zs), collapse="+"), paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
        m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0, silent = 0 )
        mm  <-  model.matrix(as.formula(f), A.temp)[,-1] 
        coefs  <-  as.matrix(coef( m)[,'Coef.'])
        negative_outcome_pred  <-  mm %*% coefs
        Ws  <-  noc.names.temp
        A.temp  <-  A.final2 %>% mutate( 
                             outcome.count  = !!rlang::sym(outcome.name),
                              W = rowSums( ( across( all_of(Ws))))
        )
        f  <-  sprintf( 'W ~ tx + %s + %s', paste(sprintf('const(%s)', Zs), collapse="+"),  paste(sprintf('%s', adjust.for), collapse="+") )
        stage1  <- glm( as.formula(f) , family = poisson(link='log'), data = A.temp)
        negative_outcome_pred2  <-  predict(stage1, type = 'link')
        A.temp  <-  A.final2 %>% mutate( 
                               outcome.count  = !!rlang::sym(outcome.name),
                                negative_outcome_pred = negative_outcome_pred,
                                negative_outcome_pred2 = negative_outcome_pred2,
        )
        f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.enrolled))  + %s+ negative_outcome_pred2 + negative_outcome_pred',  paste(sprintf('%s', adjust.for.scaled), collapse="+") )
        m  <-  glm( f  ,  data = A.temp, family = poisson(link=log))
        summary(m)
       boot.res  <-  exp(c( coef(m)['txsbrt'] ))
}

set.seed(3)
B  <-  500
poisson.nocs.two.step  <-  make.odds.ratio.df ( noc.any.count.names) 
ses  <- rep(NA, length(noc.any.count.names))
for (outcome.i in 1:( length(noc.any.count.names) - 1)){ 
    outcome.name  <-  noc.any.count.names[outcome.i]
    print(outcome.name)
    noc.names.temp  <- setdiff( noc.any.count.names, c(outcome.name))
    mout  <-  two.step.poisson.sumW(A.final, Zs,  outcome.name, noc.names.temp,adjust.for.scaled)
    est_boot <- parallel::mclapply(1:B, function(bb){
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        boot.res  <-  two.step.poisson.sumW(A.final2, Zs,  outcome.name, noc.names.temp,adjust.for.scaled)
        return(boot.res)
    }, mc.cores =8)
    est_boot
    se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
    ses[outcome.i]  <-  se
   poisson.nocs.two.step[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
   print(poisson.nocs.two.step[outcome.i,1:3])
}
rownames(poisson.nocs.two.step)  <-  gsub('_any_count', '', rownames(poisson.nocs.two.step))
g3.a  <-  make.OR.plot(poisson.nocs.two.step[-nrow(poisson.nocs.two.step),], label_list2)
ggsave(g3.a, width=7, height=2.5, filename = sprintf('figs/poisson.sumW.%s.noXinstage1.pdf', subset.name))
names(ses)  <- noc.any.count.names


################################
# Section V: Two-step proximal adjustment with poisson models for negative controls 
################################

two.step.poisson.sumW  <-  function(A.final2, Z, noc.names.temp, outcome.name, adjust.for){
        A.temp  <-  A.final2 %>% mutate( 
                             outcome.count  = !!rlang::sym(outcome.name),
                              W = rowSums( ( across( all_of(noc.names.temp))))
        )
        f  <-  sprintf( 'W ~ tx + %s + %s', Z,  paste(sprintf('%s', adjust.for), collapse="+") )
        # f  <-  sprintf( 'W ~ tx + %s ', Z )
        stage1  <- glm( as.formula(f) , family = poisson(link='log'), data = A.temp)
        if  (coef(summary(stage1))[Z,4] > 0.05) 
            warning('Z is not related to W; likely poor choice of Z, will lead to high variance.')
        negative_outcome_pred  <-  predict(stage1, type = 'link')
        A.temp  <-  A.final2 %>% mutate( 
                               outcome.count  = !!rlang::sym(outcome.name),
                                negative_outcome_pred = negative_outcome_pred)
        f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.enrolled))  + %s + negative_outcome_pred',  paste(sprintf('%s', adjust.for), collapse="+") )
        m  <-  glm( f  ,  data = A.temp, family = poisson(link=log))
        summary(m)
       boot.res  <-  exp(c( coef(m)['txsbrt'] ))
}

set.seed(3)
Z  <-  'optho_pre_count'
B  <-  1000
poisson.nocs.two.step  <-  make.odds.ratio.df ( noc.any.count.names) 
ses  <- rep(NA, length(noc.any.count.names))
for (outcome.i in 1:( length(noc.any.count.names) - 1)){ 
    outcome.name  <-  noc.any.count.names[outcome.i]
    noc.names.temp  <- setdiff( noc.pre.count.names, c(Z, sprintf('%s_pre_count', outcome.name)))
    print(outcome.name)
    mout  <-  two.step.poisson.sumW(A.final, Z, noc.names.temp, outcome.name, adjust.for)
    est_boot <- parallel::mclapply(1:B, function(bb){
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        boot.res  <-  two.step.poisson.sumW(A.final2, Z, noc.names.temp, outcome.name, adjust.for)
        return(boot.res)
    }, mc.cores =8)
    est_boot
    se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
    ses[outcome.i]  <-  se
   poisson.nocs.two.step[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
   print(poisson.nocs.two.step[outcome.i,1:3])
}
rownames(poisson.nocs.two.step)  <-  gsub('_any_count', '', rownames(poisson.nocs.two.step))
g3.a  <-  make.OR.plot(poisson.nocs.two.step[-nrow(poisson.nocs.two.step),], label_list2)
ggsave(g3.a, width=7, height=2.5, filename = sprintf('figs/poisson.sumW.%s.noXinstage1.pdf', subset.name))
names(ses)  <- noc.any.count.names

# # correlation of Z with each of the Ws, but keep the variable names
# correlations  <- lapply( noc.any.count.names, function(x) cor(A.final[[x]], A.final[[Z]])) %>% unlist
# names(correlations)  <-  noc.any.count.names
# sort(correlations)
# sort(ses)
# cor(ses,correlations, use = 'pairwise.complete.obs')



# correlations  <- lapply( noc.any.count.names, function(x) cor(A.final[[x]], A.final[['optho_pre_count']])) %>% unlist
# names(correlations)  <-  noc.any.count.names
# sort(correlations)



# # correlation of Z with each of the Ws, but keep the variable names
# correlations  <- lapply( noc.any.count.names, function(x) cor(A.final[[x]], A.final[['tx']] == 'sbrt')) %>% unlist
# names(correlations)  <-  noc.any.count.names
# W  <- 'gout_any_count'
# f  <-  sprintf( "tx =='sbrt' ~ %s + %s",  W,paste(sprintf('%s', adjust.for), collapse="+") )
# summary(lm(f, data = A.final))
# f  <-  sprintf( "optho_pre_count ~ %s + %s", W, paste(sprintf('%s', adjust.for), collapse="+") )
# summary(lm(f, data = A.final))
# f  <-  sprintf( "oral_pre_count ~ %s + %s", W, paste(sprintf('%s', adjust.for), collapse="+") )
# summary(lm(f, data = A.final))


# f  <-  sprintf( "tx =='sbrt' ~ %s + %s + treatment.year",  'gout_post_count',paste(sprintf('%s', adjust.for), collapse="+") )
# summary(lm(f, data = A.final))

# outcome.name = 'death'
# A.temp  <-  A.final %>% mutate( 
#                       outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365,
#                       outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
#     )
# summary(coxph( Surv(outcome.time, outcome.bool) ~ fall_any_count,  data = A.temp))



# # I pick my Z and W based on the U most relevant to Y. 
# # My problem is when an N is not related to Z or W. 
# # A can totally cause gout.




# set.seed(3)
# Z  <-  'oral_pre_count'
# B  <-  1000
# noc.any.count.names.2  <- setdiff(noc.any.count.names, Z)
# poisson.nocs.two.step.2  <-  make.odds.ratio.df ( noc.any.count.names.2) 
# ses2  <- rep(NA, length(noc.any.count.names.2))
# for (outcome.i in 1:( length(noc.any.count.names.2))){ 
#     outcome.name  <-  noc.any.count.names.2[outcome.i]
#     noc.names.temp  <- setdiff( noc.pre.count.names, c(Z, sprintf('%s_pre_count', outcome.name)))
#     print(outcome.name)
#     mout  <-  two.step.poisson.sumW(A.final, Z, noc.names.temp, outcome.name, adjust.for)
#      for (i in 1:100) {
#      A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
#          boot.res  <-  two.step.poisson.sumW(A.final2, Z, noc.names.temp, outcome.name, adjust.for)
#          print(boot.res)
#          if (boot.res > 10)
#              break
#      }
#      est_boot <- parallel::mclapply(1:B, function(bb){
#          A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
#          boot.res  <-  two.step.poisson.sumW(A.final2, Z, noc.names.temp, outcome.name, adjust.for)
#          return(boot.res)
#      }, mc.cores =8)
#     est_boot
#     se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
#     ses2[outcome.i]  <-  se
#    poisson.nocs.two.step.2[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
#    print(poisson.nocs.two.step.2[outcome.i,1:3])
# }
# rownames(poisson.nocs.two.step.2)  <-  gsub('_any_count', '', rownames(poisson.nocs.two.step.2))
# g3.a  <-  make.OR.plot(poisson.nocs.two.step.2[-nrow(poisson.nocs.two.step.2),], label_list2)
# ggsave(g3.a, width=7, height=2.5, filename = sprintf('figs/poisson.Z%s.sumW.%s.pdf', Z, subset.name))
# names(ses2)  <- noc.any.count.names.2


# # correlation of Z with each of the Ws, but keep the variable names
# correlations  <- lapply( noc.any.count.names, function(x) cor(A.final[[x]], A.final[[Z]])) %>% unlist
# names(correlations)  <-  noc.any.count.names
# sort(correlations)
# sort(ses2)
# cor(ses,correlations, use = 'pairwise.complete.obs')




################################
# Section VI: Cox models 
################################
source("stats/nc_ph.R")
# outcomes
set.seed(3)
Z  <-  'optho_pre_count_z'
# Z  <-  'oral_any_count'
hazard.ratios.outcomes.two.step  <-  make.odds.ratio.df ( outcome.names.temp) 
# for (outcome.i in 1:length(outcome.names.temp)){ 
 for (outcome.i in 1:1){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    print(outcome.name)
    noc.names.temp.no.Z  <- setdiff( noc.pre.count.names, c('optho_pre_count', sprintf('%s_pre_count', outcome.name)))
    A.temp  <-  A.final %>% mutate( 
                                   outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                                   outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                                   tx = ifelse(tx == 'sbrt', 1, 0),
                                   sex = ifelse( sex == 'Male', 1, 0),
                                   marital.status = ifelse( marital.status == 'Married', 1, 0),
                                   t_stage_8 = case_when( t_stage_8 == 'T1a' ~ 1,
                                                         t_stage_8 == 'T1b' ~ 2,
                                                         t_stage_8 == 'T1c' ~ 3,
                                                         TRUE ~ 0),
    ) %>% mutate(  optho_pre_count_z= scale(optho_pre_count))
    # Make variables race and histology into dummy variables using dummyVars
    dummyVarsOut  <- caret::dummyVars('~race + histology', data = A.temp, fullRank = T)
    dummy.vars  <-  predict(dummyVarsOut, newdata = A.temp)
    adjust.for.2  <- c(adjust.for.scaled[!adjust.for.scaled %in% c('race', 'histology')], colnames(dummy.vars))
    # adjust.for.2  <- c('age_z', 'treatment.year_z', 'O2accessories_pre_count_z', 'transportation_services_pre_count_z', 'CHF_pre_count_z', 'CVD_pre_count_z', 'COPD_pre_count_z', 'PUD_pre_count_z', 'RD_pre_count_z', 'METS_pre_count_z','t_stage_8','sex', 'histologyAdenocarcinoma, NOS','histologySquamous cell carcinoma' )
      # adjust.for.2  <- c('age_z', 'treatment.year_z', 'CVD_pre_count_z', 'COPD_pre_count_z', 'PUD_pre_count_z', 'RD_pre_count_z', 'METS_pre_count_z','t_stage_8','sex', 'histologyAdenocarcinoma, NOS','histologySquamous cell carcinoma' )
      # adjust.for.2  <- c('treatment.year_z')
    A.temp.mat  <-  cbind(A.temp, dummy.vars) %>% 
        select(tx, all_of(adjust.for.2), outcome.bool, outcome.time,all_of(Z), all_of(noc.names.temp.no.Z)) %>%
        # mutate_at( all_of( c('treatment.year', 'age')), list(~scale(.))) %>% 
        as.matrix()
    m <- cox_nc(A.temp.mat, A = "tx", X = adjust.for.2,
                Y = "outcome.time", Z = Z, W = noc.names.temp.no.Z,
                event = "outcome.bool", ncores = 1, B= 1)
    ### start temp
    ### end temp
    hazard.ratios.outcomes.two.step[outcome.i,1:3]  <-  exp(c( m$estimate, m$estimate - 1.96*m$se, m$estimate + 1.96*m$se))
    print(hazard.ratios.outcomes.two.step[outcome.i,1:3])
}
# saveRDS(hazard.ratios.outcomes.two.step, sprintf('data/hazard.ratios.outcomes.two.step.%s.rds', subset.name))
g3.a  <-  make.OR.plot(hazard.ratios.outcomes.two.step[-nrow(hazard.ratios.outcomes.two.step),], label_list2) + ggtitle('Adjusting for treatment year and using proximal approach ')
ggsave(g3.a, width=7, height=2.5, filename = sprintf('figs/cox.proximal.%s.pdf', subset.name))


