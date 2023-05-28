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
filename.in  <-  sprintf('data/A.final.%s.RDS', subset.name)
A.final  <-  readRDS(filename.in)  %>% filter (time.enrolled > 0 ) %>% mutate(treatment.year = year(tx.date))

label_list  <-  readRDS('data/label.list.RDS')

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


comorbidities  <-  c('DM','DMcx', 'LiverMild', 'Pulmonary', 'PVD', 'CHF', 'MI', 'Renal', 'Stroke',  'PUD', 'Rheumatic', 'Dementia', 'LiverSevere', 'Paralysis', 't_stage_8')
negative.outcomes.oi  <-  c( 'fall',  'other_injury', 'diverticular_disease', 'hernia',  'arthropathy','GU_sx',  'hpb',  'gout', 'oral', 'optho' )
outcome.names  <-  c( 'death' )
noc.any.count.names = sprintf( '%s_any_count', negative.outcomes.oi )
noc.post.count.names = sprintf( '%s_post_count', negative.outcomes.oi )
noc.pre.count.names = sprintf( '%s_pre_count', negative.outcomes.oi )
noc.tte.names = sprintf( '%s', negative.outcomes.oi )
label_list2  <-  c( label_list,
                   death = '**Death**', 
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
                   hpb = '*HPB-related*',
                   optho = '*Ophthalmic*',
                   oral = '*Oral*'
)

################################
#   Section I: Cox models for outcomes and negative controls
################################
# Raw
outcome.names.temp  <-  c(outcome.names, noc.tte.names)
hazard.ratios.outcomes  <-  make.odds.ratio.df ( outcome.names.temp) 
outcome.i  <-  1
for (outcome.i in 1:length(outcome.names.temp)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    # print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365,
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                          treatment.year = year(tx.date)
        )
    m  <-  coxph( Surv(outcome.time, outcome.bool) ~ tx,  data = A.temp)
    # print(summary(m))
    hazard.ratios.outcomes[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'], confint(m,'txsbrt'))) 
} 
hazard.ratios.outcomes
g1.a  <-  make.OR.plot(hazard.ratios.outcomes, label_list2) + ggtitle('Raw')

# Adjusting for X
 adjust.for  <-  setdiff( c('age', 'sex', 'race', 'marital.status', 'histology','treatment.year', comorbidities) , c()) # For now, removing some of the rarer comorbiditis
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

hazard.differences.outcomes  <-  make.odds.ratio.df ( outcome.names.temp) 
outcome.i  <-  1
for (outcome.i in 1:length(outcome.names.temp)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    # print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365,
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                          treatment.year = year(tx.date)
        )
    m  <-  aalen( Surv(outcome.time, outcome.bool) ~ const(tx) ,  data = A.temp, robust = 0)
    # print(summary(m))
    hazard.differences.outcomes[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
} 
hazard.differences.outcomes
g1.a  <-  make.HD.plot(hazard.differences.outcomes, label_list2) + ggtitle('Outcomes')



# adjusting
hazard.differences.outcomes.adj  <-  make.odds.ratio.df ( outcome.names.temp) 
for (outcome.i in 1:length(outcome.names.temp)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365,
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                          treatment.year = year(tx.date)
        )
    f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + %s',  paste(sprintf('const(%s)', adjust.for), collapse="+") )
    m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
    # print(summary(m))
    hazard.differences.outcomes.adj[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
} 
hazard.differences.outcomes.adj
g1.b  <-  make.HD.plot(hazard.differences.outcomes.adj, label_list2) + ggtitle('Outcomes')
g  <- g1.a / g1.b 
ggsave(g, width=7, height=5, filename = sprintf('figs/aalen.adjX.%s.pdf', subset.name))



# Adjusting,b ut using time-varying covariates
hazard.differences.outcomes.adj  <-  make.odds.ratio.df ( outcome.names.temp) 
for (outcome.i in 1:length(outcome.names.temp)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365,
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                          treatment.year = year(tx.date)
        )
    f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + %s',  paste(sprintf('(%s)', adjust.for), collapse="+") )
    m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
    # print(summary(m))
    hazard.differences.outcomes.adj[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
} 
hazard.differences.outcomes.adj
g1.b  <-  make.HD.plot(hazard.differences.outcomes.adj, label_list2, xlim=c(-0.1, 0.7)) + ggtitle('Outcomes')
ggsave(g1.b, width=7, height=2.5, filename = sprintf('figs/aalen.adjX.time.varying.%s.pdf', subset.name))


################################
# Section III: Poisson Models for negative control outcomes
################################

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
# Section IV: Two-step proximal adjustment with Aalen's additive hazards model for outcomes and negative controls 
################################

two.step.aalen.sumW  <-  function(A.final2, Z, noc.names.temp, outcome.name, adjust.for){
        A.temp  <-  A.final2 %>% mutate( 
                              outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                              outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                              W = rowSums( ( across( all_of(noc.names.temp))))
        )
        f  <-  sprintf( 'W ~ tx + %s + %s', Z,  paste(sprintf('%s', adjust.for), collapse="+") )
        stage1  <- glm( as.formula(f) , family = poisson(link='log'), data = A.temp)
        negative_outcome_pred  <-  predict(stage1, type = 'link')
        A.temp  <-  A.final2 %>% mutate( 
                              outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                              outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                                negative_outcome_pred = negative_outcome_pred)
        f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + const(negative_outcome_pred) +  %s',  paste(sprintf('const(%s)', adjust.for), collapse="+") )
        m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
        boot.res  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')])
}

set.seed(3)
Z  <-  'optho_pre_count'
B  <-  1000
hazard.differences.outcomes.two.step  <-  make.odds.ratio.df ( outcome.names.temp) 
for (outcome.i in 1:( length(outcome.names.temp) - 1)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    noc.names.temp  <- setdiff( noc.pre.count.names, c(Z, sprintf('%s_pre_count', outcome.name)))
    print(outcome.name)
    mout  <-  two.step.aalen.sumW(A.final, Z, noc.names.temp, outcome.name, adjust.for)
    est_boot <- parallel::mclapply(1:B, function(bb){
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        boot.res  <-  two.step.aalen.sumW(A.final2, Z, noc.names.temp, outcome.name, adjust.for)
        return(boot.res)
    }, mc.cores =8)
    est_boot
    se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
   hazard.differences.outcomes.two.step[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
   print(hazard.differences.outcomes.two.step[outcome.i,1:3])
}
g3.a  <-  make.HD.plot(hazard.differences.outcomes.two.step[-nrow(hazard.differences.outcomes.two.step),], label_list2, xlims = c(-0.08, 0.25))
ggsave(g3.a, width=7, height=2.5, filename = sprintf('figs/aalen.sumW.%s.pdf', subset.name))


################################
# Section V: Two-step proximal adjustment with poisson models for negative controls 
################################

two.step.poisson.sumW  <-  function(A.final2, Z, noc.names.temp, outcome.name, adjust.for){
        A.temp  <-  A.final2 %>% mutate( 
                             outcome.count  = !!rlang::sym(outcome.name),
                              W = rowSums( ( across( all_of(noc.names.temp))))
        )
        f  <-  sprintf( 'W ~ tx + %s + %s', Z,  paste(sprintf('%s', adjust.for), collapse="+") )
        stage1  <- glm( as.formula(f) , family = poisson(link='log'), data = A.temp)
        negative_outcome_pred  <-  predict(stage1, type = 'link')
        A.temp  <-  A.final2 %>% mutate( 
                               outcome.count  = !!rlang::sym(outcome.name),
                                negative_outcome_pred = negative_outcome_pred)
        f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.enrolled))  + %s + negative_outcome_pred',  paste(sprintf('%s', adjust.for), collapse="+") )
        m  <-  glm( f  ,  data = A.temp, family = poisson(link=log))
       boot.res  <-  exp(c( coef(m)['txsbrt'] )) # to use profile confidence intervals (slower), remove the .default
}
set.seed(3)
Z  <-  'optho_pre_count'
B  <-  1000
poisson.nocs.two.step  <-  make.odds.ratio.df ( noc.any.count.names) 
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
   poisson.nocs.two.step[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
   print(poisson.nocs.two.step[outcome.i,1:3])
}
rownames(poisson.nocs.two.step)  <-  gsub('_any_count', '', rownames(poisson.nocs.two.step))
g3.a  <-  make.OR.plot(poisson.nocs.two.step[-nrow(poisson.nocs.two.step),], label_list2)
ggsave(g3.a, width=7, height=2.5, filename = sprintf('figs/poisson.sumW.%s.pdf', subset.name))


################################
# Section VI: Cox models 
################################
source("stats/nc_ph.R")
# outcomes
set.seed(3)
Z  <-  'optho_pre_count'
# Z  <-  'oral_any_count'
hazard.ratios.outcomes.two.step  <-  make.odds.ratio.df ( outcome.names.temp) 

for (outcome.i in 4:length(outcome.names.temp)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    print(outcome.name)
    noc.names.temp.no.Z  <- setdiff( noc.names, c('optho_pre_count', sprintf('%s_pre_count', outcome.name)))
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
    )
    # Make variables race and histology into dummy variables using dummyVars
    dummyVarsOut  <- caret::dummyVars('~race + histology', data = A.temp, fullRank = T)
    dummy.vars  <-  predict(dummyVarsOut, newdata = A.temp)
    adjust.for.2  <- c(adjust.for[!adjust.for %in% c('race', 'histology')], colnames(dummy.vars))
    A.temp.mat  <-  cbind(A.temp, dummy.vars) %>% 
        select(tx, all_of(adjust.for.2), outcome.bool, outcome.time,all_of(Z), all_of(noc.names.temp.no.Z)) %>%
        mutate_at( all_of( c('treatment.year', 'age')), list(~scale(.))) %>% 
        as.matrix()
    m <- cox_nc((A.temp.mat), A = "tx", X = adjust.for.2,
                Y = "outcome.time", Z = Z, W = noc.names.temp.no.Z,
                event = "outcome.bool", ncores = 8, B= 1000)
    hazard.ratios.outcomes.two.step[outcome.i,1:3]  <-  exp(c( m$estimate, m$estimate - 1.96*m$se, m$estimate + 1.96*m$se))
    print(hazard.ratios.outcomes.two.step[outcome.i,1:3])
}
# saveRDS(hazard.ratios.outcomes.two.step, sprintf('data/hazard.ratios.outcomes.two.step.%s.rds', subset.name))
g3.a  <-  make.OR.plot(hazard.ratios.outcomes.two.step[-nrow(hazard.differences.outcomes.two.step),], label_list2) + ggtitle('Adjusting for treatment year and using proximal approach ')
ggsave(g3.a, width=7, height=2.5, filename = sprintf('figs/cox.proximal.%s.pdf', subset.name))
