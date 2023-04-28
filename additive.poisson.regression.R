library(dplyr)
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
A.final  <-  readRDS(filename.in)  %>% filter (time.enrolled > 0 )
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
# fall_any: date of first code 
# fall_post: date of first code after treatment
# fall_post_count: count of codes after treatment
# fall_post_date_count: count of dates with codes after treatment
# fall_any_count: count of codes at any time
# fall_any_date_count: count of dates with codes at any time


################################
#   Not adjusting for X
################################
comorbidities  <-  c('DM','DMcx', 'LiverMild', 'Pulmonary', 'PVD', 'CHF', 'MI', 'Renal', 'Stroke',  'PUD', 'Rheumatic', 'Dementia', 'LiverSevere', 'Paralysis', 't_stage_8')
negative.outcomes.oi  <-  c( 'fall',  'other_injury', 'diverticular_disease', 'hernia',  'arthropathy','GU_sx',  'hpb',  'gout', 'oral', 'optho' )
outcome.names  <-  c( 'death' )
noc.names = sprintf( '%s_any_count', negative.outcomes.oi )
label_list2  <-  c( label_list,
                   death = '**Death**', 
                   fall_any_count = '*Fall*',
                   other_injury_any_count = '*Injury*',
                   GU_sx_any_count = '*GU-related*',
                   arthropathy_any_count = '*Arthropathy*',
                   cholelithiasis_any_count = '*Cholelithiasis-related*',
                   gout_any_count = '*Gout*',
                   obstruction_any_count = '*Intestinal obstruction*',
                   hernia_any_count = '*Abdominal hernia*',
                   diverticular_disease_any_count = '*Diverticular disease*',
                   hemorrhoids_any_count = '*Hemorrhoids*',
                   hpb_any_count = '*HPB-related*',
                   optho_any_count = '*Ophthalmic*',
                   oral_any_count = '*Oral*'
)

# The outcomes, which right now are only death, will be modeled using Aalen's additive hazards model.
# Eventually, there will be more outcomes, hence the loop structure for only one variable.
hazard.ratios.outcomes  <-  make.odds.ratio.df ( outcome.names[1]) 
outcome.i  <-  1
for (outcome.i in 1:length(outcome.names)){ 
    outcome.name  <-  outcome.names[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365,
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                          treatment.year = year(tx.date)
        )
    m  <-  aalen( Surv(outcome.time, outcome.bool) ~ 1+ const(tx) + (treatment.year) ,  data = A.temp, robust = 0)
    print(summary(m))
    hazard.ratios.outcomes[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
} 
hazard.ratios.outcomes
g1.a  <-  make.HD.plot(hazard.ratios.outcomes, label_list2) + ggtitle('Outcomes')

# In prior iterations, all the negative outcomes were also included as
# time-to-event models. However, there was substantial loss of information, as
# it is useful to know if someone had 10 codes for hemorrhoids vs. just the
# date of the first episode. There is also an immortal time bias. SBRT is a
# newer technology, so most SBRTs are done in more recent years. There is thus
# less follow up and less time for SBRT patients to have events, on average,
# than resection patients.  For these reasons, negative outcomes are modeled as
# rates, using Poisson models. Specifically, the number of events is modelled,
# and the log of time enrolled in the trial is used as an offset (as here, https://stats.stackexchange.com/a/11183/103007)

odds.ratios.nocs  <-  make.odds.ratio.df ( noc.names) 
for (outcome.i in 1:length(noc.names)){ 
    outcome.name  <-  noc.names[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( outcome.count  = !!rlang::sym(outcome.name))
    summary( A.temp$time.enrolled, useNA="ifany")
    m  <-  glm( outcome.count ~ 1+ tx + offset(log(time.enrolled))   ,  data = A.temp, family = poisson(link=log))
    print(summary(m))
    odds.ratios.nocs[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'] , confint(m)['txsbrt',]))
}
odds.ratios.nocs 
g1.b  <-  make.OR.plot(odds.ratios.nocs, label_list2) + ggtitle('Negative controls')
g1  <- g1.a + g1.b + plot_layout(heights = (c(1,6))) + plot_annotation(title='Treatment effect of SBRT (unadjusted)')
ggsave(g1, width=7, height=3.5, filename = sprintf('figs/additive.raw.%s.pdf', subset.name))


################################
# Now, we adjust for various comorbidities 
################################
adjust.for  <-  setdiff( c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities) , c('race', 'histology', 'LiverSevere', 'HIV', 'Paralysis')) # For now, removing some of the rarer comorbiditis

# Outcomes
hazard.ratios.outcomes.adj  <-  make.odds.ratio.df ( outcome.names[1]) 
for (outcome.i in 1:length(outcome.names)){ 
    outcome.name  <-  outcome.names[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365,
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                          treatment.year = year(tx.date)
        )
    f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + %s + treatment.year',  paste(sprintf('const(%s)', adjust.for), collapse="+") )
    m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
    print(summary(m))
    hazard.ratios.outcomes.adj[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
} 
hazard.ratios.outcomes.adj
g2.a  <-  make.HD.plot(hazard.ratios.outcomes.adj, label_list2) + ggtitle('Outcomes')

# Negative controls
odds.ratios.nocs.adj  <-  make.odds.ratio.df ( noc.names) 
for (outcome.i in 1:length(noc.names)){ 
    outcome.name  <-  noc.names[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( outcome.count  = !!rlang::sym(outcome.name))
    summary( A.temp$time.enrolled, useNA="ifany")
    f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.enrolled))  + %s',  paste(sprintf('%s', adjust.for), collapse="+") )
    m  <-  glm( f  ,  data = A.temp, family = poisson(link=log))
    print(summary(m))
    odds.ratios.nocs.adj[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'] , confint(m)['txsbrt',]))
}
odds.ratios.nocs.adj 
g2.b  <-  make.OR.plot(odds.ratios.nocs.adj, label_list2) + ggtitle('Negative controls')
g2  <- g2.a + g2.b + plot_layout(heights = (c(1,6))) + plot_annotation(title='Treatment effect of SBRT (adjusted)')
ggsave(g2, width=7, height=3.5, filename = sprintf('figs/additive.adj.%s.pdf', subset.name))


################################
# Two step inference, where W is the sum of the counts of each negative control 
################################

# outcomes
set.seed(3)
Z  <-  'optho_any_count'
hazard.ratios.outcomes.two.step  <-  make.odds.ratio.df ( outcome.names) 
num_cores  <-  detectCores()
cl  <-  makeCluster(num_cores)
registerDoParallel(cl)
B  <-  1000
for (outcome.i in 1:length(outcome.names.temp)){ 
    outcome.name  <-  outcome.names[outcome.i]
    noc.names.temp  <- setdiff( noc.names, Z)
    print(outcome.name)
    boot_samples  <-  foreach (Bi = 1:B , .combine = 'rbind', .packages=c('dplyr', 'timereg')) %dopar% {
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        A.temp  <-  A.final2 %>% mutate( 
                              outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                              outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                              W = rowSums( ( across( all_of(noc.names))))
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
        warnings()
        return(boot.res)
    }
   hazard.ratios.outcomes.two.step[outcome.i,1:3]  <-  c( mean(boot_samples[,1]), quantile(boot_samples[,1], c(0.025, 0.975)) )
   print(hazard.ratios.outcomes.two.step[outcome.i,1:3])
}
stopCluster(cl)
hazard.ratios.outcomes.two.step
g3.a  <-  make.HD.plot(hazard.ratios.outcomes.two.step, label_list2) + ggtitle('Outcomes')



# negative controls
set.seed(3)
Z  <-  'optho_any_count'
num_cores  <-  detectCores()
cl  <-  makeCluster(num_cores)
registerDoParallel(cl)
noc.names.temp  <- setdiff( noc.names, Z)
odds.ratios.two.step.nocs  <-  make.odds.ratio.df ( noc.names.temp) 
B  <-  100
for (outcome.i in 1:length(noc.names.temp)){ 
    outcome.name  <-  noc.names.temp[outcome.i]
    noc.names.temp.2  <- setdiff(noc.names.temp, outcome.name)
    print(outcome.name)
    boot_samples  <-  foreach (Bi = 1:B , .combine = 'rbind', .packages=c('dplyr', 'timereg')) %dopar% {
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        A.temp  <-  A.final2 %>% mutate( 
                              W = rowSums( ( across( all_of(noc.names.temp.2))))
        )
        f  <-  sprintf( 'W ~ tx + %s + %s', Z,  paste(sprintf('%s', adjust.for), collapse="+") )
        stage1  <- glm( as.formula(f) , family = poisson(link='log'), data = A.temp)
        negative_outcome_pred  <-  predict(stage1, type = 'link')
        A.temp  <-  A.final2 %>% mutate( 
                                       outcome.count  = !!rlang::sym(outcome.name),
                                        negative_outcome_pred = negative_outcome_pred)
        f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.enrolled))  + %s + negative_outcome_pred',  paste(sprintf('%s', adjust.for), collapse="+") )
        m  <-  glm( f  ,  data = A.temp, family = poisson(link=log))
        print(summary(m))
       boot.res  <-  exp(c( coef(m)['txsbrt'] , confint.default(m, 'txsbrt'))) # to use profile confidence intervals (slower), remove the .default
       return(boot.res)
    }
   odds.ratios.two.step.nocs[outcome.i,1:3]  <-  c( mean(boot_samples[,1]), quantile(boot_samples[,1], c(0.025, 0.975)) )
   print(odds.ratios.two.step.nocs[outcome.i,1:3])
}
stopCluster(cl)



g3.b  <-  make.OR.plot(odds.ratios.two.step.nocs, label_list2) + ggtitle('Negative controls')
g3  <- g3.a + g3.b + plot_layout(heights = (c(1,6))) + plot_annotation(title='Treatment effect of SBRT (two step)')
ggsave(g3, width=7, height=3.5, filename = sprintf('figs/additive.two.step.%s.pdf', subset.name))





