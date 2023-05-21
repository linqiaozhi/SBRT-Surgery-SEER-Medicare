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


################################
#   Not adjusting for X
################################
comorbidities  <-  c('DM','DMcx', 'LiverMild', 'Pulmonary', 'PVD', 'CHF', 'MI', 'Renal', 'Stroke',  'PUD', 'Rheumatic', 'Dementia', 'LiverSevere', 'Paralysis', 't_stage_8')
negative.outcomes.oi  <-  c( 'fall',  'other_injury', 'diverticular_disease', 'hernia',  'arthropathy','GU_sx',  'hpb',  'gout', 'oral', 'optho' )
outcome.names  <-  c( 'death' )
noc.names = sprintf( '%s_pre_count', negative.outcomes.oi )
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

# The outcomes, which right now are only death, will be modeled using Aalen's additive hazards model.
# Eventually, there will be more outcomes, hence the loop structure for only one variable.
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
g1.a  <-  make.OR.plot(hazard.ratios.outcomes, label_list2) + ggtitle('Outcomes')
ggsave(g1.a, width=7, height=3.5, filename = sprintf('figs/cox.raw.%s.pdf', subset.name))


################################
# Now, we adjust for various comorbidities 
################################
adjust.for  <-  setdiff( c('age', 'sex', 'race', 'marital.status', 'histology','treatment.year', comorbidities) , c('race', 'histology', 'LiverSevere', 'HIV', 'Paralysis')) # For now, removing some of the rarer comorbiditis
# Outcomes
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
g2.a  <-  make.OR.plot(hazard.ratios.outcomes.adj, label_list2) + ggtitle('Adjusted for X')
ggsave(g2.a, width=7, height=3.5, filename = sprintf('figs/cox.adj.%s.pdf', subset.name))


################################
# Two step inference, where W is the sum of the counts of each negative control 
################################
source("stats/nc_ph.R")
# outcomes
set.seed(3)
Z  <-  'optho_any_count'
hazard.ratios.outcomes.two.step  <-  make.odds.ratio.df ( outcome.names.temp) 
for (outcome.i in 1:length(outcome.names.temp)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
     # noc.names.temp.no.Z  <- setdiff( noc.names[1:2], 'optho_pre_count')
      noc.names.temp.no.Z  <- setdiff( noc.names, 'optho_pre_count')
    print(outcome.name)
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
 A.temp.mat  <-  A.temp %>% 
        select(tx, all_of(adjust.for), outcome.bool, outcome.time,all_of(Z), all_of(noc.names.temp.no.Z)) %>% 
        mutate_at( all_of( c('treatment.year', 'age')), list(~scale(.))) %>% 
        # mutate_at(all_of(noc.names.temp.no.Z), list(~100*.)) %>% 
        as.matrix()
  m <- cox_nc((A.temp.mat), A = "tx", X = adjust.for,
                        Y = "outcome.time", Z = Z, W = noc.names.temp.no.Z,
                       event = "outcome.bool", ncores = 8, B= 1000)
   hazard.ratios.outcomes.two.step[outcome.i,1:3]  <-  exp(c( m$estimate, m$estimate - 1.96*m$se, m$estimate + 1.96*m$se))
   print(hazard.ratios.outcomes.two.step[outcome.i,1:3])
}

g3.a  <-  make.OR.plot(hazard.ratios.outcomes.two.step, label_list2) + ggtitle('Outcomes')
ggsave(g3.a, width=7, height=3.5, filename = sprintf('figs/cox.proximal.%s.pdf', subset.name))
# ggsave(g3.a, width=7, height=3.5, filename = sprintf('figs/cox.proximal.only2Ws.%s.pdf', subset.name))

# all the Ws
# R> hazard.ratios.outcomes.two.step
#                      estimate    low_ci  high_ci y_axis              outcome
# death                2.187159 1.6282571 2.937906      1                death
# fall                 1.562914 1.3317587 1.834190      2                 fall
# other_injury         1.260193 1.1198181 1.418164      3         other_injury
# diverticular_disease 1.066135 0.9407705 1.208206      4 diverticular_disease
# hernia               1.041045 0.8985794 1.206097      5               hernia
# arthropathy          1.002843 0.8852744 1.136025      6          arthropathy
# GU_sx                1.139847 1.0016330 1.297132      7                GU_sx
# hpb                  1.218061 1.0795600 1.374331      8                  hpb
# gout                 1.299934 0.9883083 1.709820      9                 gout
# oral                 1.097631 0.8419543 1.430949     10                 oral


# hazard.ratios.outcomes.two.step
#                       estimate    low_ci  high_ci y_axis              outcome
# death                1.8245462 1.5403381 2.161194      1                death
# fall                 1.4277058 1.2540791 1.625371      2                 fall
# other_injury         1.1494168 1.0649926 1.240533      3         other_injury
# diverticular_disease 1.0170611 0.9059234 1.141833      4 diverticular_disease
# hernia               0.9882388 0.8692978 1.123454      5               hernia
# arthropathy          0.9894468 0.9038633 1.083134      6          arthropathy
# GU_sx                1.0205371 0.9366934 1.111886      7                GU_sx
# hpb                  1.2268532 1.1088004 1.357475      8                  hpb
# gout                 1.1419000 0.9344278 1.395437      9                 gout
# oral                 1.1564056 0.9229923 1.448846     10                 oral
