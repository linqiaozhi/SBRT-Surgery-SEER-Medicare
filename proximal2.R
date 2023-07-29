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

################################
# Load data 
################################
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
# Full set of X
 X.factor  <-  c('t_stage_8','sex', 'race', 'marital.status', 'histology' ) 
 X.numeric  <-  c('age', 'treatment.year', sprintf('%s_pre_count', c(
                     # DMEs
                       'hospital_beds_and_supplies', 'wheelchairs_accessories', 'walking_aids', 'O2accessories', 'other_supplies', 'diabetic_footwear', 'transportation_services', 
                     # Diagnoses
                        'smoking', 'o2', 'other_bacterial_diseases', 'pneumonia_and_influenza', 'pressure_ulcer', 'ischemic_heart_disease', 'CHF', 'PVD', 'CVD', 'dementia', 'COPD', 'PUD', 'MILDLD', 'DIAB_UC', 'DIAB_C', 'PARA', 'RD', 'cancer_nonlung', 'MSLD', 'METS',  'mental_disorders', 'nervous_system', 'other_heart_disease', 'veins_lymphatics_other_circulatory', 'rheum',
                     # Drugs
                     'Insulin', 'Anticoags')
                         )) 
# Minimal set of X
# X.factor  <-  c('t_stage_8','sex', 'race',  'histology' ) 
# X.numeric  <-  c('age', 'treatment.year', sprintf('%s_pre_count', c(
#                     # Diagnoses
#                        'smoking', 'o2',  'ischemic_heart_disease', 'CHF', 'dementia', 'COPD', 'DIAB_C', 'cancer_nonlung', 'MSLD', 'METS')
#                         )) 

# Get number of non-zero values for each X.numeric variable
 A.final %>% select(all_of(X.numeric)) %>% summarise_all(list(~sum(. > 0, na.rm = T))) %>% t() %>% as.data.frame() %>% arrange(desc(V1))

# Use dplyr to scale all variables in the numeric.X list
scale_  <-  function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T)
# Make each X.numeric variable into a z-score, but name them with suffix _z
A.final  <- A.final %>% mutate( 
    across( all_of(X.numeric), scale_ , .names = "{.col}_z" )
)

negative.outcomes.oi  <-  c( 'fall',  'other_injury', 'diverticular_disease', 'hernia',  'arthropathy','GU_sx', 'optho' )
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
Zs  <-  c('O2accessories_pre_count_z', 'walking_aids_pre_count_z' , 'hospital_beds_and_supplies_pre_count_z' , 'wheelchairs_accessories_pre_count_z' , 'transportation_services_pre_count_z', 'other_supplies_pre_count_z', 'diabetic_footwear_pre_count_z' )


X.factor.min  <-  c('t_stage_8','sex', 'race',  'histology' ) 
X.numeric.min  <-  c('age', 'treatment.year', sprintf('%s_pre_count', c( 'smoking', 'o2',  'ischemic_heart_disease', 'CHF', 'dementia', 'COPD', 'DIAB_C', 'cancer_nonlung', 'MSLD', 'METS'))) 
adjust.for.scaled  <-  setdiff( c(sprintf('%s_z', X.numeric.min), X.factor.min) , Zs) 
# adjust.for.scaled  <-  setdiff( c(sprintf('%s_z', X.numeric), X.factor) , c()) 
B  <-  1000
for (i in 1:length(adjust.for.scaled)) cat(i, adjust.for.scaled[i], '\n')

################################
# Proximal  function 
################################
two.step  <-  function(A.final2, Zs,  outcome.name,noc.names.temp, adjust.for.scaled, Y.count = F, verbose = F){
    # Stage 1: W1 is non-cause mortality, W2 is the sum of negative control outcomes (excluding the negative outcome of interest)
    W1  <-  'death.other.cause'
    A.temp  <-  A.final2 %>% mutate( 
                                    W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt)/365,
                                    W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
                                    W2 = rowSums( ( across( all_of(noc.names.temp)))),
    )
    # W1
    f  <-  sprintf( 'Surv(W1.time, W1.bool) ~ const(tx) +%s+ %s',  
                   paste(sprintf('const(%s)', Zs), collapse="+"), 
                   paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
    m1  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0, silent = 0 )
    if (verbose) 
        print(summary(m1))
    mm  <-  model.matrix(as.formula(f), A.temp)[,-1] 
    coefs  <-  as.matrix(coef( m1)[,'Coef.'])
    negative_outcome_pred1  <-  mm %*% coefs
    # W2
    f  <-  sprintf( 'W2 ~ tx + %s + %s', 
                   paste(sprintf('const(%s)', Zs), collapse="+"),  
                   paste(sprintf('%s', adjust.for.scaled), collapse="+") )
    m2  <- glm( as.formula(f) , family = poisson(link='log'), data = A.temp)
    if (verbose) 
        print(summary(m2))
    negative_outcome_pred2  <-  predict(m2, type = 'link')
    A.temp  <-  A.temp %>% mutate( 
                                  What1 = scale(negative_outcome_pred1[,1]),
                                  What2 = scale(negative_outcome_pred2)
    )
    # Stage 2
    # Poisson for count Ys
    if (Y.count) {
        A.temp  <- A.temp %>% mutate( Y.count  = !!rlang::sym(outcome.name),)
        f  <-  sprintf( 'Y.count ~ tx + offset(log(time.enrolled))  + %s+ What1 + What2',  
        # f  <-  sprintf( 'Y.count ~ tx + offset(log(time.enrolled))  + %s',  
                       paste(sprintf('%s', adjust.for.scaled), collapse="+") )
        m  <-  glm( f  ,  data = A.temp, family = poisson(link=log))
        if (verbose) 
            print(summary(m))
        boot.res  <-  exp(c( coef(m)['txsbrt'] ))
    }else{
        # Aalen's for time-to-event Ys
        A.temp  <- A.temp %>% mutate(
                                     Y.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                                     Y.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        )
        f  <-  sprintf( 'Surv(Y.time, Y.bool) ~ const(tx) + const(What1) + const(What2) +  %s',  
                       paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
        m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0, silent = 0)
        if (verbose) 
            print(summary(m))
        boot.res  <-  coef(m)['const(tx)sbrt', 'Coef.']
    }
    return(boot.res)
}

################################
# Time-to-event outcomes 
################################

hazard.differences.outcomes  <-  make.odds.ratio.df ( outcome.names) 
hazard.differences.outcomes.adj  <-  make.odds.ratio.df ( outcome.names) 
hazard.differences.outcomes.proximal  <-  make.odds.ratio.df ( outcome.names) 
for (outcome.i in 1:length(outcome.names)){ 
    outcome.name  <-  outcome.names[outcome.i]
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365,
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        )
    # Raw
    m  <-  aalen( Surv(outcome.time, outcome.bool) ~ const(tx) ,  data = A.temp, robust = 0)
    hazard.differences.outcomes[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
    # Adjust for X
    f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + %s',  paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
    m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
    hazard.differences.outcomes.adj[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
    # Proximal
    mout  <-  two.step(A.final, Zs,  outcome.name, noc.pre.count.names, adjust.for.scaled)
    est_boot <- parallel::mclapply(1:B, function(bb){
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        mout  <-  two.step(A.final2, Zs,  outcome.name, noc.pre.count.names, adjust.for.scaled)
        return(mout)
    }, mc.cores =8)
    se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
   hazard.differences.outcomes.proximal[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
   print(hazard.differences.outcomes.proximal[outcome.i,1:3])
} 

################################
#  Counts
################################

odds.ratios.nocs  <-  make.odds.ratio.df ( noc.pre.count.names) 
odds.ratios.nocs.adj  <-  make.odds.ratio.df ( noc.pre.count.names) 
odds.ratios.nocs.proximal  <-  make.odds.ratio.df ( noc.pre.count.names) 
for (outcome.i in 1:length(noc.pre.count.names)){ 
    outcome.name  <-  noc.pre.count.names[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( outcome.count  = !!rlang::sym(outcome.name))
    # Raw
    f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.enrolled))' )
    m  <-  glm( f  ,  data = A.temp, family = poisson(link=log))
    odds.ratios.nocs[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'] , confint(m)['txsbrt',]))
    # Adjust for X
    f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.enrolled))  + %s',  paste(sprintf('%s', adjust.for.scaled), collapse="+") )
    m  <-  glm( f  ,  data = A.temp, family = poisson(link=log))
    odds.ratios.nocs.adj[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'] , confint(m)['txsbrt',]))
    # Proximal
    noc.names.temp  <- setdiff( noc.pre.count.names, c(outcome.name))
    mout  <-  two.step(A.final, Zs,  outcome.name, noc.names.temp, adjust.for.scaled,Y.count =T, verbose=F)
    est_boot <- parallel::mclapply(1:B, function(bb){
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        boot.res  <-  two.step(A.final2, Zs,  outcome.name, noc.names.temp,adjust.for.scaled, Y.count=T, verbose=F)
        return(boot.res)
    }, mc.cores =8)
    se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
    odds.ratios.nocs.proximal[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
    print(odds.ratios.nocs.proximal[outcome.i,1:3])
}




################################
# Plot 
################################

g1.a  <-  make.HD.plot(hazard.differences.outcomes, label_list2)
g1.b  <-  make.OR.plot(odds.ratios.nocs, label_list2)
g1  <-  g1.a / g1.b+ plot_layout(heights = (c(1,2))) + plot_annotation(title="Raw")
ggsave(g1, width=7, height=3, filename = sprintf('figs/raw.pdf'))
g2.a  <-  make.HD.plot(hazard.differences.outcomes.adj, label_list2)
g2.b  <-  make.OR.plot(odds.ratios.nocs.adj, label_list2)
g2  <-  g2.a / g2.b+ plot_layout(heights = (c(1,2)))+ plot_annotation(title="Adj")
ggsave(g2, width=7, height=3, filename = sprintf('figs/adj.pdf'))
g3.a  <-  make.HD.plot(hazard.differences.outcomes.proximal, label_list2)
g3.b  <-  make.OR.plot(odds.ratios.nocs.proximal, label_list2)
g3  <-  g3.a / g3.b + plot_layout(heights = (c(1,2)))+ plot_annotation(title="Proximal")
ggsave(g3, width=7, height=3, filename = sprintf('figs/proximal.pdf'))


#################################
## Assess quality of our negative outcomes 
#################################
#noc.i  <-  9
#for (noc.i in 1:length(noc.pre.count.names)){ 
#    noc.name  <-  noc.pre.count.names[noc.i]
#    Y  <-  'death'
#    A.temp  <-  A.final %>% mutate( 
#                          outcome.time  = if_else (nna(!!rlang::sym(Y)), as.numeric( !!rlang::sym(Y) - tx.date, units = "days" ), tt)/ 365,
#                          outcome.bool = ifelse( nna(!!rlang::sym(Y)), T, F),
#                          noc = !!rlang::sym(noc.name)
#        )
#    # Raw
#    m  <-  aalen( Surv(outcome.time, outcome.bool) ~ const(noc) +treatment.year_z,  data = A.temp, robust = 0)
#    print(noc.name)
#    print( c( coef(m)['const(noc)', c('Coef.', 'lower2.5%', 'upper97.5%')]) )
#    # Adjust for X
#    f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(noc) + %s',  paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
#    m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
#    print( c( coef(m)['const(noc)', c('Coef.', 'lower2.5%', 'upper97.5%')]) )
#}


#################################
## Figure out confidence intervals 
#################################
#A.final %>% select(all_of(noc.pre.count.names)) %>% summarise_all(list(~sum(. > 0, na.rm = T))) %>% t() %>% as.data.frame() %>% arrange(desc(V1))

#outcome.i  <-  2
#outcome.name  <-  noc.pre.count.names[outcome.i]
#print(outcome.name)
#A.temp  <-  A.final %>% mutate( outcome.count  = !!rlang::sym(outcome.name))
#print(table( A.temp2$outcome.count, useNA="ifany"))
#f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.enrolled))  + %s',  paste(sprintf('%s', adjust.for.scaled), collapse="+") )
#m  <-  glm( f  ,  data = A.temp, family = poisson(link=log))
#point.est  <-  coef(m)['txsbrt']
#ci  <-  confint(m)['txsbrt',]

#est_boot <- parallel::mclapply(1:B, function(bb){
#    A.temp2  <-  A.temp[sample(nrow(A.temp),replace=T ),]
#    m  <-  glm( f  ,  data = A.temp2, family = poisson(link=log))
#    return(coef(m)['txsbrt'])
#}, mc.cores =8)
#se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
#ci2 <-  c( point.est - 1.96*se, point.est + 1.96*se )
#print(ci)
#print(ci2)



#var(A.temp$outcome.count)
#mean(A.temp$outcome.count)

#f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.enrolled))  + %s',  paste(sprintf('%s', adjust.for.scaled), collapse="+") )
#m  <-  MASS::glm.nb( f  ,  data = A.temp)
#point.est  <-  coef(m)['txsbrt']
#ci  <-  confint(m)['txsbrt',]

#est_boot <- parallel::mclapply(1:B, function(bb){
#    A.temp2  <-  A.temp[sample(nrow(A.temp),replace=T ),]
#    # print(table( A.temp2$pancreatic_pre_count > 0, useNA="ifany"))
#    m  <-  MASS::glm.nb( f  ,  data = A.temp2)
#    return(coef(m)['txsbrt'])
#}, mc.cores =8)
#se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
#ci2 <-  c( point.est - 1.96*se, point.est + 1.96*se )
#print(ci)
#print(ci2)

