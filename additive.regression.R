library(dplyr)
library(foreach)
library(doParallel)
library(timereg)
library('ggtext')
library(arsenal)
library(ggplot2)
source('utilities.R')
set.seed(3)
#subset.name <- 'age.gte.80'
subset.name <- 'all'
filename.in  <-  sprintf('data/A.final.%s.2.RDS', subset.name)
A.final  <-  readRDS(filename.in) %>% filter ( race != 'Unknown')
summary( A.final$age)
label_list  <-  readRDS('data/label.list.RDS')
table( nna(A.final$hemorrhoids), useNA="ifany")
################################
#   Logistic regression 
################################
comorbidities  <-  c('DM','DMcx', 'LiverMild', 'Pulmonary', 'PVD', 'CHF', 'MI', 'Renal', 'Stroke',  'PUD', 'Rheumatic', 'Dementia', 'LiverSevere', 'Paralysis',  'Smoking', 'o2')
negative.outcomes.oi  <-  c( 'fall',  'other_injury', 'diverticular_disease', 'hernia',  'arthropathy','GU_sx',  'hpb', 'optho' )
outcome.names  <-  c( 'death',negative.outcomes.oi ) 
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
                    optho = '*Ophthalmic*'
)

tblcontrol <- tableby.control(numeric.stats = c('Nmiss', 'meansd'), numeric.simplify = T, cat.simplify =T, digits = 1,total = T,test = F)
f  <-  sprintf( 'tx ~ %s', paste( sprintf('%s', outcome.names), collapse = "+") )
A.final.toprint  <- A.final %>% mutate( across(outcome.names, ~ nna(.x)))
labels(A.final.toprint)  <-  label_list2
tt <- tableby(as.formula(f), data=A.final.toprint, control = tblcontrol)
summary(tt) %>% write2html(sprintf( '%s/outcomes.htm', tools::file_path_as_absolute('tbls') ))


odds.ratios  <-  make.odds.ratio.df ( outcome.names) 
# Contributing person time until 1) time of outcome, 2) death, or 3) end of
# follow up. The latter two are in tt
outcome.i  <-  1
for (outcome.i in 1:length(outcome.names)){ 
    outcome.name  <-  outcome.names[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365,
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F))
    m  <-  aalen( Surv(outcome.time, outcome.bool) ~ 1+ const(tx)  ,  data = A.temp, robust = 0)
    #print(summary(m))
    odds.ratios[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
}
odds.ratios
g1  <-  make.HD.plot(odds.ratios, label_list2)+ ggtitle('A) Treatment effect of SBRT (unadjusted)')
ggsave(g1, width=7, height=2, filename = sprintf('figs/additive.raw.%s.pdf', subset.name))

# adjusting

odds.ratios.adj  <-  make.odds.ratio.df ( outcome.names) 
adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities)
adjust.for  <-  setdiff( c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities) , c('race', 'histology', 'LiverSevere', 'HIV', 'Paralysis'))
print('WARNING: Removing columns')
for (outcome.i in 1:length(outcome.names)){ 
    outcome.name  <-  outcome.names[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F))
    f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + %s',  paste(sprintf('const(%s)', adjust.for), collapse="+") )
    m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
    odds.ratios.adj[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
}
odds.ratios.adj
g2  <-  make.HD.plot(odds.ratios.adj, label_list2)+ ggtitle('B) Treatment effect of SBRT (adjusted)')
ggsave(g2 , width=7, height=2, filename =sprintf('figs/additive.adj.%s.pdf', subset.name))




################################
# two step sandbox, with W compiled from multiple 
# Bootstrap edition!
################################
adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities)
Z  <-  'optho'
outcome.names.temp  <-  setdiff(outcome.names, Z)
odds.ratios.adj  <-  make.odds.ratio.df ( outcome.names.temp) 
outcome.i = 1 
num_cores  <-  detectCores()
cl  <-  makeCluster(num_cores)
registerDoParallel(cl)
B  <-  1000
for (outcome.i in 1:length(outcome.names.temp)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    outcome.names.temp2  <-  setdiff(outcome.names.temp, outcome.name)
    print(outcome.name)
    boot_samples  <-  foreach (Bi = 1:B , .combine = 'rbind', .packages=c('dplyr', 'timereg')) %dopar% {
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        A.temp  <-  A.final2 %>% mutate( 
                              outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                              outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                              W = rowSums( !is.na( across( all_of(outcome.names.temp2))))
        )
        f  <-  sprintf( 'W ~ tx + nna(%s) + %s', Z,  paste(sprintf('const(%s)', adjust.for), collapse="+") )
        table( A.temp$W, useNA="ifany")
        stage1  <- glm( as.formula(f) , family = poisson(link='log'), data = A.temp)
        negative_outcome_pred  <-  predict(stage1, type = 'link')
        A.temp  <-  A.final2 %>% mutate( 
                              outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                              outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F))
        f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + negative_outcome_pred +  %s',  paste(sprintf('const(%s)', adjust.for), collapse="+") )
        m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
        boot.res  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')])
        return(boot.res)
    }
   odds.ratios.adj[outcome.i,1:3]  <-  c( mean(boot_samples[,1]), quantile(boot_samples[,1], c(0.025, 0.975)) )
   print(odds.ratios.adj[outcome.i,1:3])
}
stopCluster(cl)
odds.ratios.adj
g2  <-  make.HD.plot(odds.ratios.adj, label_list2)+ ggtitle(sprintf('C) Treatment effect of SBRT (proximal)\n Z: %s', label_list2[Z]))
ggsave(g2 , width=7, height=2, filename =sprintf('figs/additive.prox3.%s.pdf', subset.name))


#################################
## two step sandbox, with W compiled from multiple 
##No bootstrap
#################################
#adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities)
#Z  <-  'optho'
#outcome.names.temp  <-  setdiff(outcome.names, Z)
#odds.ratios.adj  <-  make.odds.ratio.df ( outcome.names.temp) 
#outcome.i = 1 
#for (outcome.i in 1:length(outcome.names.temp)){ 
#    outcome.name  <-  outcome.names.temp[outcome.i]
#    outcome.names.temp2  <-  setdiff(outcome.names.temp, outcome.name)
#    print(outcome.name)
#    print(outcome.names.temp2)
#    #A  <- A.final %>% mutate(W = rowSums( !is.na( across( all_of(outcome.names.temp)))))
#    A.temp  <-  A.final %>% mutate( 
#                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
#                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
#                          W = rowSums( !is.na( across( all_of(outcome.names.temp2))))
#    )
#    f  <-  sprintf( 'W ~ tx + nna(%s) + %s', Z,  paste(sprintf('const(%s)', adjust.for), collapse="+") )
#    table( A.temp$W, useNA="ifany")
#    stage1  <- glm( as.formula(f) , family = poisson(link='log'), data = A.temp)
#    negative_outcome_pred  <-  predict(stage1, type = 'link')
#    A.temp  <-  A.final %>% mutate( 
#                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
#                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F))
#    f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + negative_outcome_pred +  %s',  paste(sprintf('const(%s)', adjust.for), collapse="+") )
#    m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
#    #fofo  <-  predict.timereg(m, resample.iid=T, newdata = A.temp)
#    odds.ratios.adj[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
#    print(odds.ratios.adj[outcome.i,1:3])
#}
#odds.ratios.adj
#g2  <-  make.HD.plot(odds.ratios.adj, label_list2)+ ggtitle(sprintf('C) Treatment effect of SBRT (proximal)\n Z: %s', label_list2[Z]))
#ggsave(g2 , width=7, height=2, filename =sprintf('figs/additive.prox2.%s.pdf', subset.name))

