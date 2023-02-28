library(dplyr)
library('ggtext')
library(arsenal)
library(ggplot2)
source('utilities.R')
#subset.name <- 'age.gte.80'
subset.name <- 'all'
filename.in  <-  sprintf('data/A.final.%s.2.RDS', subset.name)
A.final  <-  readRDS(filename.in)
summary( A.final$age)
label_list  <-  readRDS('data/label.list.RDS')
table( nna(A.final$hemorrhoids), useNA="ifany")
################################
#   Logistic regression 
################################
comorbidities  <-  c('DM','DMcx', 'LiverMild', 'Pulmonary', 'PVD', 'CHF', 'MI', 'Renal', 'Stroke',  'PUD', 'Rheumatic', 'Dementia', 'LiverSevere', 'Paralysis', 'HIV', 'Smoking', 'o2')
#negative.outcomes.oi  <-  c( 'fall', 'cholelithiasis', 'diverticular_disease', 'hernia', 'hemorrhoids', 'GU_sx', 'arthropathy', 'hpb' )
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
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F))
    m  <- glm( (outcome.bool) ~ tx +  offset( log(outcome.time) ) , data = A.temp, family = poisson(link=log))
    #print(summary(m))
    odds.ratios[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'], confint(m,'txsbrt'))) 
}
odds.ratios
g1  <-  make.OR.plot(odds.ratios, label_list2)+ ggtitle('A) Treatment effect of SBRT (unadjusted)')
ggsave(g1, width=7, height=2, filename = sprintf('figs/regression.raw.%s.pdf', subset.name))


# adjusting

odds.ratios.adj  <-  make.odds.ratio.df ( outcome.names) 
adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities)
adjust.for  <-  setdiff( c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities) , c('race', 'histology', 'LiverSevere', 'HIV', 'Paralysis'))
print('WARNING: Removing columns')
for (outcome.i in 1:length(outcome.names)){ 
    outcome.name  <-  outcome.names[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F))
    f  <-  sprintf( 'outcome.bool ~ tx + %s + offset( log(outcome.time) )',  paste(adjust.for, collapse="+") )
    f
    m  <- glm( as.formula(f), data = A.temp, family = poisson(link=log))
    print(summary(m))
    odds.ratios.adj[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'], confint(m,'txsbrt'))) 
    #A.temp.colsfiltered  <-  A.temp[, lapply(1:ncol(A.temp), FUN = function(x) min( table(as.data.frame(A.temp)[,x], as.data.frame(A.temp)[,'outcome.bool'])  )) >= 1]
    #print(setdiff(adjust.for, colnames(A.temp.colsfiltered)))
}
odds.ratios.adj
g2  <-  make.OR.plot(odds.ratios.adj, label_list2)+ ggtitle('B) Treatment effect of SBRT (adjusted)')
ggsave(g2 , width=7, height=2, filename =sprintf('figs/regression.adj.%s.pdf', subset.name))


################################
# two step sandbox 
################################
adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities)
#adjust.for  <-  setdiff( c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities) , c('race', 'histology', 'LiverSevere', 'HIV', 'Paralysis'))
#print('WARNING: Removing columns')
#devtools::install_url('https://cran.r-project.org/src/contrib/Archive/instruments/instruments_0.1.0.tar.gz')
#outcome.names.temp  <-  setdiff(outcome.names, c('fall', 'hernia'))
#outcome.names.temp  <-  setdiff(outcome.names, c('fall', 'hernia'))
Z  <-  'optho'
W  <-   'GU_sx'
outcome.names.temp  <-  setdiff(outcome.names, c(W, Z))
odds.ratios.adj  <-  make.odds.ratio.df ( outcome.names.temp) 
for (outcome.i in 1:length(outcome.names.temp)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F))
    #stage1  <- lm( 'nna(fall) ~ tx + nna(hernia)', data = A.temp)
    #stage1  <- lm( 'nna(diverticular_disease) ~ tx + nna(arthropathy)', data = A.temp)
    #stage1  <- lm( 'nna(arthropathy) ~ tx + nna(optho)', data = A.temp)
    #stage1  <- lm( sprintf('nna(diverticular_disease) ~ tx +%s + nna(optho)', paste(adjust.for, collapse="+")) , data = A.temp)
    stage1  <- lm( sprintf('nna(%s) ~ tx + nna(%s)', W, Z) , data = A.temp)
    #stage1  <- lm( sprintf('nna(%s) ~ tx + nna(%s)+ %s', W, Z, paste(adjust.for, collapse="+")) , data = A.temp)
    #stage1  <- lm( 'nna(hpb) ~ tx + nna(optho)', data = A.temp)
    negative_outcome_pred  <-  predict(stage1)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F))
    f  <-  sprintf( 'outcome.bool ~ tx + negative_outcome_pred + %s + offset( log(outcome.time) )',  paste(adjust.for, collapse="+") )
    f
    m  <- glm( as.formula(f), data = cbind(A.temp, negative_outcome_pred), family = poisson(link=log))
    print(summary(m))
    odds.ratios.adj[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'], confint(m,'txsbrt'))) 
    print(odds.ratios.adj[outcome.i,1:3])
}
odds.ratios.adj
g2  <-  make.OR.plot(odds.ratios.adj, label_list2)+ ggtitle(sprintf('C) Treatment effect of SBRT (proximal)\nW: %s, Z: %s', label_list2[W],label_list2[Z]))
ggsave(g2 , width=7, height=2, filename =sprintf('figs/regression.prox.%s.pdf', subset.name))


################################
# Using just two variables
################################
outcome.names.temp  <-  setdiff(outcome.names, c('diverticular_disease', 'optho'))
odds.ratios.expadj  <-  make.odds.ratio.df ( outcome.names.temp) 
for (outcome.i in 1:length(outcome.names.temp)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
     #adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities, sprintf('nna(%s_pre)', c('GU_sx', 'optho' ))
    #adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities, sprintf('nna(%s_pre)', negative.outcomes.oi ))
    adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities, sprintf('nna(%s)', c('diverticular_disease', 'optho') ))
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F))
    f  <-  sprintf( 'outcome.bool ~ tx + %s + offset( log(outcome.time) )',  paste(adjust.for, collapse="+") )
    f
    m  <- glm( as.formula(f), data = A.temp, family = poisson(link=log))
    print(summary(m))
    odds.ratios.expadj[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'], confint(m,'txsbrt'))) 
}
odds.ratios.expadj
g3  <-  make.OR.plot(odds.ratios.expadj, label_list2)+ ggtitle('C) Treatment effect of SBRT (expanded adjustment)') 
ggsave(g3 , width=7, height=2, filename =sprintf('figs/regression.expadj.2.%s.pdf', subset.name))






#    A.temp  <-  A.final %>% mutate( 
#                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
#                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
#                          neg_outcome = nna(diverticular_disease),
#                          neg_exposure = nna(arthropathy), 
#                          tx = ifelse( tx == 'sublobar', T, F))
#    
#    A.temp.2  <-  A.temp[,c('outcome.bool', 'neg_outcome', 'tx', 'neg_exposure', 'outcome.time')]
#    sum(is.na(A.temp.2))
#    table( A.temp$diverticular_disease, useNA="ifany")
#    f2  <-  sprintf( 'neg_outcome ~ tx +  neg_exposure',paste(adjust.for, collapse="+"))
#    stage1  <- lm( f2, A.temp.2)
#    neg_outcome_pred  <- predict(stage1)
#    f  <-  sprintf( 'outcome.bool ~ tx + neg_outcome_pred  ',  paste(adjust.for, collapse="+") )
#    #stage2  <- glm( as.formula(f),  data = A.temp.2, family = poisson(link=log))
#    stage2  <- glm( as.formula(f),  data = cbind( A.temp.2, neg_outcome_pred), family = 'binomial')
#    summary(stage2)
#    #m  <-  ivglm(estmethod = 'ts', fitX.LZ = stage1, fitY.LX = stage2, data = (A.temp.2), ctrl = T )
#    A.temp.2[8213:8216,]
#
#
#    traceback()
#
#    #find_instruments( f,f2)
#
#    print(summary(m))
#    odds.ratios.adj[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'], confint(m,'txsbrt'))) 
#    #A.temp.colsfiltered  <-  A.temp[, lapply(1:ncol(A.temp), FUN = function(x) min( table(as.data.frame(A.temp)[,x], as.data.frame(A.temp)[,'outcome.bool'])  )) >= 1]
#    #print(setdiff(adjust.for, colnames(A.temp.colsfiltered)))




# Diagnostics to remove columns which have <3 observations for a level

################################
# Leaving one out 
################################
odds.ratios.expadj  <-  make.odds.ratio.df ( outcome.names) 
for (outcome.i in 1:length(outcome.names)){ 
    outcome.name  <-  outcome.names[outcome.i]
     adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities, sprintf('nna(%s_pre)', negative.outcomes.oi[ negative.outcomes.oi != outcome.name] ))
    #adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities, sprintf('nna(%s_pre)', negative.outcomes.oi ))
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F))
    f  <-  sprintf( 'outcome.bool ~ tx + %s + offset( log(outcome.time) )',  paste(adjust.for, collapse="+") )
    f
    m  <- glm( as.formula(f), data = A.temp, family = poisson(link=log))
    print(summary(m))
    odds.ratios.expadj[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'], confint(m,'txsbrt'))) 
}
odds.ratios.expadj
g3  <-  make.OR.plot(odds.ratios.expadj, label_list2)+ ggtitle('C) Treatment effect of SBRT (expanded adjustment)') 
ggsave(g3 , width=7, height=2, filename =sprintf('figs/regression.expadj.%s.pdf', subset.name))


library(patchwork)
g  <-  g1 / g2 / g3
ggsave(g , width=7, height=6, filename =sprintf('figs/regression.%s.pdf', subset.name))

table( nna(A.final$death), useNA="ifany")



