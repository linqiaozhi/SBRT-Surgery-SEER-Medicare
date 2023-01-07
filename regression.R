library(tidyverse)
library(arsenal)
source('utilities.R')
subset.name <- 'gte.80'
subset.name <- 'all'
filename.in  <-  sprintf('data/A.final.%s.RDS', subset.name)
A.final  <-  readRDS(filename.in)
label_list  <-  readRDS('data/label.list.RDS')
################################
#   Logistic regression #TODO: (will go in separate file)
################################
negative.outcomes.oi  <-  c( 'fall', 'cholelithiasis', 'diverticular_disease', 'hernia', 'hemorrhoids', 'GU_sx', 'arthropathy', 'hpb' )
outcome.names  <-  c( 'death',negative.outcomes.oi ) 
label_list2  <-  c( label_list,
                    death = 'Death', 
                    fall = 'Fall',
                    GU_sx = 'Genitourinary-related',
                    arthropathy = 'Arthropathy',
                    cholelithiasis = 'Cholelithiasis-related',
                    gout = 'Gout',
                    obstruction = 'Intestinal obstruction',
                    hernia = 'Abdominal hernia',
                    diverticular_disease = 'Diverticular disease',
                    hemorrhoids = 'Hemorrhoids',
                    hpb = 'HPB-related'
)

tblcontrol <- tableby.control(numeric.stats = c('Nmiss', 'meansd'), numeric.simplify = T, cat.simplify =T, digits = 1,total = T,test = F)
f  <-  sprintf( 'tx ~ %s', paste( sprintf('%s', outcome.names), collapse = "+") )
A.final.toprint  <- A.final %>% mutate( across(outcome.names, ~ nna(.x)))
labels(A.final.toprint)  <-  label_list2
tt <- tableby(as.formula(f), data=A.final.toprint, control = tblcontrol)
summary(tt) %>% write2html('/PHShome/gcl20/Research_Local/SEER-Medicare/tbls/outcomes.htm')


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
g  <-  make.OR.plot(odds.ratios, label_list2)
ggsave(g, width=7, height=2, filename = sprintf('figs/regression.raw.%s.pdf', subset.name))


# adjusting

odds.ratios.adj  <-  make.odds.ratio.df ( outcome.names) 
adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities)
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
}
odds.ratios.adj
g  <-  make.OR.plot(odds.ratios.adj, label_list2)
ggsave(g + ggtitle('Treatment effect of SBRT (adjusted)'), width=7, height=2, filename =sprintf('figs/regression.adj.%s.pdf', subset.name))


################################
# Leaving one out 
################################
odds.ratios.expadj  <-  make.odds.ratio.df ( outcome.names) 
for (outcome.i in 1:length(outcome.names)){ 
    outcome.name  <-  outcome.names[outcome.i]
    # adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities, sprintf('nna(%s_pre)', negative.outcomes.oi[ negative.outcomes.oi != outcome.name] ))
    adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities, sprintf('nna(%s_pre)', negative.outcomes.oi ))
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
g  <-  make.OR.plot(odds.ratios.expadj, label_list2)
ggsave(g + ggtitle('Treatment effect of SBRT (expanded adjustment)'), width=7, height=2, filename =sprintf('figs/regression.expadj.%s.pdf', subset.name))

