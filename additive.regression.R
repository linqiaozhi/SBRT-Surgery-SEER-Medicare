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
filename.in  <-  sprintf('data/A.final.%s.5.RDS', subset.name)
A.final  <-  readRDS(filename.in) %>% filter ( race != 'Unknown')
A.final %>% count(tx)
table( nna(readRDS(filename.in)$optho), useNA="ifany")
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
                    optho2 = '*Ophthalmic*'
)

tblcontrol <- tableby.control(numeric.stats = c('Nmiss', 'meansd'), numeric.simplify = T, cat.simplify =T, digits = 1,total = T,test = F)
f  <-  sprintf( 'tx ~ %s', paste( sprintf('%s', outcome.names), collapse = "+") )
A.final.toprint  <- A.final %>% mutate( across(outcome.names, ~ nna(.x)))
labels(A.final.toprint)  <-  label_list2
tt <- tableby(as.formula(f), data=A.final.toprint, control = tblcontrol)
summary(tt) %>% write2html(sprintf( '%s/outcomes.htm', tools::file_path_as_absolute('tbls') ))

max(A.final$death, na.rm = T)

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
set.seed(3)
adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities)
# Z  <-  'optho2_any_date_count'
Z  <-  'optho_any_date_count'
#Z  <-  'gout_post_count'
outcome.names.temp  <-  setdiff(outcome.names, c('optho', Z))
odds.ratios.adj  <-  make.odds.ratio.df ( outcome.names.temp) 
outcome.i = 3 
num_cores  <-  detectCores()
cl  <-  makeCluster(num_cores)
registerDoParallel(cl)
B  <-  1000
for (outcome.i in 1:length(outcome.names.temp)){ 
    outcome.name  <-  outcome.names.temp[outcome.i]
    #outcome.names.temp2  <-  setdiff(outcome.names.temp, c(outcome.name)
    outcome.names.temp2  <-  setdiff(outcome.names.temp, c(outcome.name, 'fall', 'other_injury'))
    print(outcome.name)
    boot_samples  <-  foreach (Bi = 1:B , .combine = 'rbind', .packages=c('dplyr', 'timereg')) %dopar% {
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        A.temp  <-  A.final2 %>% mutate( 
                              outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                              outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                              W = rowSums( !is.na( across( all_of(outcome.names.temp2))))
        )
        #f  <-  sprintf( 'W ~ tx + nna(%s) + %s', Z,  paste(sprintf('const(%s)', adjust.for), collapse="+") )
        f  <-  sprintf( 'W ~ tx + %s + %s', Z,  paste(sprintf('%s', adjust.for), collapse="+") )
        #f  <-  sprintf( 'W ~ tx +optho2_post_count + gout_post_count + %s',  paste(sprintf('const(%s)', adjust.for), collapse="+") )
        # f  <-  sprintf( 'W ~ tx + %s', Z   )
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
   odds.ratios.adj[outcome.i,1:3]  <-  c( mean(boot_samples[,1]), quantile(boot_samples[,1], c(0.025, 0.975)) )
   print(odds.ratios.adj[outcome.i,1:3])
}
stopCluster(cl)
odds.ratios.adj



g2  <-  make.HD.plot(odds.ratios.adj, label_list2)+ ggtitle(sprintf('C) Treatment effect of SBRT (proximal)\n Z: %s', label_list2[Z]))
ggsave(g2 , width=7, height=2, filename =sprintf('figs/additive.prox4.%s.pdf', subset.name))




################################
# Diagnostic 
################################
# Diagnostics
adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities)
Z  <-  'optho_any_date_count'
outcome.i = 4 
outcome.name  <-  outcome.names.temp[outcome.i]
A.final2  <-  A.final
A.temp  <-  A.final2 %>% mutate( 
                      outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                      outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                      W = rowSums( !is.na( across( all_of(outcome.names.temp2))))
)
f  <-  sprintf( 'W ~ tx + %s + %s', Z,  paste(sprintf('%s', adjust.for), collapse="+") )
diagnostic_1  <- glm( as.formula(f) , family = poisson(link='log'), data = A.temp)
summary(diagnostic_1)
f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + const(%s) +  %s',Z,  paste(sprintf('const(%s)', adjust.for), collapse="+") )
diagnostic_2  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
summary(diagnostic_2)


                     # estimate   low_ci high_ci y_axis              outcome
# death                  0.0861  0.06950  0.1030      1                death
# fall                   0.0457  0.03120  0.0602      2                 fall
# other_injury           0.0958  0.06170  0.1300      3         other_injury
# diverticular_disease   0.0273  0.01120  0.0434      4 diverticular_disease
# hernia                 0.0213  0.00691  0.0357      5               hernia
# arthropathy            0.0608  0.03490  0.0867      6          arthropathy
# GU_sx                  0.0193 -0.00834  0.0469      7                GU_sx
# hpb                    0.0359  0.01810  0.0537      8                  hpb
# optho                  0.0943  0.05760  0.1310      9                optho



   





# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -5.1418  -0.8396  -0.0241   0.6782   2.9201  

# Coefficients:
#                                     Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                        0.2065572  0.0792065   2.608 0.009112 ** 
# txsbrt                             0.0168701  0.0193091   0.874 0.382288    
# optho_any_date_count               0.0076995  0.0003475  22.158  < 2e-16 ***
# age                                0.0040889  0.0009730   4.202 2.64e-05 ***
# sexMale                           -0.0092509  0.0142637  -0.649 0.516622    
# raceOther                         -0.1030548  0.0432367  -2.384 0.017149 *  
# raceWhite                          0.1044707  0.0277294   3.768 0.000165 ***
# marital.statusNever married       -0.0669835  0.0235593  -2.843 0.004466 ** 
# marital.statusOther                0.0161391  0.0151780   1.063 0.287635    
# marital.statusUnknown             -0.0532879  0.0323421  -1.648 0.099428 .  
# histologyAdenocarcinoma, lepidic   0.2885491  0.0427576   6.748 1.49e-11 ***
# histologyAdenocarcinoma, mixed     0.1690542  0.0386639   4.372 1.23e-05 ***
# histologyAdenocarcinoma, mucinous  0.1486808  0.0461694   3.220 0.001280 ** 
# histologyAdenocarcinoma, NOS       0.1033787  0.0295983   3.493 0.000478 ***
# histologyNon-small cell carcinoma  0.1891434  0.0485000   3.900 9.62e-05 ***
# histologyOther                     0.0921840  0.0331160   2.784 0.005375 ** 
# histologySquamous cell carcinoma   0.0879934  0.0312514   2.816 0.004868 ** 
# DMTRUE                             0.0918648  0.0153981   5.966 2.43e-09 ***
# DMcxTRUE                          -0.0174762  0.0211996  -0.824 0.409732    
# LiverMildTRUE                      0.1432609  0.0159646   8.974  < 2e-16 ***
# PulmonaryTRUE                      0.2434469  0.0170166  14.306  < 2e-16 ***
# PVDTRUE                            0.0644470  0.0149919   4.299 1.72e-05 ***
# CHFTRUE                            0.0467230  0.0167765   2.785 0.005352 ** 
# MITRUE                            -0.0122037  0.0190446  -0.641 0.521655    
# RenalTRUE                          0.0456425  0.0172732   2.642 0.008232 ** 
# StrokeTRUE                         0.0326358  0.0150500   2.168 0.030121 *  
# PUDTRUE                            0.0457131  0.0259265   1.763 0.077870 .  
# RheumaticTRUE                      0.0882238  0.0212130   4.159 3.20e-05 ***
# DementiaTRUE                      -0.0237293  0.0354298  -0.670 0.503012    
# LiverSevereTRUE                    0.0589387  0.0565313   1.043 0.297140    
# ParalysisTRUE                      0.0197639  0.0437057   0.452 0.651121    
# SmokingTRUE                       -0.0105334  0.0157698  -0.668 0.504167    
# o2TRUE                            -0.0567377  0.0270778  -2.095 0.036139 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# (Dispersion parameter for poisson family taken to be 1)

#     Null deviance: 11735  on 8154  degrees of freedom
# Residual deviance: 10068  on 8122  degrees of freedom
# AIC: 31418

# Number of Fisher Scoring iterations: 5


                                         # P-val lower2.5% upper97.5%
# const(tx)sbrt                                0  0.058100   0.129000
# const(optho_any_date_count)                  0  0.001190   0.002410
# const(age)                                   0 -0.001490   0.000956
# const(sex)Male                               0 -0.012200   0.020800
# const(race)Other                             0 -0.004620   0.072200
# const(race)White                             0  0.049000   0.103000
# const(marital.status)Never married           0  0.000248   0.054000
# const(marital.status)Other                   0 -0.001790   0.034200
# const(marital.status)Unknown                 0 -0.023600   0.049000
# const(histology)Adenocarcinoma, lepidic      0 -0.040600   0.045700
# const(histology)Adenocarcinoma, mixed        0 -0.028400   0.055400
# const(histology)Adenocarcinoma, mucinous     0 -0.044000   0.054800
# const(histology)Adenocarcinoma, NOS          0 -0.020000   0.039900
# const(histology)Non-small cell carcinoma     0 -0.054800   0.088200
# const(histology)Other                        0 -0.040400   0.025900
# const(histology)Squamous cell carcinoma      0 -0.043400   0.021600
# const(DM)TRUE                                0  0.018300   0.059100
# const(DMcx)TRUE                              0 -0.004880   0.065700
# const(LiverMild)TRUE                         0  0.011900   0.060100
# const(Pulmonary)TRUE                         0  0.077900   0.113000
# const(PVD)TRUE                               0  0.018400   0.060400
# const(CHF)TRUE                               0  0.033100   0.086900
# const(MI)TRUE                                0 -0.005200   0.054000
# const(Renal)TRUE                             0  0.022700   0.077900
# const(Stroke)TRUE                            0  0.023400   0.067600
# const(PUD)TRUE                               0 -0.025800   0.057600
# const(Rheumatic)TRUE                         0  0.035200   0.109000
# const(Dementia)TRUE                          0  0.028500   0.177000
# const(LiverSevere)TRUE                       0 -0.017300   0.203000
# const(Paralysis)TRUE                         0 -0.019500   0.142000
# const(Smoking)TRUE                           0  0.003150   0.039300
# const(o2)TRUE                                0  0.002680   0.108000



################################
# fofo 
################################

outcome.name
A.temp  <-  A.final2 %>% mutate( 
                      outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                      outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                      W = rowSums( !is.na( across( all_of(outcome.names.temp2))))
)
f  <-  sprintf( 'W ~ tx +optho2_post_count + gout_post_count + %s',  paste(sprintf('const(%s)', adjust.for), collapse="+") )
f  <-  sprintf( 'W ~ tx+optho2_post_count + gout_post_count' )
stage1  <- glm( as.formula(f) , family = poisson(link='log'), data = A.temp)
summary(stage1)


# R> odds.ratios.adj
#                      estimate   low_ci high_ci y_axis              outcome
# death                  0.0861  0.06950  0.1030      1                death
# fall                   0.0457  0.03120  0.0602      2                 fall
# other_injury           0.0958  0.06170  0.1300      3         other_injury
# diverticular_disease   0.0273  0.01120  0.0434      4 diverticular_disease
# hernia                 0.0213  0.00691  0.0357      5               hernia
# arthropathy            0.0608  0.03490  0.0867      6          arthropathy
# GU_sx                  0.0193 -0.00834  0.0469      7                GU_sx
# hpb                    0.0359  0.01810  0.0537      8                  hpb
# optho                  0.0943  0.05760  0.1310      9                optho

#################################
## two step sandbox, with W compiled from multiple 
## Bootstrap edition!
## No parallel
#################################
#set.seed(3)
#adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities)
#Z  <-  'optho'
#outcome.names.temp  <-  setdiff(outcome.names, Z)
#odds.ratios.adj  <-  make.odds.ratio.df ( outcome.names.temp) 
#outcome.i = 1 
#B  <-  1000
#for (outcome.i in 1:length(outcome.names.temp)){ 
#    outcome.name  <-  outcome.names.temp[outcome.i]
#    #outcome.names.temp2  <-  setdiff(outcome.names.temp, c(outcome.name)
#    outcome.names.temp2  <-  setdiff(outcome.names.temp, c(outcome.name, 'fall', 'other_injury'))
#    print(outcome.name)
#    boot_samples  <-  foreach (Bi = 1:B , .combine = 'rbind', .packages=c('dplyr', 'timereg')) %do% {
#        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
#        A.temp  <-  A.final2 %>% mutate( 
#                              outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
#                              outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
#                              W = rowSums( !is.na( across( all_of(outcome.names.temp2))))
#        )
#        f  <-  sprintf( 'W ~ tx + nna(%s) + %s', Z,  paste(sprintf('const(%s)', adjust.for), collapse="+") )
#        table( A.temp$W, useNA="ifany")
#        stage1  <- glm( as.formula(f) , family = poisson(link='log'), data = A.temp)
#        negative_outcome_pred  <-  predict(stage1, type = 'link')
#        A.temp  <-  A.final2 %>% mutate( 
#                              outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
#                              outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F))
#        f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + negative_outcome_pred +  %s',  paste(sprintf('const(%s)', adjust.for), collapse="+") )
#        m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
#        boot.res  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')])
#        return(boot.res)
#    }
#   odds.ratios.adj[outcome.i,1:3]  <-  c( mean(boot_samples[,1]), quantile(boot_samples[,1], c(0.025, 0.975)) )
#   print(odds.ratios.adj[outcome.i,1:3])
#}
#stopCluster(cl)
#odds.ratios.adj
#g2  <-  make.HD.plot(odds.ratios.adj, label_list2)+ ggtitle(sprintf('C) Treatment effect of SBRT (proximal)\n Z: %s', label_list2[Z]))
#ggsave(g2 , width=7, height=2, filename =sprintf('figs/additive.prox4.%s.pdf', subset.name))

##################################
### two step sandbox, with W compiled from multiple 
###No bootstrap
##################################
##adjust.for  <-  c('age', 'sex', 'race', 'marital.status', 'histology', comorbidities)
##Z  <-  'optho'
##outcome.names.temp  <-  setdiff(outcome.names, Z)
##odds.ratios.adj  <-  make.odds.ratio.df ( outcome.names.temp) 
##outcome.i = 1 
##for (outcome.i in 1:length(outcome.names.temp)){ 
##    outcome.name  <-  outcome.names.temp[outcome.i]
##    outcome.names.temp2  <-  setdiff(outcome.names.temp, outcome.name)
##    print(outcome.name)
##    print(outcome.names.temp2)
##    #A  <- A.final %>% mutate(W = rowSums( !is.na( across( all_of(outcome.names.temp)))))
##    A.temp  <-  A.final %>% mutate( 
##                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
##                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
##                          W = rowSums( !is.na( across( all_of(outcome.names.temp2))))
##    )
##    f  <-  sprintf( 'W ~ tx + nna(%s) + %s', Z,  paste(sprintf('const(%s)', adjust.for), collapse="+") )
##    table( A.temp$W, useNA="ifany")
##    stage1  <- glm( as.formula(f) , family = poisson(link='log'), data = A.temp)
##    negative_outcome_pred  <-  predict(stage1, type = 'link')
##    A.temp  <-  A.final %>% mutate( 
##                          outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
##                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F))
##    f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + negative_outcome_pred +  %s',  paste(sprintf('const(%s)', adjust.for), collapse="+") )
##    m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
##    #fofo  <-  predict.timereg(m, resample.iid=T, newdata = A.temp)
##    odds.ratios.adj[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
##    print(odds.ratios.adj[outcome.i,1:3])
##}
##odds.ratios.adj
##g2  <-  make.HD.plot(odds.ratios.adj, label_list2)+ ggtitle(sprintf('C) Treatment effect of SBRT (proximal)\n Z: %s', label_list2[Z]))
##ggsave(g2 , width=7, height=2, filename =sprintf('figs/additive.prox2.%s.pdf', subset.name))

