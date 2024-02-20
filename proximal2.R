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
 subset.name <- 'all.gte.65'
 # subset.name <- 'sens1'

filename.in  <-  sprintf('data/A.final22.%s.RDS', subset.name)
A.final  <-  readRDS(filename.in)  %>% 
    mutate(treatment.year = year(tx.date),
            death.90.day = if_else ( ninety.day.mortality, death, as.Date(NA_character_)),
            death.cause.specific = if_else ( cause.specific.mortality == 'Death', death, as.Date(NA_character_)),
            death.other.cause = if_else ( other.cause.mortality == 'Death', death, as.Date(NA_character_)),
            death.other.cause.gt90day = if_else ( other.cause.mortality == 'Death' & tt > 90, death, as.Date(NA_character_)),
            pre.tx.days = pre.tx.months * 30.5,
    )
A.final$tx  <-  droplevels(A.final$tx)
table( A.final$tx, useNA="ifany")

table( A.final$treatment.year, useNA="ifany")
# Define the analysis
 analysis.name  <-  'primary'
# analysis.name  <-  'sens1'
 # A.final  <-  A.final %>% 
  # filter(age >= 65 & age <= 80) 
# A.final %>% filter (tx == 'sublobar') %>% group_by(tnm.n >0) %>% summarise(n = n()) %>% mutate( freq = n / sum(n) )


# Define X
 X.factor  <-  c('sex', 'race', 'marital.status', 'histology' ) 
X.numeric  <-  c('age', 'size', 'treatment.year', sprintf('%s_pre_count', c(
                     # Diagnoses
                        'smoking', 'o2',  'pneumonia_and_influenza','asthma', 'COPD','interstitial_lung', 'pressure_ulcer', 'ischemic_heart_disease', 'CHF', 'PVD', 'CVD', 'dementia',   'MILDLD','MSLD', 'DIAB_UC', 'DIAB_C',  'RD', 'mental_disorders', 'nervous_system',   'dialysis', 'echo',
                     # Drugs
                     'Insulin', 'Anticoags')
                         )) 
X.factor.min  <-  c('t_stage_8','sex', 'race',  'histology' ) 
adjust.for.scaled  <-  ( c(sprintf('%s_z', X.numeric), X.factor) )  # To use full set, change to X.factor and X.numeric


negative.exposures  <- c('O2accessories', 'walking_aids' , 'hospital_beds_and_supplies' , 'wheelchairs_accessories' , 'transportation_services', 'other_supplies', 'diabetic_footwear' )
ne.any.count.names = sprintf( '%s_any_count', negative.exposures )
ne.post.count.names = sprintf( '%s_post_count', negative.exposures )
ne.pre.count.names = sprintf( '%s_pre_count', negative.exposures )
negative.outcomes.oi  <-  c( 'fall',  'other_injury', 'diverticular_disease', 'hernia',  'arthropathy','GU_sx', 'optho2' )
Zs  <-   sprintf( '%s_z', ne.pre.count.names )
outcome.names  <-  c( 'death', 'death.cause.specific', 'death.90.day', 'death.other.cause' , 'death.other.cause.gt90day')
noc.any.count.names = sprintf( '%s_any_count', negative.outcomes.oi )
noc.post.count.names = sprintf( '%s_post_count', negative.outcomes.oi )
noc.pre.count.names = sprintf( '%s_pre_count', negative.outcomes.oi )
noc.count.names  <-  noc.pre.count.names

B  <-  1000
print('X')
for (i in 1:length(adjust.for.scaled)) cat(i, adjust.for.scaled[i], '\n')
print('Z')
for (i in 1:length(Zs)) cat(i, Zs[i], '\n')
print('W')
for (i in 1:length(noc.count.names)) cat(i, noc.count.names[i], '\n')
cat(i+1,'other.cause.mortality\n')


A.final <- A.final %>% mutate( time.offset = pre.tx.days)
A.final  <- A.final %>% mutate( 
    across( all_of(c(X.numeric,ne.pre.count.names)), scale_ , .names = "{.col}_z" )
)

################################
# Table 1 
################################
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")
names(label_list)[names(label_list) %in% c(X.numeric, X.factor)]
library(arsenal)
tblcontrol <- tableby.control(numeric.stats = c('Nmiss', 'meansd'), numeric.simplify = T, cat.simplify =T, digits = 1,total = T,test = F)
f  <-  sprintf( 'tx ~ %s', paste( c(X.numeric, X.factor), collapse = "+"))
labels(A.final)  <-  label_list
tt <- tableby(as.formula(f), data=A.final, control = tblcontrol)
summary(tt) %>% write2html(sprintf('/Users/george/Research_Local/SEER-Medicare/tbls/table1_2_%s.htm', analysis.name))


##################################
### Inspect negative outcomes 
##################################
#for (outcome.i in 1:length(noc.count.names)){ 
#    outcome.name  <-  noc.count.names[outcome.i]
#    outcome.name  <-  'optho_pre_count'
#    print(outcome.name)
#    A.temp  <-  A.final %>% mutate( outcome.count  = !!rlang::sym(outcome.name))
#    print(table( A.temp$outcome.count, useNA="ifany"))
#    table(A.temp$outcome.count)
#    hist(A.temp %>% filter (outcome.count > 0) %>% pull(outcome.count), breaks = 100)
#}
# g  <-  A.final %>% select( all_of(noc.count.names)) %>% gather() %>% filter(value >0) %>%  ggplot(aes(value)) + geom_histogram() + facet_wrap(~key, scales = 'free')  + scale_y_sqrt() + cowplot::theme_cowplot()  # ggsave(g, file='figs/noc.counts.pdf', width = 10, height = 10, dpi = 300)

################################
# Proximal  function 
################################
two.step  <-  function(A.final2, Zs,  outcome.name,noc.names.temp, adjust.for.scaled, Y.count = F, verbose = F){
    # Stage 1: W1 is non-cause mortality, W2 is the sum of negative control outcomes (excluding the negative outcome of interest)
    W1  <-  'death.other.cause.gt90day'
    A.temp  <-  A.final2 %>% mutate( 
                                    W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt)/365,
                                    W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
                                    W2 = rowSums( ( across( all_of(noc.names.temp)))),
    )
    # W1
    f  <-  sprintf( 'Surv(W1.time, W1.bool) ~ const(tx) +%s+ %s',  
                   paste(sprintf('const(%s)', Zs), collapse="+"), 
                   # paste(sprintf('const(tx:%s)', Zs), collapse="+"), 
                   paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
    m1  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0, silent = 0 )
    if (verbose) 
        print(summary(m1))
    mm  <-  model.matrix(as.formula(f), A.temp)[,-1] 
    coefs  <-  as.matrix(coef( m1)[,'Coef.'])
    negative_outcome_pred1  <-  mm %*% coefs
    # W2
     f  <-  sprintf( 'W2 ~ tx + %s + %s + offset(log(time.offset))', 
     # f  <-  sprintf( 'W2 ~ tx + %s + %s', 
                   paste(sprintf('%s', Zs), collapse="+"),  
                   paste(sprintf('%s', adjust.for.scaled), collapse="+") )
    m2  <- glm( as.formula(f) , family = quasipoisson(link='log'), data = A.temp)
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
         f  <-  sprintf( 'Y.count ~ tx + offset(log(time.offset))  + %s+ What1 + What2',  
                       paste(sprintf('%s', adjust.for.scaled), collapse="+") )
         m  <-  glm( f  ,  data = A.temp, family = quasipoisson(link=log))
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
    mout  <-  two.step(A.final, Zs,  outcome.name, noc.count.names, adjust.for.scaled)
    est_boot <- parallel::mclapply(1:B, function(bb){
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        mout  <-  two.step(A.final2, Zs,  outcome.name, noc.count.names, adjust.for.scaled)
        return(mout)
    }, mc.cores =8)
    se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
   hazard.differences.outcomes.proximal[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
   print(hazard.differences.outcomes.proximal[outcome.i,1:3])
} 


################################
#  Counts
################################

odds.ratios.nocs  <-  make.odds.ratio.df ( noc.count.names) 
odds.ratios.nocs.adj  <-  make.odds.ratio.df ( noc.count.names) 
odds.ratios.nocs.proximal  <-  make.odds.ratio.df ( noc.count.names) 
for (outcome.i in 1:length(noc.count.names)){ 
    outcome.name  <-  noc.count.names[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( outcome.count  = !!rlang::sym(outcome.name))
    # Raw
    f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.offset))' )
    m  <-  glm( f  ,  data = A.temp, family = quasipoisson(link=log))
    odds.ratios.nocs[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'] , confint(m)['txsbrt',]))
    # Adjust for X
    f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.offset))  + %s',  paste(sprintf('%s', adjust.for.scaled), collapse="+") )
     m  <-  glm( f  ,  data = A.temp, family = quasipoisson(link=log))
    odds.ratios.nocs.adj[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'] , confint(m)['txsbrt',]))
    # Proximal
    noc.names.temp  <- setdiff( noc.count.names, c(outcome.name))
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


A.temp %>% count(tx)
################################
# Individual plots
################################
saveRDS(hazard.differences.outcomes, sprintf('data/%s.hazard.differences.outcomes.rds', analysis.name))
saveRDS(hazard.differences.outcomes.adj, sprintf('data/%s.hazard.differences.outcomes.adj.rds', analysis.name))
saveRDS(hazard.differences.outcomes.proximal, sprintf('data/%s.hazard.differences.outcomes.proximal.rds', analysis.name))

height  <-  2
Y.toplot  <-  c('death', 'death.90.day', 'death.other.cause.gt90day', 'death.cause.specific')
hazard.differences.outcomes.toplot  <-  hazard.differences.outcomes[Y.toplot,]
hazard.differences.outcomes.toplot$y_axis  <-  1:nrow(hazard.differences.outcomes.toplot)
g1.a  <-  make.HD.plot(hazard.differences.outcomes.toplot, label_list2) 
g1.b  <-  make.OR.plot(odds.ratios.nocs, label_list2)
g1  <-  g1.a / g1.b+ plot_layout(heights = (c(1,height))) + plot_annotation(title="Raw (Unadjusted)")
ggsave(g1, width=6, height=2.65, filename = sprintf('figs/%s.raw.pdf', analysis.name))
hazard.differences.outcomes.adj.toplot  <-  hazard.differences.outcomes.adj[Y.toplot,]
hazard.differences.outcomes.adj.toplot$y_axis  <-  1:nrow(hazard.differences.outcomes.adj.toplot)
g2.a  <-  make.HD.plot(hazard.differences.outcomes.adj.toplot, label_list2)
g2.b  <-  make.OR.plot(odds.ratios.nocs.adj, label_list2)
 g2  <-  g2.a / g2.b+ plot_layout(heights = (c(1,height)))+ plot_annotation(title="Adjusted")
 ggsave(g2, width=6, height=2.65, filename = sprintf('figs/%s.adj.pdf', analysis.name))
hazard.differences.outcomes.proximal.toplot  <-  hazard.differences.outcomes.proximal[c('death.cause.specific'),]
hazard.differences.outcomes.proximal.toplot$y_axis  <-  1:nrow(hazard.differences.outcomes.proximal.toplot)
g3.a  <-  make.HD.plot(hazard.differences.outcomes.proximal.toplot, label_list2)
g3.b  <-  make.OR.plot(odds.ratios.nocs.proximal, label_list2)
 g3  <-  g3.a / g3.b + plot_layout(heights = (c(0.25,height)))+ plot_annotation(title="Proximal")
ggsave(g3, width=6, height=2.65, filename = sprintf('figs/%s.proximal.pdf', analysis.name))



################################
# Abstract plot 
################################

G  <-  (g2.a  / g2.b+ plot_layout(heights = (c(1,2)))) / (g3.a / g3.b+ plot_layout(heights = (c(1,2))))
G  <- (g2) / (g3) + plot_layout(ncol=1, nrow=4, heights = c(0.25,2, 0.25,2))
G  <-   (( g2.a / g2.b) + plot_annotation('(A) Standard multivariable models')) / (( g3.a  / g3.b)  + plot_annotation('(B) Proximal models')) + plot_layout(ncol=1, nrow=4, heights = c(1,2, 0.25,2))

G  <-    (g2.a+ ggtitle('(A) Standard multivariable models')) / g2.b  /plot_spacer()/  (g3.a + ggtitle('(B) Proximal causal inference models')) / g3.b   + plot_layout(ncol=1, nrow=5, heights = c(1,2, 0.1, 0.25,2))
# G  <-   (( g2.a / g2.b) + plot_annotation('(A) Standard multivariable models')) / (( g3.a  / g3.b)  + plot_annotation('(B) Proximal models')) + plot_layout(ncol=1, nrow=5, heights = c(1,2,0.1, 0.25,2))
# G  <-   (( g2.a / g2.b) + plot_annotation('(A) Standard multivariable models')+ plot_layout(ncol=1, nrow=2, heights = c(1,2)) ) | (( g3.a  / g3.b)  + plot_annotation('(B) Proximal models') + plot_layout(ncol=1, nrow=2, heights=c(0.25,2)))

# + plot_layout(ncol=1, nrow=4)
                                                                                                                 # , heights = c(1,2,0.1, 0.25,2))
ggsave(G, width=7, height=5, filename = ('figs/grant.png'))


################################
# Signle plot
################################

# # height  <-  2
# # Y.toplot  <-  c('death', 'death.90.day', 'death.other.cause.gt90day', 'death.cause.specific')
# # hazard.differences.outcomes.toplot  <-  hazard.differences.outcomes[Y.toplot,]
# # hazard.differences.outcomes.toplot$y_axis  <-  1:nrow(hazard.differences.outcomes.toplot)
# # g1.a  <-  make.HD.plot(hazard.differences.outcomes.toplot, label_list2) + ggtitle('Raw (unadjusted)')
# # g1.b  <-  make.OR.plot(odds.ratios.nocs, label_list2)
# # g1  <-  g1.a / g1.b+ plot_layout(heights = (c(1,height))) #+ plot_annotation(title="Raw")
# #  ggsave(g1, width=5, height=2.55, filename = sprintf('figs/%s.raw.pdf', analysis.name))
# hazard.differences.outcomes.adj.toplot  <-  hazard.differences.outcomes.adj[Y.toplot,]
# hazard.differences.outcomes.adj.toplot$y_axis  <-  1:nrow(hazard.differences.outcomes.adj.toplot)
#  g2.a  <-  make.HD.plot(hazard.differences.outcomes.adj.toplot, label_list2) + ggtitle('Adjusted')
#  g2.b  <-  make.OR.plot(odds.ratios.nocs.adj, label_list2)
# g2  <-  g2.a / g2.b+ plot_layout(heights = (c(1,height)))#+ plot_annotation(title="Adj")
# hazard.differences.outcomes.proximal.toplot  <-  hazard.differences.outcomes.proximal[c('death.cause.specific'),]
# hazard.differences.outcomes.proximal.toplot$y_axis  <-  1:nrow(hazard.differences.outcomes.proximal.toplot)
#  g3.a  <-  make.HD.plot(hazard.differences.outcomes.proximal.toplot, label_list2) + ggtitle('Proximal')
#  g3.b  <-  make.OR.plot(odds.ratios.nocs.proximal, label_list2)
#   g3  <-  g3.a / g3.b + plot_layout(heights = (c(1,height)))#+ plot_annotation(title="Proximal")
# # G  <-  (g1.a/g1.b+ plot_layout(heights = (c(1,2))))| (g2.a / g2.b+ plot_layout(heights = (c(1,2)))) | (g3.a / g3.b+ plot_layout(heights = (c(1,2))))  
# G  <-  (g1.a/g1.b+ plot_layout(heights = (c(1,2))))| (g2.a / g2.b+ plot_layout(heights = (c(1,2)))) | (g3.a / g3.b+ plot_layout(heights = (c(1,2))))  
#  # ggsave(G, width=8.5, height=3, filename = ('figs/grant.pdf'))



################################
# Check effect of F 
################################

# two.step.resid  <-  function(A.final2, Zs,  outcome.name,noc.names.temp, adjust.for.scaled, Y.count = F, verbose = F, resids = NULL){
#     # Stage 1: W1 is non-cause mortality, W2 is the sum of negative control outcomes (excluding the negative outcome of interest)
#     W1  <-  'death.other.cause'
#     A.temp  <-  A.final2 %>% mutate( 
#                                     W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt)/365,
#                                     W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
#                                     W2 = rowSums( ( across( all_of(noc.names.temp)))),
#     )
#     # W1
#     f  <-  sprintf( 'Surv(W1.time, W1.bool) ~ const(tx) +%s+ %s',  
#                    paste(sprintf('const(%s)', Zs), collapse="+"), 
#                    # paste(sprintf('const(tx:%s)', Zs), collapse="+"), 
#                    paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
#     m1  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0, silent = 0 )
#     if (verbose) 
#         print(summary(m1))
#     mm  <-  model.matrix(as.formula(f), A.temp)[,-1] 
#     coefs  <-  as.matrix(coef( m1)[,'Coef.'])
#     negative_outcome_pred1  <-  mm %*% coefs
#     # W2
#     f  <-  sprintf( 'W2 ~ tx + %s + %s + offset(log(time.offset))', 
#                    # f  <-  sprintf( 'W2 ~ tx + %s + %s', 
#                    paste(sprintf('%s', Zs), collapse="+"),  
#                    paste(sprintf('%s', adjust.for.scaled), collapse="+") )
#     m2  <- glm( as.formula(f) , family = quasipoisson(link='log'), data = A.temp)
#     if (verbose) 
#         print(summary(m2))
#     negative_outcome_pred2  <-  predict(m2, type = 'link')
#     A.temp  <-  A.temp %>% mutate( 
#                                   What1 = scale(negative_outcome_pred1[,1]),
#                                   What2 = scale(negative_outcome_pred2)
#     )
#     if (!is.null(resids)) {
#         A.temp$resids  <-  resids
#     }
#     # Stage 2
#     # Poisson for count Ys
#     if (Y.count) {
#         A.temp  <- A.temp %>% mutate( Y.count  = !!rlang::sym(outcome.name),)
#         if (is.null(resids)) {
#             f  <-  sprintf( 'Y.count ~ tx + offset(log(time.offset))  + %s+ What1 + What2',  paste(sprintf('%s', adjust.for.scaled), collapse="+") )
#         } else {
#             f  <-  sprintf( 'Y.count ~ tx + offset(log(time.offset))  + %s+ What1 + What2 + resids',  paste(sprintf('%s', adjust.for.scaled), collapse="+") )
#         }
#         m  <-  glm( f  ,  data = A.temp, family = quasipoisson(link=log))
#         if (verbose) 
#             print(summary(m))
#         if (!is.null(resids)) boot.res  <-  exp(c( coef(m)['resids'] ))
#     }else{
#         # Aalen's for time-to-event Ys
#         A.temp  <- A.temp %>% mutate(
#                                      Y.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
#                                      Y.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
#         )
#         if (is.null(resids)) {
#             f  <-  sprintf( 'Surv(Y.time, Y.bool) ~ const(tx) + const(What1) + const(What2) +  %s',  paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
#         }else {
#             f  <-  sprintf( 'Surv(Y.time, Y.bool) ~ const(tx) + const(What1) + const(What2) +  %s + const(resids)',  paste(sprintf('const(%s)', adjust.for.scaled), collapse="+") )
#         }
#         m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0, silent = 0)
#         if (verbose) 
#             print(summary(m))
#         if (!is.null(resids)) boot.res  <-  coef(m)['const(resids)', 'Coef.']
#     }
#     if (is.null(resids)) {
#         return(resid(m))
#     }else {
#         return(boot.res)
#     }
# }



# outcome.name  <- 'fall_pre_count'
# noc.names.temp  <- setdiff( noc.count.names, c(outcome.name))
# rout  <-  two.step.resid(A.final, Zs,  outcome.name, noc.names.temp, adjust.for.scaled,Y.count =T, verbose=T, resids = NULL)
# outcome.name  <- 'death.cause.specific'
# noc.names.temp  <- setdiff( noc.count.names, c(outcome.name))
# mean.out  <-  two.step.resid(A.final, Zs,  outcome.name, noc.count.names, adjust.for.scaled, resids = rout, verbose=T)
# mean.out

# est_boot <- parallel::mclapply(1:B, function(bb){
#     A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
#     outcome.name  <- 'fall_pre_count'
#     noc.names.temp  <- setdiff( noc.count.names, c(outcome.name))
#     rout  <-  two.step.resid(A.final2, Zs,  outcome.name, noc.names.temp, adjust.for.scaled,Y.count =T, verbose=F, resids = NULL)
#     outcome.name  <- 'death.cause.specific'
#     noc.names.temp  <- setdiff( noc.count.names, c(outcome.name))
#     mout  <-  two.step.resid(A.final2, Zs,  outcome.name, noc.count.names, adjust.for.scaled, resids = rout, verbose=F)
#     return(mout)
# }, mc.cores =8)
# se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
# c( mean.out, mean.out - 1.96*se, mean.out + 1.96*se )








#################################
## Figure out confidence intervals 
#################################
#A.final %>% select(all_of(noc.pre.count.names)) %>% summarise_all(list(~sum(. > 0, na.rm = T))) %>% t() %>% as.data.frame() %>% arrange(desc(V1))



# outcome.i  <-  1
# outcome.name  <-  noc.count.names[outcome.i]
# print(outcome.name)
# A.temp  <-  A.final %>% mutate( outcome.count  = !!rlang::sym(outcome.name))

# f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.offset))  + %s',  paste(sprintf('%s', adjust.for.scaled), collapse="+") )
# m  <-  glm( f  ,  data = A.temp, family = poisson(link=log))
# AER::dispersiontest(m,trafo=1)
# point.est  <-  coef(m)['txsbrt']
# ci  <-  confint(m)['txsbrt',]
# m  <-  glm( f  ,  data = A.temp, family = quasipoisson(link=log))
# point.est  <-  coef(m)['txsbrt']
# ci2  <-  confint(m)['txsbrt',]
# est_boot <- parallel::mclapply(1:B, function(bb){
#     A.temp2  <-  A.temp[sample(nrow(A.temp),replace=T ),]
#     m  <-  glm( f  ,  data = A.temp2, family = poisson(link=log))
#     return(coef(m)['txsbrt'])
# }, mc.cores =8)
# se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
# ci3 <-  c( point.est - 1.96*se, point.est + 1.96*se )

#m  <-  glm( f  ,  data = A.temp, family = quasipoisson(link=log))
# point.est  <-  coef(m)['txsbrt']
# ci2  <-  confint(m)['txsbrt',]
# est_boot  <- rep(NA, B)
# for (i in 1:B) {
#      w  <- rexp( nrow(A.temp), 1)
#      boot.res  <-  glm( f  ,  data = A.temp, family = quasipoisson(link=log), weights = w)
#      est_boot[i]  <-  coef(boot.res)['txsbrt']
# }
#  # est_boot <- parallel::mclapply(1:B, function(bb){
#  #     w  <- rexp( nrow(A.temp), 1)
#  #     print(w[1])
#  #     boot.res  <-  glm( f  ,  data = A.temp, family = quasipoisson(link=log), weights = w)
#  #     print(coef(boot.res)['txsbrt'])
#  #     return(coef(boot.res)['txsbrt'])
#  # }, mc.cores =1)
# se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
# ci4 <-  c( point.est - 1.96*se, point.est + 1.96*se )

# print(ci)
# print(ci2)
# print(ci3)
# print(ci4)


# ################################
##  figure out outlier
# ################################

#    outcome.name  <-  noc.count.names[outcome.i]
#    print(outcome.name)
#    A.temp  <-  A.final %>% mutate( outcome.count  = !!rlang::sym(outcome.name))
#    %>%
#        filter(outcome.count <= 50)
#    table( A.temp$outcome.count, useNA="ifany")
#    # Raw
#    # Proximal
#    noc.names.temp  <- setdiff( noc.count.names, c(outcome.name))
#    mout  <-  two.step(A.final, Zs,  outcome.name, noc.names.temp, adjust.for.scaled,Y.count =T, verbose=F)
#    est_boot <- parallel::mclapply(1:B, function(bb){
#        A.temp2  <-  A.temp[sample(nrow(A.temp),replace=T ),]
#        boot.res  <-  two.step(A.temp2, Zs,  outcome.name, noc.names.temp,adjust.for.scaled, Y.count=T, verbose=F)
#        return(boot.res)
#    }, mc.cores =8)
#    se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
#    print(c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se ))

