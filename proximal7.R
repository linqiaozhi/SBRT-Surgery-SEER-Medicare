library(dplyr)
library(MASS)
source("stats/nc_ph.R")
library(patchwork)
library(lubridate)
library(foreach)
library(doParallel)
library(timereg)
library(addhazard)
library('ggtext')
library(arsenal)
library(ggplot2)
source('utilities.R')
source('loglin.R')
set.seed(3)


################################
# Load data 
################################
subset.name <- 'all.gte.65'

filename.in  <-  sprintf('data/A.final24.%s.RDS', subset.name)
A.final  <-  readRDS(filename.in)  %>% 
    mutate(treatment.year = year(tx.date),
            death.90.day = if_else ( ninety.day.mortality, death, as.Date(NA_character_)),
            death.cause.specific = if_else ( cause.specific.mortality == 'Death', death, as.Date(NA_character_)),
            death.other.cause = if_else ( other.cause.mortality == 'Death', death, as.Date(NA_character_)),
            death.other.cause.gt90day = if_else ( other.cause.mortality == 'Death' & tt > 90, death, as.Date(NA_character_)),
            pre.tx.days = pre.tx.months * 30.5,
    ) 
A.final$tx  <-  droplevels(A.final$tx)

analysis.name  <-  'nb'
# analysis.name  <-  'sens1'
# A.final  <-  A.final %>% 
 # filter(age >= 65 & age <= 80) 
# A.final %>% filter (tx == 'sublobar') %>% group_by(tnm.n >0) %>% summarise(n = n()) %>% mutate( freq = n / sum(n) )

# Define X
X.factor  <-  c('sex', 'race',  'histology' ) 
X.numeric  <-  c('age', 'size', 'treatment.year') 
# X.count  <-  sprintf('%s_pre_count', c( 'smoking', 'o2',  'pneumonia_and_influenza','asthma', 'COPD','interstitial_lung'))
# Xs  <-  ( c(sprintf('%s_z', X.numeric),sprintf('%s_s', X.count), X.factor) )  # To use full set, change to X.factor and X.numeric
Xs  <-  c(sprintf('%s_z', X.numeric),  X.factor)

Z.count  <- c('O2accessories', 'walking_aids' , 'hospital_beds_and_supplies' , 'wheelchairs_accessories' , 'transportation_services', 'other_supplies', 'diabetic_footwear' ,  'pressure_ulcer', 'ischemic_heart_disease', 'CHF', 'PVD', 'CVD', 'dementia',   'MILDLD','MSLD', 'DIAB_UC', 'DIAB_C',  'RD', 'mental_disorders', 'nervous_system',   'dialysis', 'echo', 'Insulin', 'Anticoags',  'smoking', 'o2',  'pneumonia_and_influenza','asthma', 'COPD','interstitial_lung')
Z.count.unscaled = sprintf( '%s_pre_count', Z.count )
Zs  <-   sprintf( '%s_pre_count_s', Z.count )

W.count  <-  c( 'fall',  'other_injury', 'diverticular_disease', 'hernia',  'arthropathy','GU_sx', 'optho2' )
Ws = sprintf( '%s_pre_count', W.count )

outcome.names  <-  c( 'death', 'death.cause.specific', 'death.90.day', 'death.other.cause' , 'death.other.cause.gt90day')



A.final <- A.final %>% mutate( time.offset = pre.tx.days)

A.final  <- A.final %>% mutate( 
     across( all_of(c(Z.count.unscaled)), function(x) scale_(x/time.offset), .names = "{.col}_s" )
 )
A.final  <- A.final %>% mutate( 
    across( all_of(c(X.numeric)), scale_, .names = "{.col}_z" )
)


B  <-  500
print('X')
for (i in 1:length(Xs)) cat(i, Xs[i], '\n')
print('Z')
for (i in 1:length(Zs)) cat(i, Zs[i], '\n')
print('W')
for (i in 1:length(Ws)) cat(i, Ws[i], '\n')
cat(i+1,'other.cause.mortality\n')


#################################
## Table 1 
#################################
#Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")
#names(label_list)[names(label_list) %in% c(X.numeric, X.factor)]
#library(arsenal)
#tblcontrol <- tableby.control(numeric.stats = c('Nmiss', 'meansd'), numeric.simplify = T, cat.simplify =T, digits = 1,total = T,test = F)
#f  <-  sprintf( 'tx ~ %s', paste( c(X.numeric, X.factor), collapse = "+"))
#labels(A.final)  <-  label_list
#tt <- tableby(as.formula(f), data=A.final, control = tblcontrol)
#summary(tt) %>% write2html(sprintf('/Users/george/Research_Local/SEER-Medicare/tbls/table1_2_%s.htm', analysis.name))

#################################
## Inspect each count negative outcome 
#################################

#set.seed(3)
## Change this from 1 to 7 to see each one
#outcome.i <- 5
#outcome.name  <-  Ws[outcome.i]
#A.temp  <-  A.final %>% mutate( outcome.count  = !!rlang::sym(outcome.name))
#g  <- ggplot(A.temp, aes(x=outcome.count)) + geom_histogram(binwidth=1) + ggtitle(outcome.name)
#ggsave(g, file=sprintf('figs/vars/%s_hist.pdf', outcome.name), width = 5, height = 4, dpi = 300)
#f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.offset))  + %s',  paste(sprintf('%s', Xs), collapse="+") )
#m.nb <-  glm.nb(f, data = A.temp, control = glm.control(maxit=100))
#summary(m.nb) %>% write2html(sprintf('/Users/george/Research_Local/SEER-Medicare/figs/vars/%s_nb_model.htm', outcome.name), quiet=T)
#simulationOutput.nb <- simulateResiduals(m.nb, n=1000)
## pdf to file
#pdf(sprintf('figs/vars/%s_nb_dharma.pdf', outcome.name), width = 10, height = 7)
#plot(simulationOutput.nb)
#dev.off()


################################
# Proximal  function 
################################
 verbose=F
 Y.count.bool=T
A.final2 = A.final
two.step  <-  function(A.final2, Zs,  outcome.name,noc.names.temp, Xs, Y.count.bool = F, verbose = F){
      X.subset  <-1:length(Xs)
    # Stage 1: W1 is non-cause mortality, W2 is the sum of negative control outcomes (excluding the negative outcome of interest)
    W1  <-  'death.other.cause.gt90day'
    A.temp  <-  A.final2 %>% mutate( 
                                    W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt)/365,
                                    W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
                                    W2 = rowSums( ( across( all_of(noc.names.temp)))),
    )
    # %>% filter (CHF_pre_count_s < 4 & MILDLD_pre_count_s < 4 & MSLD_pre_count_s < 4)
    # W1
    f  <-  sprintf( 'Surv(W1.time, W1.bool) ~ const(tx) +%s+ %s',  
                   paste(sprintf('const(%s)', Zs), collapse="+"), 
                   paste(sprintf('const(%s)', Xs[X.subset]), collapse="+") )
    m1  <-  aalen( as.formula(f) ,  data = A.temp, robust = T, silent = 0 )
     if (verbose) {
         print('W1 model')
        print(summary(m1))
     }
    mm  <-  model.matrix(as.formula(f), A.temp)[,-1] 
    coefs  <-  as.matrix(coef( m1)[,'Coef.'])
    negative_outcome_pred1  <-  mm %*% coefs
    # W2
    f  <-  sprintf( 'W2 ~ tx + %s + %s + offset(log(time.offset))', 
                  paste(sprintf('%s', Zs), collapse="+"),  
                  paste(sprintf('%s', Xs), collapse="+") )
    # m2  <-  glm.nb(as.formula(f), data = A.temp, control = glm.control(maxit=500, trace=0), init.theta = 0.1)
    # if (verbose) {
    #      print('W2 model')
    #    print(summary(m2))
    # }
    # negative_outcome_pred2  <-  predict(m2, type = 'link')
    W2.nb_fit.out  <- negbin_fit(A.temp$W2, mm, offset = log(A.temp$time.offset), variance = F)
    # head(cbind(coef(m2),W2.nb_fit.out$ESTIMATE[-1]))
    negative_outcome_pred2  <-  W2.nb_fit.out$PREDICT
    A.temp  <-  A.temp %>% mutate( 
                                  What1 = negative_outcome_pred1[,1],
                                  What2 = negative_outcome_pred2
    )
    # Stage 2
    if (Y.count.bool) {
       A.temp  <- A.temp %>% mutate( Y.count  = !!rlang::sym(outcome.name),)
         f  <-  sprintf( 'Y.count ~ tx + offset(log(time.offset))  + %s+ What1 + What2',  
                       paste(sprintf('%s', Xs[X.subset]), collapse="+") )
        mm  <-  model.matrix(as.formula(f), A.temp) 
        # m  <-  glm( f  ,  data = A.temp, family = quasipoisson(link=log))
          # m  <-  glm.nb(f, data = A.temp, control = glm.control(maxit=500, trace=0), init.theta = 0.2)
         nb_fit.out  <- negbin_fit(A.temp$Y.count, mm, offset = log(A.temp$time.offset), variance = F)
         boot.res  <-  nb_fit.out$ESTIMATE[3]
         # cbind(coef(m),nb_fit.out$ESTIMATE[-1])
        if (verbose) {
            # print('Stage 2 moel')
            # print(summary(m))
            # print(car::vif(m))
        }
    }
    # summary(A.temp$CHF_pre_count_s)
    if(!Y.count.bool){
        # Aalen's for time-to-event Ys
        A.temp  <- A.temp %>% mutate(
                                     Y.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                                     Y.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        )
        f  <-  sprintf( 'Surv(Y.time, Y.bool) ~ const(tx) + const(What1) + const(What2) +  %s',  
                       paste(sprintf('const(%s)', Xs), collapse="+") )
        m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0, silent = 0)
        if (verbose) {
            print(summary(m))
        }
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
                          # outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365, 
                           outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365 + runif(dim(A.final)[1] ,0,1)*1e-8,
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        )
    # Raw
    f  <- as.formula('Surv(outcome.time, outcome.bool) ~ const(tx)')
    m  <-  aalen( f ,  data = A.temp, robust = 0)
    f  <- as.formula('Surv(outcome.time, outcome.bool) ~ tx')
    m2  <-  ah( f ,  data = A.temp, robust = 0, ties=0, verbose=T)
    hazard.differences.outcomes[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
    # Adjust for X
    f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + %s',  paste(sprintf('const(%s)', c(Xs, Zs)), collapse="+") )
    m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
    hazard.differences.outcomes.adj[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
    # Proximal
    mout  <-  two.step(A.final, Zs,  outcome.name, Ws, Xs)
    est_boot <- parallel::mclapply(1:B, function(bb){
        A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
        mout  <-  two.step(A.final2, Zs,  outcome.name, Ws, Xs)
        return(mout)
    }, mc.cores =8)
   hazard.differences.outcomes.proximal[outcome.i,1:3]  <-  c(mout[1], quantile(unlist(est_boot), 0.025), quantile(unlist(est_boot), 0.975))
   print(hazard.differences.outcomes.proximal[outcome.i,1:3])
} 


################################
#  Counts
################################

B
odds.ratios.nocs  <-  make.odds.ratio.df ( Ws) 
odds.ratios.nocs.adj  <-  make.odds.ratio.df ( Ws) 
odds.ratios.nocs.proximal  <-  make.odds.ratio.df ( Ws) 
outcome.i <- 1
for (outcome.i in 1:length(Ws[-7])){ 
    outcome.name  <-  Ws[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( outcome.count  = !!rlang::sym(outcome.name))
    # Raw
    f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.offset))' )
    mm  <- model.matrix(as.formula(f), data=A.temp)
    m  <-  glm.nb(f, data = A.temp, control = glm.control(maxit=100), init.theta = 0.1)
    odds.ratios.nocs[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'] , confint(m)['txsbrt',]))
    # Adjust for X
    f  <-  sprintf( 'outcome.count ~ tx + offset(log(time.offset))  + %s',  paste(sprintf('%s', c(Xs, Zs)), collapse="+") )
    mm  <- model.matrix(as.formula(f), data=A.temp)
    nb_fit.out  <- negbin_fit(A.temp$outcome.count, mm, offset = log(A.temp$time.offset))
    est  <-  nb_fit.out$ESTIMATE[3]
    se  <- nb_fit.out$SE[3]
    # print(sprintf('Outcome name: %s. Adj CI: (%.2f, %.2f)', outcome.name, exp(est-1.96*se), exp(est+1.96*se))) 
    # odds.ratios.nocs.adj[outcome.i,1:3]  <-  exp(c( coef(m)['txsbrt'] , confint(m)['txsbrt',]))
    odds.ratios.nocs.adj[outcome.i,1:3]  <-  exp(c( est, est-1.96*se, est+1.96*se))
    # Proximal
    noc.names.temp  <- setdiff( Ws, c(outcome.name))
    mout  <-  two.step(A.final, Zs,  outcome.name, noc.names.temp, Xs,Y.count.bool =T, verbose=F)
    est_boot <- parallel::mclapply(1:B, function(bb){
                                       A.final2  <-  A.final[sample(nrow(A.final),replace=T ),]
                                       boot.res  <-  two.step(A.final2, Zs,  outcome.name, noc.names.temp,Xs, Y.count.bool=T, verbose=F)
                                       return(boot.res)
    }, mc.cores =8)
    # se  <-  sd(unlist(lapply(est_boot, function(x) x[1])))
    # odds.ratios.nocs.proximal[outcome.i,1:3]  <-  c( mout[1], mout[1] - 1.96*se, mout[1] + 1.96*se )
    odds.ratios.nocs.proximal[outcome.i,1:3]  <- c( exp(mout[1]), quantile(exp(unlist(est_boot)), 0.025), quantile(exp(unlist(est_boot)), 0.975))
    print(odds.ratios.nocs.adj[outcome.i,1:3])
    print(odds.ratios.nocs.proximal[outcome.i,1:3])
    cat('==========\n\n\n')
}

odds.ratios.nocs.proximal

## table( A.temp$arthropathy_pre_count, useNA="ifany")
#A.temp %>% count(tx)
#################################
## Individual plots
#################################
#saveRDS(hazard.differences.outcomes, sprintf('data/%s.hazard.differences.outcomes.rds', analysis.name))
#saveRDS(hazard.differences.outcomes.adj, sprintf('data/%s.hazard.differences.outcomes.adj.rds', analysis.name))
#saveRDS(hazard.differences.outcomes.proximal, sprintf('data/%s.hazard.differences.outcomes.proximal.rds', analysis.name))

height  <-  2
Y.toplot  <-  c('death', 'death.90.day', 'death.other.cause.gt90day', 'death.cause.specific')
hazard.differences.outcomes.toplot  <-  hazard.differences.outcomes[Y.toplot,]
hazard.differences.outcomes.toplot$y_axis  <-  1:nrow(hazard.differences.outcomes.toplot)
g1.a  <-  make.HD.plot(hazard.differences.outcomes.toplot, label_list2) 
g1.b  <-  make.OR.plot(odds.ratios.nocs, label_list2)
g1  <-  g1.a / g1.b+ plot_layout(heights = (c(1,height))) + plot_annotation(title="Raw (Unadjusted)")
hazard.differences.outcomes.adj.toplot  <-  hazard.differences.outcomes.adj[Y.toplot,]
hazard.differences.outcomes.adj.toplot$y_axis  <-  1:nrow(hazard.differences.outcomes.adj.toplot)
g2.a  <-  make.HD.plot(hazard.differences.outcomes.adj.toplot, label_list2)
g2.b  <-  make.OR.plot(odds.ratios.nocs.adj, label_list2)
 g2  <-  g2.a / g2.b+ plot_layout(heights = (c(1,height)))+ plot_annotation(title="Adjusted")
hazard.differences.outcomes.proximal.toplot  <-  hazard.differences.outcomes.proximal[c('death.cause.specific'),]
hazard.differences.outcomes.proximal.toplot$y_axis  <-  1:nrow(hazard.differences.outcomes.proximal.toplot)
g3.a  <-  make.HD.plot(hazard.differences.outcomes.proximal.toplot, label_list2)
g3.b  <-  make.OR.plot(odds.ratios.nocs.proximal, label_list2)
 g3  <-  g3.a / g3.b + plot_layout(heights = (c(0.25,height)))+ plot_annotation(title="Proximal")
G  <-  (g1.a/g1.b+ plot_layout(heights = (c(1,2)))) / (g2.a / g2.b+ plot_layout(heights = (c(1,2)))) / (g3.a / g3.b+ plot_layout(heights = (c(1,2))))  
# ggsave(G, height=8.5, width=5, filename = sprintf('figs/%s.pdf', analysis.name))
ggsave(g1, width=6, height=2.65, filename = sprintf('figs/%s.raw.pdf', analysis.name))
ggsave(g2, width=6, height=2.65, filename = sprintf('figs/%s.adj.pdf', analysis.name))
ggsave(g3, width=6, height=2.65, filename = sprintf('figs/%s.proximal.pdf', analysis.name))

