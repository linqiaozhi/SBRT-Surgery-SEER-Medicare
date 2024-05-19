library(pci2s)
library(dplyr)
library(MASS)
library(patchwork)
library(lubridate)
library(timereg)
library(addhazard)
library(survival)
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

filename.in  <-  sprintf('data/A.final.%s.RDS', subset.name)
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

# preprocessing
table( A.final$treatment.year, A.final$tx, useNA="ifany")
table( A.final$race, useNA="ifany")
table( A.final$histology, A.final$tx, useNA="ifany")
A.final <- A.final %>% mutate( race2 = ifelse (race == 'White' , 'White', 'Other'),
                              treatment.year2 = as.character(treatment.year),
                              treatment.year2 = (ifelse (treatment.year2 %in% c('2019', '2020'), '2019-2020', treatment.year2)),
                                histology2 = ifelse (grepl('Adeno', histology), 'Adenocarcinoma', 'Squamous cell'))
# analysis.name  <-  'sens.with.rad'
# analysis.name  <-  'sens1'
 analysis.name  <-  'nb2020'
# A.final  <-  A.final %>% 
 # filter(age >= 65 & age <= 80) 
# A.final %>% filter (tx == 'sublobar') %>% group_by(tnm.n >0) %>% summarise(n = n()) %>% mutate( freq = n / sum(n) )

# Define X
X.factor  <-  c('sex', 'race',  'histology2', 'treatment.year2' ) 
X.numeric  <-  c('age', 'size') 
Xs  <-  c(sprintf('%s_z', X.numeric),  X.factor)

 Z.count  <- c('O2accessories', 'walking_aids' , 'hospital_beds_and_supplies' , 'wheelchairs_accessories' , 'transportation_services', 'other_supplies', 'diabetic_footwear' ,  'pressure_ulcer', 'ischemic_heart_disease', 'CHF', 'PVD', 'CVD', 'dementia',   'MILDLD','MSLD', 'DIAB_UC', 'DIAB_C',  'RD', 'mental_disorders', 'nervous_system',   'dialysis', 'echo', 'Insulin', 'Anticoags',  'smoking', 'o2',  'pneumonia_and_influenza','asthma', 'COPD','interstitial_lung')
Z.count.unscaled = sprintf( '%s_pre_month_count', Z.count )
Zs  <-   sprintf( '%s_pre_month_count_s', Z.count )

 W.count  <-  c( 'fall',  'other_injury', 'diverticular_disease', 'hernia',  'arthropathy','GU_sx')
Ws = sprintf( '%s_pre_month_count', W.count )

outcome.names  <-  c( 'death', 'death.cause.specific', 'death.90.day',  'death.other.cause.gt90day')


# use cut and quantile to create group the vector
quartile  <- function(x) {
    scaled  <-  as.numeric(x)
    breaks  <- c(0,quantile(scaled[scaled!=0], probs = c(0, 0.25, 0.5, 0.75), na.rm = T, names=F), max(scaled, na.rm = T))
    scaled  <-  cut(x, breaks, include.lowest = T, labels = c('0', '1', '2', '3', '4'), right=F) %>% as.character %>% as.numeric
    return(scaled)
}

A.final <- A.final %>% mutate( time.offset = pre.tx.months)

A.final  <- A.final %>% mutate( 
     across( all_of(c(Z.count.unscaled, Ws)), function(x) (x >0), .names = "{.col}_unbinned" ),
     across( all_of(c(Z.count.unscaled)), function(x) quartile(x/time.offset), .names = "{.col}_s" )
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

table( A.final$histology, useNA="ifany")

#################################
## Table 1 
#################################
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")
#names(label_list)[names(label_list) %in% c(X.numeric, X.factor)]
library(arsenal)
tblcontrol <- tableby.control(numeric.stats = c('Nmiss', 'meansd'), numeric.simplify = T, cat.simplify =T, digits = 1,total = T,test = F)
f  <-  sprintf( 'tx ~ %s', paste( c(X.numeric, X.factor,'histology', gsub('_count_s', '_count_unbinned', Zs), gsub('_count', '_count_unbinned', c(Ws))), collapse = "+"))
labels(A.final)  <-  label_list
tt <- tableby(as.formula(f), data=A.final, control = tblcontrol)
summary(tt) %>% write2html(sprintf('/Users/george/Research_Local/SEER-Medicare/tbls/table1_2_%s.htm', analysis.name))

table( A.final$SEQUENCE_NUMBER, useNA="ifany")
################################
# Proximal  function 
################################
# verbose=F
# Y.count.bool=F
# Y.count.bool=T
# A.final2 = A.final
# noc.names.temp  <-  Ws
two.step  <-  function(A.final2, Zs,  outcome.name,noc.names.temp, Xs, Y.count.bool = F, verbose = F, B=500, skip.W1 = F){
    W1  <-  'death.other.cause.gt90day'
    A.temp  <-  A.final2 %>% mutate( 
                                    W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt)/365,
                                    W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
                                    # W2 = rowSums( ( across( all_of(noc.names.temp)))),
    )
    A_ = (A.temp$tx == 'sbrt')*1.0
    X_  <-  model.matrix(as.formula(sprintf('~ %s', paste(Xs, collapse = '+'))),  A.temp)[,-1]
    W_ <- cbind( A.temp$W1.time, A.temp[,noc.names.temp])
    D2_  <- A.temp$W1.bool*1.0
    Z_  <- A.temp[,Zs]
    N_  <- dim(A.temp)[1]
    offset_  <- log(A.temp$time.offset)
        Xw = append( list(  as.matrix(cbind(  A_, X_, Z_ ))),
                replicate (length(noc.names.temp),   list(as.matrix(cbind(1, A_, X_, Z_ )))))
    if (Y.count.bool) {
        A.temp  <- A.temp %>% mutate( Y.count  = !!rlang::sym(outcome.name),)
        Y_ = A.temp$Y.count
        out.model  <- pcinb2s(Y = Y_, offset= offset_,  A = A_, X = X_,
                              W = W_, Z = Z_, nboot = B,
                              Xw = Xw,        
                              nco_type = c("ah", rep("negbin", length(noc.names.temp))),
                              nco_args = append( list(list(offset = rep(0, N_), event = D2_)),
                                                replicate (length(noc.names.temp),  list(offset = offset_, init = NA), simplify=F)),
                              se_method= 'analytic')
    }
    if(!Y.count.bool){
        A.temp  <- A.temp %>% mutate(
                                     Y.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                                     Y.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        )
        Y_ = A.temp$Y.time
        D_ = A.temp$Y.bool*1.0
        if (skip.W1) {
            print('Not using other cause mortality in the W set')
            Xw =  replicate (length(noc.names.temp),   list(as.matrix(cbind(1, A_, X_, Z_ ))))
            W_ <- cbind(  A.temp[,noc.names.temp])
            out.model  <- pciah2s(Y = Y_, D = D_,  A = A_, X = X_,
                                  W = W_, Z = Z_, nboot = B,
                                  Xw = Xw,        
                                  nco_type =  rep("negbin", length(noc.names.temp)),
                                  nco_args =  replicate (length(noc.names.temp),  list(offset = offset_, init = NA), simplify=F),
                                  se_method= 'analytic')
        }else {
            out.model  <- pciah2s(Y = Y_, D = D_,  A = A_, X = X_,
                                  W = W_, Z = Z_, nboot = B,
                                  Xw = Xw,        
                                  nco_type = c("ah", rep("negbin", length(noc.names.temp))),
                                  nco_args = append( list(list(offset = rep(0, N_), event = D2_)),
                                                    replicate (length(noc.names.temp),  list(offset = offset_, init = NA), simplify=F)),
                                  se_method= 'analytic')
        }
    }
    return(out.model)
}


################################
#  Counts
################################
# options(warn=2)
# options(error = recover)
odds.ratios.nocs  <-  make.odds.ratio.df ( Ws) 
odds.ratios.nocs.adj  <-  make.odds.ratio.df ( Ws) 
odds.ratios.nocs.proximal  <-  make.odds.ratio.df ( Ws) 
for (outcome.i in 1:length(Ws)){ 
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
    odds.ratios.nocs.adj[outcome.i,1:3]  <-  exp(c( est, est-1.96*se, est+1.96*se))
    # Proximal
    noc.names.temp  <- setdiff( Ws, c(outcome.name))
    mout  <-  two.step(A.final, Zs,  outcome.name, noc.names.temp, Xs,Y.count.bool =T, verbose=F)
    est <- mout$ESTIMATE[2]
    se  <- mout$SE[2]
    odds.ratios.nocs.proximal[outcome.i,1:3]  <-  exp(c( est, est- 1.96*se, est + 1.96*se ))
    print(odds.ratios.nocs[outcome.i,1:3])
    print(odds.ratios.nocs.adj[outcome.i,1:3])
    print(odds.ratios.nocs.proximal[outcome.i,1:3])
    cat('==========\n\n\n')
}

################################
# Time-to-event outcomes 
################################
hazard.differences.outcomes  <-  make.odds.ratio.df ( outcome.names) 
hazard.differences.outcomes.adj  <-  make.odds.ratio.df ( outcome.names) 
hazard.differences.outcomes.proximal  <-  make.odds.ratio.df ( outcome.names) 
table( A.final$tt, useNA="ifany")
for (outcome.i in 1:length(outcome.names)){ 
    outcome.name  <-  outcome.names[outcome.i]
    A.temp  <-  A.final %>% mutate( 
                           outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/ 365+ runif(dim(A.final)[1] ,0,1)*1e-8,,
                          outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        ) %>% filter (tnm.n ==0)
    # Raw
    f  <- as.formula('Surv(outcome.time, outcome.bool) ~ const(tx)')
    m  <-  aalen( f ,  data = A.temp, robust = 0)
    hazard.differences.outcomes[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
    # Adjust for X
    f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + %s',  paste(sprintf('const(%s)', c(Xs, Zs)), collapse="+") )
    m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
    hazard.differences.outcomes.adj[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
    # Proximal
    if (outcome.name != 'death.cause.specific') {
        mout  <-  two.step( A.final, Zs,  outcome.name, Ws, Xs, skip.W1=T)
    }else {
        mout  <-  two.step( A.final, Zs,  outcome.name, Ws, Xs, skip.W1=F)
    }
    est <- mout$ESTIMATE[1]
     se  <- mout$SE[1]
   hazard.differences.outcomes.proximal[outcome.i,1:3]  <-  c(est, est - 1.96*se, est + 1.96*se)
   print(hazard.differences.outcomes.proximal[outcome.i,1:3])
} 



## table( A.temp$arthropathy_pre_count, useNA="ifany")
A.final %>% count(tx)
#################################
## Individual plots
#################################
saveRDS(hazard.differences.outcomes, sprintf('data/%s.hazard.differences.outcomes.rds', analysis.name))
saveRDS(hazard.differences.outcomes.adj, sprintf('data/%s.hazard.differences.outcomes.adj.rds', analysis.name))
fofo  <- readRDS(sprintf('data/%s.hazard.differences.outcomes.adj.rds', analysis.name))
saveRDS(hazard.differences.outcomes.proximal, sprintf('data/%s.hazard.differences.outcomes.proximal.rds', analysis.name))
saveRDS(odds.ratios.nocs, sprintf('data/%s.odds.ratios.nocs.rds', analysis.name))
saveRDS(odds.ratios.nocs.adj, sprintf('data/%s.odds.ratios.nocs.adj.rds', analysis.name))
saveRDS(odds.ratios.nocs.proximal, sprintf('data/%s.odds.ratios.nocs.proximal.rds', analysis.name))
fifi  <- readRDS(sprintf('data/%s.hazard.differences.outcomes.proximal.rds', analysis.name))

height  <-  2
Y.toplot  <-  c('death', 'death.90.day', 'death.other.cause.gt90day', 'death.cause.specific')
hazard.differences.outcomes.toplot  <-  hazard.differences.outcomes[Y.toplot,]
hazard.differences.outcomes.toplot$y_axis  <-  1:nrow(hazard.differences.outcomes.toplot)
label_list3  <-  label_list2
names(label_list3)  <- gsub( '_unbinned', '', names(label_list3))
g1.a  <-  make.HD.plot(hazard.differences.outcomes.toplot, label_list3) 
g1.b  <-  make.OR.plot(odds.ratios.nocs, label_list3)
g1  <-  g1.a / g1.b+ plot_layout(heights = (c(1,height))) + plot_annotation(title="Raw (Unadjusted)")
hazard.differences.outcomes.adj.toplot  <-  hazard.differences.outcomes.adj[Y.toplot,]
hazard.differences.outcomes.adj.toplot$y_axis  <-  1:nrow(hazard.differences.outcomes.adj.toplot)
g2.a  <-  make.HD.plot(hazard.differences.outcomes.adj.toplot, label_list3)
g2.b  <-  make.OR.plot(odds.ratios.nocs.adj, label_list3)
 g2  <-  g2.a / g2.b+ plot_layout(heights = (c(1,height)))+ plot_annotation(title="Adjusted")
hazard.differences.outcomes.proximal.toplot  <-  hazard.differences.outcomes.proximal[c('death.cause.specific'),]
hazard.differences.outcomes.proximal.toplot$y_axis  <-  1:nrow(hazard.differences.outcomes.proximal.toplot)
g3.a  <-  make.HD.plot(hazard.differences.outcomes.proximal.toplot, label_list3)
g3.b  <-  make.OR.plot(odds.ratios.nocs.proximal, label_list3)
 g3  <-  g3.a / g3.b + plot_layout(heights = (c(0.25,height)))+ plot_annotation(title="Proximal")
G  <-  (g1.a/g1.b+ plot_layout(heights = (c(1,2)))) / (g2.a / g2.b+ plot_layout(heights = (c(1,2)))) / (g3.a / g3.b+ plot_layout(heights = (c(1,2))))  
# ggsave(G, height=8.5, width=5, filename = sprintf('figs/%s.pdf', analysis.name))
ggsave(g1, width=6, height=2.65, filename = sprintf('figs/%s.raw.pdf', analysis.name))
ggsave(g2, width=6, height=2.65, filename = sprintf('figs/%s.adj.pdf', analysis.name))
ggsave(g3, width=6, height=2.65, filename = sprintf('figs/%s.proximal.pdf', analysis.name))




toprint  <- hazard.differences.outcomes.adj
toprint  <- odds.ratios.nocs.proximal
for (i in 1:nrow(toprint)) {
    cat( sprintf( '%s| %.2f (%.2f, %.2f)\n', label_list3[rownames(toprint)[i]], toprint[i,1], toprint[i,2], toprint[i,3]) ) }

table( A.final$tnm.n, useNA="ifany")

sum(A.final$tx == 'sublobar' & A.final$tnm.n %in% c('1','2','3')) / sum(A.final$tx == 'sublobar')
