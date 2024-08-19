 devtools::load_all('../pci2s_gcl/pci2s')
# devtools::load_all('../pci2s')
 # library(pci2s)
library(glmnet)
library(arsenal)
library(dplyr)
library(tidyr)
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
source('two.step.variable.selection.R')
 source('baseline.hazard.R')
seed  <-  3
alpha  <- 1
set.seed(seed)
get.selected.columns  <-  function(fit, verbose=F) {
    coefs  <-  coef(fit, s = 'lambda.min')[,1] != 0
    selected.columns  <- names(coefs)[c(-1,-2)][coefs[c(-1,-2)]]
    if (verbose ==T) print(sprintf('Excluded %s', paste(names(coefs)[c(-1,-2)][!coefs[c(-1,-2)]], collapse = ', ')))
    return(selected.columns)
}
quartile  <- function(x) {
    scaled  <-  as.numeric(x)
    breaks  <- c(0,quantile(scaled[scaled!=0], probs = c(0, 0.25, 0.5, 0.75), na.rm = T, names=F), max(scaled, na.rm = T))
    scaled  <-  cut(x, breaks, include.lowest = T, labels = c('0', '1', '2', '3', '4'), right=F) %>% as.character %>% as.numeric
    return(scaled)
}
check.if.present  <- function(needles,haystack) sapply(needles, function(x) any(grepl(sprintf('^%s',x), haystack)))


################################
# Load data 
################################
 subset.name <- 'sens1'
 variable.selection  <- T
 analysis.name  <-  sprintf('sens1.lasso.seed%d', seed)

 # subset.name <- 'sens1'
 # variable.selection  <- F
 # analysis.name  <-  sprintf('sens1.nolasso.seed%d', seed)


# subset.name <- 'all.gte.65'
# variable.selection  <- F
# analysis.name  <-  sprintf('nolasso.seed%d', seed)

# subset.name <- 'all.gte.65'
# variable.selection  <- T
# analysis.name  <-  sprintf('overall.seed%d', seed)

filename.in  <-  sprintf('data/A.final3.%s.RDS', subset.name)
A.final  <-  readRDS(filename.in)  %>% 
    mutate(treatment.year = year(tx.date),
            death.90.day = if_else ( ninety.day.mortality, death, as.Date(NA_character_)),
            death.cause.specific = if_else ( cause.specific.mortality == 'Death', death, as.Date(NA_character_)),
            death.other.cause = if_else ( other.cause.mortality == 'Death', death, as.Date(NA_character_)),
            death.other.cause.gt90day = if_else ( other.cause.mortality == 'Death' & tt > 90, death, as.Date(NA_character_)),
            pre.tx.days = pre.tx.months * 30.5,
    )

table( A.final$tx, useNA="ifany")
sum(A.final$tnm.n != 0 & A.final$tx == 'sublobar')/ sum(A.final$tx == 'sublobar')
sum(A.final$tx == 'sublobar')
# A.final %>% filter(tt == 0) %>% t 
#%>% filter (tt > 0)
A.final$tx  <-  droplevels(A.final$tx)

# preprocessing
A.final <- A.final %>% mutate( race2 = ifelse (race == 'White' , 'White', 'Other'),
                              treatment.year2 = as.character(treatment.year),
                              treatment.year2 = (ifelse (treatment.year2 %in% c('2019', '2020'), '2019-2020', treatment.year2)),
                              histology2 = ifelse (grepl('Adeno', histology), 'Adenocarcinoma', 'Squamous cell'))

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
A.final <- A.final %>% mutate( time.offset = pre.tx.months)

A.final  <- A.final %>% mutate( 
     across( all_of(c(Z.count.unscaled, Ws)), function(x) (x >0), .names = "{.col}_unbinned" ),
     across( all_of(c(Z.count.unscaled)), function(x) quartile(x/time.offset), .names = "{.col}_s" )
 )
A.final  <- A.final %>% mutate( 
    across( all_of(c(X.numeric)), scale_, .names = "{.col}_z" )
)


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
 Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")
 tblcontrol <- tableby.control(numeric.stats = c('Nmiss', 'meansd'), numeric.simplify = T, cat.simplify =T, digits = 1,total = T,test = F)
 f  <-  sprintf( 'tx ~ %s', paste( c(X.numeric, X.factor,'histology', gsub('_count_s', '_count_unbinned', Zs), gsub('_count', '_count_unbinned', c(Ws))), collapse = "+"))
 labels(A.final)  <-  label_list
 tt <- tableby(as.formula(f), data=A.final, control = tblcontrol)
 summary(tt) %>% write2html(sprintf('/Users/george/Research_Local/SEER-Medicare/tbls/table1_2_%s.htm', analysis.name))

A_ = (A.final$tx == 'sbrt')*1.0
X_  <-  model.matrix(as.formula(sprintf('~ %s', paste(Xs, collapse = '+'))),  A.final)[,-1]
Z_  <- A.final[,Zs]
W1  <-  'death.other.cause.gt90day'
if (variable.selection) {
    ################################
    # Variable selection: Stage 1. For each W,  first fit a lasso to determine
    # which X and Z will be included, which will make up the Xw matrix. Then, refit
    # the model without penalization to construct W_hat.
    ################################
    nfolds  <- 5
    # Prepare the lists and matrices.
    Xw =  replicate (length(Ws)+1,   list())
    names(Xw)  <- c(W1, Ws)
    W_hat =  matrix(NA, nrow = dim(A.final)[1], ncol = length(W1) + length(Ws))
    colnames(W_hat)  <- c(W1, Ws)

    # First, for W1
    ### Select with the Cox model
    A.temp  <-  A.final %>% mutate( 
                                    W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt),
                                    W1.time  = if_else (W1.time == 0, 0.5, W1.time)/365,
                                    W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
    ) 
    fit  <-  cv.glmnet(as.matrix(cbind(A_, X_, Z_)), Surv(A.temp$W1.time, A.temp$W1.bool), family = 'cox', alpha = 1, nfolds = nfolds)
    selected.columns  <- get.selected.columns(fit, verbose=T)
    Xw[[1]]  <-  as.matrix(cbind( A_, cbind(X_, Z_ )[,selected.columns]))
    ## Refit model without penalization
    W_hat[,W1]  <-  lin_ah( time = A.temp$W1.time, 
                            event = A.temp$W1.bool, 
                            covariates = Xw[[1]]
                            )$PREDICT

    # Next, for the rest of the Ws (the counts)
    for (i in 1:length(Ws)) {
        ## Select with a poisson L1 model
        W_i  <- A.final[,Ws[i]]
        # use a glmnet to regress W_i onto A, X, Zin a lasso
        fit  <-  cv.glmnet(as.matrix(cbind(A_, X_, Z_)), as.matrix(W_i), offset= log(A.final$time.offset), family = 'poisson', alpha = 1, nfolds = nfolds)
        selected.columns  <- get.selected.columns(fit, verbose=T)
        Xw[[Ws[i]]]  <-  as.matrix(cbind(1, A_,cbind( X_, Z_ )[,selected.columns]))
        ## Refit model without penalization
        W_hat[,Ws[i]]  <- negbin_fit(y = as.matrix(W_i), 
                                     x = Xw[[Ws[i]]], 
                                     offset = log(A.final$time.offset))$PREDICT
     }

    ################################
    # Variable selection: Stage 2. For each outcome, fit a lasso to determine which
    # of the Ws and Xs to include.
    ################################
    # Next, find the variables to include in each stage 2 model. 
    # Next, for the count Ys
    selected.Ws  <- list()
    selected.Xs  <- list()
    for (outcome.i in 1:length(Ws)){ 
        outcome.name  <-  Ws[outcome.i]
        A.temp  <-  A.final %>% mutate( outcome.count  = !!rlang::sym(outcome.name))
        noc.names.temp  <- setdiff( Ws, c(outcome.name))
        #TODO: Should the time offset be in years?
        fit  <-  cv.glmnet(cbind (A_, X_, W_hat[,c('death.other.cause.gt90day', noc.names.temp)] ), as.matrix(A.temp$outcome.count), offset= log(A.temp$time.offset), family = 'poisson', alpha = 1, nfolds = nfolds)
        selected.columns  <- get.selected.columns(fit, verbose=T)
        selected.Ws[[outcome.name]] <-  Ws[check.if.present(Ws,selected.columns)]
        selected.Xs[[outcome.name]] <-  Xs[check.if.present(Xs,selected.columns)]
    }

    for (outcome.i in 1:length(outcome.names)){ 
        outcome.name  <-  outcome.names[outcome.i]
        A.temp  <- A.temp %>% mutate(
                                     Y.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                                     Y.time  = if_else (Y.time == 0, 0.5, Y.time)/365,
                                     Y.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        )
        if (outcome.name =='death.other.cause.gt90day' ) {
            W.hat.temp  <- W_hat[,c(noc.names.temp)]
        }else{
            W.hat.temp  <- W_hat[,c('death.other.cause.gt90day', noc.names.temp)]
        }
        fit  <-  cv.glmnet(as.matrix(cbind(A_, X_, W.hat.temp)), 
                           Surv(A.temp$Y.time, A.temp$Y.bool), 
                           family = 'cox', 
                           alpha = 1, 
                           nfolds = nfolds)
        coefs  <-  coef(fit, s = 'lambda.min')[,1] != 0
        selected.columns  <- names(coefs)[c(-1,-2)][coefs[c(-1,-2)]]
        selected.Ws[[outcome.name]] <-  Ws[check.if.present(Ws,selected.columns)]
        selected.Xs[[outcome.name]] <-  Xs[check.if.present(Xs,selected.columns)]
    }
    # Print out the variables selected at each stage
    sink(sprintf('data/variable.selection.seed%f.alpha%f.nfolds%f.txt', seed,alpha, nfolds))
    # First, print the Xs and Zs for each W
    for (W_i in c(W1, Ws)) {
        cat( sprintf('\n\nW: %s, with the following Xs and Zs:', W_i))
        # past the colnames
        cat(paste(colnames(Xw[[W_i]]), collapse = '\n'))
    }
    print('====================')
    for (outcome.i in 1:length(outcome.names)){ 
        outcome.name  <-  outcome.names[outcome.i]
        cat( sprintf('\n\nOutcome: %s, with the following Ws and Xs:', outcome.name))
        for (W_i in selected.Ws[[outcome.name]]) {
            print( sprintf('-->W: %s', W_i))
        }
        for (X_i in selected.Xs[[outcome.name]]) {
            print( sprintf('-->X: %s', X_i))
        }
    }
    for (outcome.i in 1:length(Ws)){ 
        outcome.name  <-  Ws[outcome.i]
        cat( sprintf('\n\nOutcome: %s, with the following Ws and Xs:', outcome.name))
        for (W_i in selected.Ws[[outcome.name]]) {
            print( sprintf('-->W: %s', W_i))
        }
        for (X_i in selected.Xs[[outcome.name]]) {
            print( sprintf('-->X: %s', X_i))
        }
    }
    sink()

}else {
    Xw = append( list(  as.matrix(cbind(  A_, X_, Z_ ))), replicate (length(Ws),   list(as.matrix(cbind(1, A_, X_, Z_ )))))
    names(Xw)  <- c(W1, Ws)
    selected.Xs  <- list()
    selected.Ws  <- list()
    for (outcome.i in 1:length(Ws)){ 
        selected.Ws[[Ws[outcome.i]]]  <-   setdiff( Ws, c(Ws[outcome.i]))
        selected.Xs[[Ws[outcome.i]]]  <-  Xs
    }
    for (outcome.i in 1:length(outcome.names)){ 
        selected.Ws[[outcome.names[outcome.i]]]  <-  Ws
        selected.Xs[[outcome.names[outcome.i]]]  <-  Xs
    }
}
################################
# Obtain estimates using the selected variables 
################################



# Time-to-event outcomes 
hazard.differences.outcomes  <-  make.odds.ratio.df ( outcome.names) 
hazard.differences.outcomes.adj  <-  make.odds.ratio.df ( outcome.names) 
hazard.differences.outcomes.proximal  <-  make.odds.ratio.df ( outcome.names) 
mouts  <-  list()
# for (outcome.i in 2:2){ 
for (outcome.i in 1:length(outcome.names)){ 
    outcome.name  <-  outcome.names[outcome.i]
    print(outcome.name)
    A.temp  <-  A.final %>% mutate( 
                                   outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                                   outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                                   outcome.time  = if_else (outcome.time == 0, 0.5, outcome.time)/ 365+ runif(dim(A.final)[1] ,0,1)*1e-8
                                   ) 
    # Raw
    f  <- as.formula('Surv(outcome.time, outcome.bool) ~ const(tx)')
    m  <-  aalen( f ,  data = A.temp, robust = 0)
    hazard.differences.outcomes[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
    # Adjust for X
    f  <-  sprintf( 'Surv(outcome.time, outcome.bool) ~ const(tx) + %s',  paste(sprintf('const(%s)', c(Xs, Zs)), collapse="+") )
    m  <-  aalen( as.formula(f) ,  data = A.temp, robust = 0)
    hazard.differences.outcomes.adj[outcome.i,1:3]  <-  c( coef(m)['const(tx)sbrt', c('Coef.', 'lower2.5%', 'upper97.5%')]) 
    # Proximal
    mout  <-  two.step.variable.selection(
                                          data.mat = A.final, 
                                          outcome.name = outcome.name,
                                          X.selected.Y  = selected.Xs[[outcome.name]],
                                          W.selected.Y  = selected.Ws[[outcome.name]],
                                          Xw = Xw,
                                          Y.count.bool =F, 
                                          verbose=T,
                                          skip.W1 = outcome.name != 'death.cause.specific'
    )
    mouts[[outcome.name]]  <-  mout
    est <- mout$ESTIMATE[1]
    se  <- mout$SE[1]
    hazard.differences.outcomes.proximal[outcome.i,1:3]  <-  c(est, est - 1.96*se, est + 1.96*se)
    print(hazard.differences.outcomes.proximal[outcome.i,1:3])
} 

# Plot cumulative hazards
mout  <- mouts[['death.cause.specific']]
largs  <-  mout$lin_ah.args
idx  <- largs$covariates[,1] ==0
cum.hazard.out.A0  <-  baseline.cum.hazard(time =largs$t1[idx], event = largs$d1[idx], covariates = largs$covariates[idx,], ESTIMATE = mout$ESTIMATE)
idx  <- largs$covariates[,1] ==1
cum.hazard.out.A1  <-  baseline.cum.hazard(time =largs$t1[idx], event = largs$d1[idx], covariates = largs$covariates[idx,], ESTIMATE = mout$ESTIMATE)
cum.hazard.toplot  <-as.data.frame(rbind(cbind(cumhaz = cum.hazard.out.A0$cumhaz, 
                                                   cum.unadj.haz = cum.hazard.out.A0$cum.unadj.haz,
                                                   t_=cum.hazard.out.A0$t_, A=0),
                               cbind(cumhaz = cum.hazard.out.A1$cumhaz, 
                                     cum.unadj.haz = cum.hazard.out.A1$cum.unadj.haz,
                                        t_=cum.hazard.out.A1$t_, A=1) ))

toplot  <-  as.data.frame(cum.hazard.toplot) %>% 
    mutate(tx = as.factor(ifelse( A == 0, 'sublobar', 'sbrt')))
g1  <-  ggplot(toplot, aes(x = t_, y = cum.unadj.haz, color = tx)) + geom_line() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))  + ggtitle('Unadjusted') + theme_bw() + labs (x = 'Time (years)', y = 'Cumulative hazard')
g2  <- ggplot(toplot, aes(x = t_, y = cumhaz, color = tx)) + geom_line() +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle('Proximal') + theme_bw()+ labs (x = 'Time (years)', y = 'Cumulative hazard')
g  <-  g1 + g2 + plot_layout(guides = "collect", axis_titles = "collect") & theme(legend.position = "bottom") 
ggsave(g, width=8, height=4, filename = sprintf('figs/%s.cumhaz.pdf', analysis.name))

g1  <-  ggplot(toplot, aes(x = t_, y = exp(-1*cum.unadj.haz), color = tx)) + geom_line() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))  + ggtitle('Unadjusted') + theme_bw() + labs (x = 'Time (years)', y = 'Cumulative hazard')
g2  <- ggplot(toplot, aes(x = t_, y = exp(-1*cumhaz), color = tx)) + geom_line() +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle('Proximal') + theme_bw()+ labs (x = 'Time (years)', y = 'Survival')
g  <-  g1 + g2 + plot_layout(guides = "collect", axis_titles = "collect") & theme(legend.position = "bottom") 
ggsave(g, width=8, height=4, filename = sprintf('figs/%s.surv.pdf', analysis.name))





# Counts
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
    mout  <-  two.step.variable.selection(
                                          data.mat = A.final, 
                                          outcome.name = outcome.name,
                                          X.selected.Y  = selected.Xs[[outcome.name]],
                                          W.selected.Y  = selected.Ws[[outcome.name]],
                                          Xw = Xw,
                                          Y.count.bool =T, 
                                          verbose=F)
    est <- mout$ESTIMATE['A']
    se  <- mout$SE['A']
    odds.ratios.nocs.proximal[outcome.i,1:3]  <-  exp(c( est, est- 1.96*se, est + 1.96*se ))
    print(odds.ratios.nocs[outcome.i,1:3])
    print(odds.ratios.nocs.adj[outcome.i,1:3])
    print(odds.ratios.nocs.proximal[outcome.i,1:3])
    cat('==========\n\n\n')
}


# old.method  <-  readRDS('data/nb2020.hazard.differences.outcomes.proximal.rds')
# hazard.differences.outcomes.proximal[,'high_ci'] - hazard.differences.outcomes.proximal[,'low_ci']
# old.method[,'high_ci'] - old.method[,'low_ci']
# old.method  <-  readRDS('data/nb2020.odds.ratios.nocs.proximal.rds')
# odds.ratios.nocs.proximal[,'high_ci'] - odds.ratios.nocs.proximal[,'low_ci']
# old.method[,'high_ci'] - old.method[,'low_ci']



## table( A.temp$arthropathy_pre_count, useNA="ifany")
A.final %>% count(tx)
#################################
## Individual plots
#################################
saveRDS(hazard.differences.outcomes, sprintf('data/%s.hazard.differences.outcomes.raw.rds', analysis.name))
saveRDS(hazard.differences.outcomes.adj, sprintf('data/%s.hazard.differences.outcomes.adj.rds', analysis.name))
fofo  <- readRDS(sprintf('data/%s.hazard.differences.outcomes.adj.rds', analysis.name))
saveRDS(hazard.differences.outcomes.proximal, sprintf('data/%s.hazard.differences.outcomes.proximal.rds', analysis.name))
saveRDS(odds.ratios.nocs, sprintf('data/%s.odds.ratios.nocs.raw.rds', analysis.name))
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
for (i in 1:nrow(toprint)) {
    cat( sprintf( '%s| %.3f (%.3f, %.3f)\n', label_list3[rownames(toprint)[i]], toprint[i,1], toprint[i,2], toprint[i,3]) ) }

toprint  <- odds.ratios.nocs.proximal
for (i in 1:nrow(toprint)) {
    cat( sprintf( '%s| %.2f (%.2f, %.2f)\n', label_list3[rownames(toprint)[i]], toprint[i,1], toprint[i,2], toprint[i,3]) ) }


