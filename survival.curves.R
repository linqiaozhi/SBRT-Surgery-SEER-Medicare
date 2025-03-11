# library(pci2s)
devtools::load_all('/Users/george/Research_Local/pci2s_gcl/pci2s')
# devtools::load_all('/Users/george/Research_Local/pci2s/')
library(ggfortify)
library(survival)
library(glmnet)
library(ahaz)
library(arsenal)
library(tidyr)
library(dplyr)
library(patchwork)
library(lubridate)
library(timereg)
library('ggtext')
library(ggplot2)
source('utilities.R')
alpha  <- 1
################################
# Preprocessing 
################################
set.seed(3)

subset.name  <- 'all.gte.65'
W1s  <-  list(
              'death' = 'death.other.cause.gt90day', 
              'death.cause.specific' = 'death.other.cause.gt90day', 
              'death.other.cause.gt90day' = 'death.other.cause.gt90day', 
              'death.copd' = 'death.noncopd', 
              'death.heart' = 'death.nonheart',
              'death.stroke' = 'death.nonstroke',
              'death.noncopd.nonheart.nonstroke' = 'death.copd.heart.stroke' 
)
outcome.names  <- names(W1s)
tt.min <- 90
filename.in  <-  sprintf('data/A.final11.%s.RDS', subset.name)
A.final  <-  readRDS(filename.in)  %>% 
    mutate(treatment.year = year(tx.date),
           death.90.day = if_else ( ninety.day.mortality, death, as.Date(NA_character_)),
           death.cause.specific = if_else ( cause.specific.mortality == 'Death', death, as.Date(NA_character_)),
           # death.other.cause = if_else ( other.cause.mortality == 'Death', death, as.Date(NA_character_)),
           death.other.cause.gt90day = if_else ( other.cause.mortality == 'Death' & tt > tt.min, death, as.Date(NA_character_)),
           death.copd = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE == '50130' & tt > tt.min ,  death, as.Date(NA_character_)),
           death.heart = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE == '50060'& tt > tt.min , death, as.Date(NA_character_)),
           death.other = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE == '50300' & tt > tt.min , death, as.Date(NA_character_)),
           death.nonother = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50300'& tt > tt.min , death, as.Date(NA_character_)),
           death.stroke = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE == '50080'& tt > tt.min , death, as.Date(NA_character_)),
           death.nonstroke = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50080'& tt > tt.min , death, as.Date(NA_character_)),
           death.noncopd = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50130' & tt > tt.min ,  death, as.Date(NA_character_)),
           death.nonheart = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50060'& tt > tt.min , death, as.Date(NA_character_)),
           death.noncopd.nonheart = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50130' & COD_TO_SITE_RECODE != '50060'& tt > tt.min  , death, as.Date(NA_character_)),
           death.noncopd.nonheart.nonstroke = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50130' & COD_TO_SITE_RECODE != '50060' & COD_TO_SITE_RECODE != '50080'& tt > tt.min  , death, as.Date(NA_character_)),
           death.copd.heart.stroke = if_else ( other.cause.mortality == 'Death' & ( COD_TO_SITE_RECODE == '50130' | COD_TO_SITE_RECODE == '50060' | COD_TO_SITE_RECODE == '50080' ) & tt > tt.min  , death, as.Date(NA_character_)),
           pre.tx.days = pre.tx.months * 30.5,
           post.tx.months = tt/30.5
    )
A.final <- A.final %>% mutate( race2 = ifelse (race == 'White' , 'White', 'Other'),
                              treatment.year2 = as.character(treatment.year),
                              treatment.year2 = (ifelse (treatment.year2 %in% c('2019', '2020'), '2019_2020', treatment.year2)),
                              treatment.year2 = (ifelse (treatment.year2 %in% c('2010', '2011'), '2010_2011', treatment.year2)),
                              histology2 = ifelse (grepl('Adeno', histology), 'Adenocarcinoma', 'Squamous cell'))

# Define X
X.factor  <-  c('sex', 'race2',  'histology2', 'treatment.year2' ) # X.factor  <-  c('sex', 'race',  'histology2') 
X.numeric  <-  c('age', 'size') 
Xs  <-  c(sprintf('%s_z', X.numeric),  X.factor)

Z.count  <- c('O2accessories', 'walking_aids' ,  'wheelchairs_accessories' , 'transportation_services', 'other_supplies',   'pressure_ulcer', 'ischemic_heart_disease', 'CHF', 'PVD', 'CVD',    'MILDLD','MSLD', 'DIAB_UC', 'DIAB_C',  'RD', 'mental_disorders', 'nervous_system',    'echo',  'Anticoags',  'smoking', 'o2',  'pneumonia_and_influenza','asthma', 'COPD','interstitial_lung')
Z.count.unscaled = sprintf( '%s_pre_12months_count', Z.count )
# Zs  <-   sprintf( '%s_pre_12months_count_bool', Z.count )
Zs  <-   sprintf( '%s_pre_12months_count', Z.count )
Q.count  <-  c( 'fall',  'other_injury', 'diverticular_disease', 'hernia',  'arthropathy','GU_sx')
Qs.unscaled = sprintf( '%s_pre_12months_count', Q.count )

A.final  <- A.final %>% mutate( 
                               time.offset = pre.tx.months, 
                               across( all_of(c(Z.count.unscaled, Qs.unscaled)), function(x) (x >0), .names = "{.col}_bool" ),
                               across( all_of(c(X.numeric)), scale_, .names = "{.col}_z" ))

#################################
# # variable selection for W1
################################
exclude.from.stage1  <-  c('size_z', 'histology2')
X_  <-  model.matrix(as.formula(sprintf('~ %s', paste(Xs[!(Xs) %in% exclude.from.stage1], collapse = '+'))),  A.final)[,-1]
Z_  <- A.final[,Zs]
A_ = (A.final$tx == 'sbrt')*1.0
design.mat  <-  as.matrix(cbind(A_, X_, Z_))
selected.columns <- list()
for (W1i in 1:length(W1s)) {
    outcome.name  <- W1s[[W1i]]
    A.temp  <-  A.final %>% mutate( 
                                   outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                                   outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                                   outcome.time  = if_else (outcome.time == 0, 0.5, outcome.time)/ 365+ runif(dim(A.final)[1] ,0,1)*1e-8
    ) 
    mm  <-   tune.ahazpen( Surv(A.temp$outcome.time, A.temp$outcome.bool), design.mat )
    selected.columns[[outcome.name]]  <- get.selected.columns.ahaz(mm, s = 'lambda.min', colnames(design.mat), verbose=T, min.vars = 1)
    print(sprintf( 'Included variables in stage 1 for %s: %s', outcome.name, paste( selected.columns[[outcome.name]], collapse = ', ')))
}

################################
# Calculate survival curves
################################
# Time-to-event outcomes 
outcome.i  <-  2
outcome.name  <-  outcome.names[outcome.i]
W1  <- W1s[[outcome.i]]
W1.selected.variables = selected.columns[[W1]]
print(outcome.name)
A.temp  <-  A.final %>% mutate( 
                               outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                               outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                               outcome.time  = if_else (outcome.time == 0, 0.5, outcome.time)/ 365+ runif(dim(A.final)[1] ,0,1)*1e-8
) 
A.temp  <-  A.temp %>% mutate( 
                              W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt),
                              W1.time  = if_else (W1.time == 0, 0.5, W1.time)/365,
                              W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
                              Y.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                              Y.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
)
A_ = (A.temp$tx == 'sbrt')*1.0
X_  <-  model.matrix(as.formula(sprintf('~ %s', paste(Xs, collapse = '+'))),  A.temp)[,-1]
W_ <- as.matrix( A.temp$W1.time)
D2_  <- A.temp$W1.bool*1.0
Z_  <- A.temp[,Zs]
N_  <- dim(A.temp)[1]
Xw =  list(  as.matrix(cbind(  A_, X_[,colnames(X_) %in% W1.selected.variables], Z_[,colnames(Z_) %in% W1.selected.variables] )))
Y_ = A.temp$Y.time
D_ = A.temp$Y.bool*1.0


################################
# Point estimate 
################################

Y_ = A.temp$Y.time
D_ = A.temp$Y.bool*1.0
mout  <- p2sls.ah(Y = Y_, D = D_,  A = A_, X = X_,
                      W = W_, 
                      Xw = Xw,        
                      nco_type = c("ah"),
                      nco_args = list(list(offset = rep(0, N_), event = D2_)) )
est <- mout$ESTIMATE[1]
se  <- mout$SE[1]
c(est, est - 1.96*se, est + 1.96*se)

################################
# Survival function 
################################
A.temp %>% print(width=Inf)

A.temp <- A.temp %>% mutate(
                            cause = case_when (
                                               cause.specific.mortality == 'Death' ~ 0,
                                               other.cause.mortality == 'Death' ~ 1,
                                               T ~ -1
                                               )
                            )


A.temp %>% count(cause, cause.specific.mortality, other.cause.mortality)

# times = Y_
# cause = A.temp$cause
# A = A_
#  a = 1
#  X = X_
#  Z = Z_
#  nc_time = 90/365
#  nt = 1000


survfunc_a0 <- p2sls.cprisk.nc(times = Y_, cause = A.temp$cause, A = A_, a = 0, X = X_, Z = Z_, nc_time = 90/365)
survfunc_a1 <- p2sls.cprisk.nc(times = Y_, cause = A.temp$cause, A = A_, a = 1, X = X_, Z = Z_, nc_time = 90/365)
# p2sls_rslt_a1 <- p2sls.cprisk(times = Y_, cause = A.temp$cause, A = A_, a = 1, X = X_, Z = Z_)
# survfunc_a1  <- p2sls.ah.survfunc( Y = Y_, D = D_,  A = A_, a = 1, X = X_,
#                                   W = W_, Z = Z_,
#                                   Xw = Xw,        
#                                   nco_type = c("ah"),
#                                   nco_args = list(list(offset = rep(0, N_), event = D2_)) )
# survfunc_a0  <- p2sls.ah.survfunc( Y = Y_, D = D_,  A = A_, a = 0, X = X_,
#                                   W = W_, Z = Z_,
#                                   Xw = Xw,        
#                                   nco_type = c("ah"),
#                                   nco_args = list(list(offset = rep(0, N_), event = D2_)) )
survfunc_a0$strata <- 'Surgery'
survfunc_a1$strata <- 'SBRT'
surv.curves.out  <-  rbind (a1 = survfunc_a1, a0 = survfunc_a0)
surv.curves.out <- surv.curves.out %>% mutate(time = t, surv=survfunc)
surv.curves.out$strata  <- factor(surv.curves.out$strata, levels = c('Surgery', 'SBRT'))
surv.curv.2  <-  function( surv.curves.out ) {
    ggplot (surv.curves.out, aes(x = t, y = cif0, color = strata)) + geom_line() + xlab("Time") + ylab("Cumulative incidence") + scale_color_manual(values = c("#E7B800", "#2E9FDF")) + theme_minimal() + ggtitle('Cumulative incidence function (Proximally-adjusted)') + theme(legend.title=element_blank()) + ylim(0,1)
}
gp  <- surv.curv.2(surv.curves.out)
ggsave(plot=gp, filename='figs/cifnctime90.pdf', width=7, height=5)

survfunc_a0 <- p2sls.cprisk.nc(times = Y_, cause = A.temp$cause, A = A_, a = 0, X = X_, Z = Z_, nc_time = 1/365)
survfunc_a1 <- p2sls.cprisk.nc(times = Y_, cause = A.temp$cause, A = A_, a = 1, X = X_, Z = Z_, nc_time = 1/365)
# p2sls_rslt_a1 <- p2sls.cprisk(times = Y_, cause = A.temp$cause, A = A_, a = 1, X = X_, Z = Z_)
# survfunc_a1  <- p2sls.ah.survfunc( Y = Y_, D = D_,  A = A_, a = 1, X = X_,
#                                   W = W_, Z = Z_,
#                                   Xw = Xw,        
#                                   nco_type = c("ah"),
#                                   nco_args = list(list(offset = rep(0, N_), event = D2_)) )
# survfunc_a0  <- p2sls.ah.survfunc( Y = Y_, D = D_,  A = A_, a = 0, X = X_,
#                                   W = W_, Z = Z_,
#                                   Xw = Xw,        
#                                   nco_type = c("ah"),
#                                   nco_args = list(list(offset = rep(0, N_), event = D2_)) )
survfunc_a0$strata <- 'Surgery'
survfunc_a1$strata <- 'SBRT'
surv.curves.out  <-  rbind (a1 = survfunc_a1, a0 = survfunc_a0)
surv.curves.out <- surv.curves.out %>% mutate(time = t, surv=survfunc)
surv.curves.out$strata  <- factor(surv.curves.out$strata, levels = c('Surgery', 'SBRT'))
surv.curv.2  <-  function( surv.curves.out ) {
    ggplot (surv.curves.out, aes(x = t, y = cif0, color = strata)) + geom_line() + xlab("Time") + ylab("Survival probability") + scale_color_manual(values = c("#E7B800", "#2E9FDF")) + theme_minimal() + ggtitle('Cumulative incidence function (Proximally-adjusted)') + theme(legend.title=element_blank()) + ylim(0,1)
}
gp  <- surv.curv.2(surv.curves.out)
ggsave(plot=gp, filename='figs/cifnctime1.pdf', width=7, height=5)








library(survminer)
surv.curv  <-  function(fit, dframe , title, risk.table=T){
    custom.risk.table.theme <- theme_survminer() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),  axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), legend.position = "none", legend.title=element_blank()) 
    ggsurvplot(fit, data = dframe,
     # surv.median.line = "hv", # Add medians survival
    # Change legends: title & labels
     title =title,
     legend.labs = c("Surgery", "SBRT"),
     # legend = ifelse (risk.table, "top", "none"),
      legend = "none",
     censor = F,
     # Add p-value and tervals
     pval = F,
    conf.int = TRUE,
     # Add risk table
     risk.table = risk.table,
     tables.height = 0.25,
      tables.theme = custom.risk.table.theme,
     fontsize = 3,
     risk.table.title = '',
    # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
     # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
     palette = c("#E7B800", "#2E9FDF"),
     ggtheme = theme_bw(), # Change ggplot2 theme
     break.time.by = 2.5
    ) 
}
dframe  <-  data.frame(time = A.final$tt/365, status = nna(A.final$death), A=A.temp$tx)
fit  <- survfit(Surv(time, status) ~ A, data = dframe)
# Customized survival curves
g1  <- surv.curv(fit, dframe, 'Overall mortality', risk.table =F)
dframe  <-  data.frame(time = A.temp$tt/365, status = nna(A.temp$death.cause.specific), A=A.temp$tx)
fit  <- survfit(Surv(time, status) ~ A, data = dframe)
g2  <- surv.curv(fit, dframe, 'Cause-specific mortality', risk.table =F)
# Customized survival curves
dframe  <-  data.frame(time = A.temp$tt/365, status = nna(A.temp$other.cause.mortality), A=A.temp$tx)
fit  <- survfit(Surv(time, status) ~ A, data = dframe)
g3  <- surv.curv(fit, dframe, 'Other cause mortality')
# Customized survival curves
# g <- arrange_ggsurvplots(list(g1,g2,g3), nrow=1, ncol=3, print=F)
# ggsave(plot=g, filename='figs/raw.curves.pdf', width=12, height=5)
g <- arrange_ggsurvplots(list(g1,g2,g3), nrow=3, ncol=1, print=F)
ggsave(plot=g, filename='figs/raw.curves.pdf', width=5, height=12)
