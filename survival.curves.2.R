# library(pci2s)
# devtools::load_all('/Users/george/Research_Local/pci2s_gcl/pci2s')
devtools::load_all('/Users/george/Research_Local/pci2s/')
library(ggfortify)
library(survival)
library(tidycmprsk)
library(ggsurvfit)
library(glmnet)
# library(tidyr)
library(dplyr)
library(patchwork)
library(lubridate)
library(timereg)
# library('ggtext')
library(ggplot2)
source('utilities.R')
################################
# Preprocessing 
################################
set.seed(3)

subset.name  <- 'all.gte.65'
nc_time_days = 90
tt.min <- nc_time_days
W1s  <-  list(
            # 'death' = 'death.other.cause.gt90day', 
            'death.cause.specific' = 'death.other.cause.gt90day', 
             'death.other.cause.gt90day' = 'death.other.cause.gt90day', 
            'death.copd' = 'death.noncopd', 
            'death.heart' = 'death.nonheart',
            'death.stroke' = 'death.nonstroke',
             'death.noncopd.nonheart.nonstroke' = 'death.copd.heart.stroke' 
            )
outcome.names  <- names(W1s)
# After tt.min (days), W is a valid negative outcome
analysis.name  <- sprintf('v12.tte.nc_time_ss_%d_%s', nc_time_days, subset.name)
print('====================')
print(analysis.name)
print('====================')

################################
# Load data 
################################
filename.in  <-  sprintf('data/A.final12.%s.RDS', subset.name)
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
A.final %>%   
    mutate(COD = ifelse( COD_TO_SITE_RECODE %in% cod.df$COD_TO_SITE_RECODE, COD_TO_SITE_RECODE, '50300')) %>%
    group_by(COD) %>% summarise(  n = n(), p = n()/nrow(.) ) %>% arrange(desc(n)) %>%   left_join(cod.df, by = c('COD'='COD_TO_SITE_RECODE'))  %>% mutate (np = sprintf('%d (%.1f%%)', n ,p*100 ))%>% select( COD, Name,  np) %>% print(n =10)
A.final %>% summarise(sum(other.cause.mortality == 'Death'),  mean(other.cause.mortality == 'Death') *100) 

table( A.final$tx, useNA="ifany")
# For the sensitivity analysis, some node positive patients are included
A.final$tnm.n[is.na(A.final$tnm.n)]  <- 'X'
A.final  <- A.final  %>% filter (tnm.n %in% c('0', '1', '2')) 
A.final  <- A.final %>% filter (tx == 'sbrt' |
                                ( tx == 'sublobar' & tnm.n == '0' ) | 
                                (tx == 'sublobar' & tnm.n != '0' & REGIONAL_NODES_EXAMINED_1988 != '00' ) )
print(sprintf('%.3f%% of the sublobar patients are N+', 100*sum(A.final$tnm.n != 0 & A.final$tx == 'sublobar')/ sum(A.final$tx == 'sublobar')))
A.final$tx  <-  droplevels(A.final$tx)

# preprocessing
A.final <- A.final %>% mutate( race2 = ifelse (race == 'White' , 'White', 'Other'),
                              treatment.year2 = as.character(treatment.year),
                              treatment.year2 = (ifelse (treatment.year2 %in% c('2019', '2020'), '2019_2020', treatment.year2)),
                              treatment.year2 = (ifelse (treatment.year2 %in% c('2010', '2011'), '2010_2011', treatment.year2)),
                              histology2 = ifelse (grepl('Adeno', histology), 'Adenocarcinoma', 'Squamous cell'))

# Define X
X.factor  <-  c('sex', 'race2',  'histology2', 'treatment.year2' ) 
X.numeric  <-  c('age', 'size') 
Xs  <-  c(sprintf('%s_z', X.numeric),  X.factor)

Z.count  <- c('O2accessories', 'walking_aids' ,  'wheelchairs_accessories' , 'transportation_services', 'other_supplies',   'pressure_ulcer', 'ischemic_heart_disease', 'CHF', 'PVD', 'CVD',    'MILDLD','MSLD', 'DIAB_UC', 'DIAB_C',  'RD', 'mental_disorders', 'nervous_system',    'echo',  'Anticoags',  'smoking', 'o2',  'pneumonia_and_influenza','asthma', 'COPD','interstitial_lung')
Z.count.unscaled = sprintf( '%s_pre_12months_count', Z.count )
 Zs  <-   sprintf( '%s_pre_12months_count', Z.count )

A.final  <- A.final %>% mutate( 
                     time.offset = pre.tx.months, 
                               across( all_of(c(Z.count.unscaled)), function(x) (x >0), .names = "{.col}_bool" ),
                               # across( all_of(c(Z.count.unscaled)), function(x) quartile(x), .names = "{.col}_s" ),
                               across( all_of(c(X.numeric)), scale_, .names = "{.col}_z" ))

################################
# Overall survival curve
################################
# Without proximal adjustment
library(survminer)
surv.curv  <-  function(fit, dframe , title, risk.table=T){
    custom.risk.table.theme <- theme_survminer() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),  axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), legend.position = "none", legend.title=element_blank()) 
    ggsurvplot(fit, data = dframe,
               title =title,
               legend.labs = c("Surgery", "SBRT"),
               legend = "none",
               censor = F,
               pval = F,
               conf.int = F,
               # Add risk table
               risk.table = risk.table,
               tables.height = 0.25,
               tables.theme = custom.risk.table.theme,
               fontsize = 3,
               risk.table.title = '',
               palette = c("#E7B800", "#2E9FDF"),
               ggtheme = theme_bw(), # Change ggplot2 theme
               break.time.by = 2.5
    ) 
}
dframe  <-  data.frame(time = A.final$tt/365, status = nna(A.final$death), A=A.final$tx)
fit  <- survfit(Surv(time, status) ~ A, data = dframe)
# Customized survival curves
g1  <- surv.curv(fit, dframe, 'Overall mortality', risk.table =T)
risk.table  <- g1$table

################################
# Cause-specific and other-cause CIFs
################################
# Time-to-event outcomes 
outcome.name  <-  'death.cause.specific'
W1  <- 'death.other.cause.gt90day'
A.temp  <-  A.final %>% mutate( 
                               W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt),
                               W1.time  = if_else (W1.time == 0, 0.5, W1.time)/365,
                               W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
                               Y.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                               Y.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
                               cause = case_when (
                                                  nna(!!rlang::sym(outcome.name))~ 0, # primary event
                                                  nna(!!rlang::sym(W1))~ 1,
                                                  T ~ -1
                               )
)

A_ = (A.temp$tx == 'sbrt')*1.0
X_  <-  model.matrix(as.formula(sprintf('~ %s', paste(Xs, collapse = '+'))),  A.temp)[,-1]
Z_  <- A.temp[,Zs]
Y_ = A.temp$Y.time

################################
# Cumulative incidence functions 
################################
dframe  <-  data.frame( time = A.temp$tt/365, status = A.temp$cause, A=A.temp$tx )  
dframe$A  <-  case_match(dframe$A, 
                         'sbrt' ~ 'SBRT',
                         'sublobar' ~ 'Surgery')
dframe$status  <-  case_match(
                              dframe$status, 
                              0 ~ 'Cause-specific death', 
                              1 ~ 'Other-cause death',
                              -1 ~ 'Censored'
                              ) %>% as.factor
dframe$status  <- factor(dframe$status, levels=c('Censored', 'Cause-specific death', 'Other-cause death'))
cuminc.obj <- cuminc(Surv(time, status) ~ A, data = dframe) 
g2   <- ggcuminc(cuminc.obj, outcome = "Cause-specific death", linewidth = 1)  + ggtitle ("Cause-specific death")+ g1$plot$theme + theme(legend.position = 'none') + ylim(0, 0.6)  + scale_colour_manual( values=c("#E7B800", "#2E9FDF"))
g3  <-  ggcuminc(cuminc.obj, outcome = "Other-cause death", linewidth = 1)+ ggtitle ("Other-cause death") + ylim(0, 0.6) + g1$plot$theme + theme(legend.position = 'bottom') + scale_colour_manual( values=c("#E7B800", "#2E9FDF")) 
g  <-  g1$plot / g2 /g3 + plot_layout(heights = c(1, 1, 1), axes='collect')
ggsave(plot=g, filename='figs/raw.curves.pdf', width=5, height=9)
ggsave(risk.table, filename='figs/risk.table.pdf', width=5, height=1)

g  <-  g1$plot / g2 /g3 / risk.table + plot_layout(heights = c(1, 1, 1, 0.4), axes='collect')
g

################################
# Proximally adjusted CIFs
################################

survfunc_a0 <- p2sls.cprisk.nc.cif(times = Y_, cause = A.temp$cause, A = A_, a = 0, X = X_, Z = Z_, nc_time = 90/365)
survfunc_a1 <- p2sls.cprisk.nc.cif(times = Y_, cause = A.temp$cause, A = A_, a = 1, X = X_, Z = Z_, nc_time = 90/365)

survfunc_a0$strata <- 'Surgery'
survfunc_a1$strata <- 'SBRT'
surv.curves.out  <-  rbind (a1 = survfunc_a1, a0 = survfunc_a0)
surv.curves.out <- surv.curves.out %>% mutate(time = t, surv=survfunc)

surv.curves.out$strata  <- factor(surv.curves.out$strata, levels = c('SBRT','Surgery' ))
surv.curv.2  <-  function( surv.curves.out ) {
    ggplot (surv.curves.out, aes(x = t, y = cif0, color = strata)) + geom_line(linewidth=1) + xlab("Time") + ylab("Cumulative incidence") + scale_color_manual(values = c("#E7B800", "#2E9FDF")) + theme_minimal() + ggtitle('Cumulative incidence function (Proximally-adjusted)') + theme(legend.position = "bottom", legend.title=element_blank()) + ylim(0,0.6)
}
gp  <- surv.curv.2(surv.curves.out)
ggsave(plot=gp, filename='figs/cifnctime90.pdf', width=5, height=4)

