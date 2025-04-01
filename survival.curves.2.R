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
source('load.data.R')
# library('ggtext')
library(ggplot2)
source('utilities.R')
################################
# Preprocessing 
################################
set.seed(3)


subset.name  <- 'all.gte.65'
nc_time_days = 90
sm  <- load.data(subset.name, nc_time_days = nc_time_days)

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
dframe  <-  data.frame(time = sm$A.final$tt/365, status = nna(sm$A.final$death), A=sm$A.final$tx)
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
A.temp  <-  sm$A.final %>% mutate( 
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
X_  <-  model.matrix(as.formula(sprintf('~ %s', paste(sm$Xs, collapse = '+'))),  A.temp)[,-1]
Z_  <- A.temp[,sm$Zs]
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

survfunc_a0 <- p2sls.cprisk.nc.cif(times = Y_, cause = A.temp$cause, A = A_, a = 0, X = X_, Z = Z_, nc_time = nc_time_days/365)
survfunc_a1 <- p2sls.cprisk.nc.cif(times = Y_, cause = A.temp$cause, A = A_, a = 1, X = X_, Z = Z_, nc_time = nc_time_days/365)

survfunc_a0$strata <- 'Surgery'
survfunc_a1$strata <- 'SBRT'
surv.curves.out  <-  rbind (a1 = survfunc_a1, a0 = survfunc_a0)
surv.curves.out <- surv.curves.out %>% mutate(time = t, surv=survfunc)

surv.curves.out$strata  <- factor(surv.curves.out$strata, levels = c('SBRT','Surgery' ))
surv.curv.2  <-  function( surv.curves.out ) {
    ggplot (surv.curves.out, aes(x = t, y = cif0, color = strata)) + geom_line(linewidth=1) + xlab("Time") + ylab("Cumulative incidence") + scale_color_manual(values = c("#E7B800", "#2E9FDF")) + theme_minimal() + ggtitle('Cumulative incidence function (Proximally-adjusted)') + theme(legend.position = "bottom", legend.title=element_blank()) + ylim(0,0.6)
}
gp  <- surv.curv.2(surv.curves.out) + g1$plot$theme
ggsave(plot=gp, filename='figs/cifnctime90.pdf', width=5, height=4)

