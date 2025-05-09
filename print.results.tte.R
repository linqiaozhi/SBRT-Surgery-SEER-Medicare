
library(tidyverse)
source('utilities.R')
library(patchwork)
library('ggplot2')
variable.selection  <- 'automatic'
 subset.names  <- list('all.gte.65', 'sens1')
display.outcomes  <- c(  'death.copd','death.heart',    'death.stroke', 'death.noncopd.nonheart.nonstroke', sprintf('death.other.%d', 1:4))
 lambda.ss  <- list('lambda.min')
nc_time_days = 90
for (subset.name in subset.names ){
    # analysis.name  <- sprintf('v10.tte.ttmin0.%s', subset.name)
    analysis.name  <- sprintf('data35.mc.v1.vsF.nctime%d_%s', nc_time_days, subset.name)
    # analysis.name  <- sprintf('v42.tte.nc_time_ss_no_vs_%d_%s', nc_time_days, subset.name)
    plot.list  <- list()
    for (lambda.s in lambda.ss) {
        # analysis.name  <- sprintf('v3.tte.nocopdinW1.%s', subset.name)
        print(analysis.name)
        for (type in c('raw','adj', 'proximal') ) {
            print(sprintf(">%s - %s", analysis.name, type))
            HDs  <- readRDS(sprintf('data/%s.hazard.differences.outcomes.%s.rds', analysis.name, type))
            toprint  <- rbind(HDs)
            for (i in 1:nrow(toprint)) {
                cat( sprintf( '%s| %.4f (%.4f, %.4f)\n', rownames(toprint)[i], toprint[i,1], toprint[i,2], toprint[i,3]) ) 
            }
            print('=====================')
            type.titles  <-  list('raw'='Unadjusted', 'adj'='Adjusted', 'proximal'='Proximal')
            height  <-  2.0
            if (type == 'proximal') {
                Y.toplot  <-  c( display.outcomes) 
            }else {
                 Y.toplot  <-  c('death', 'death.lc.specific', 'death.other.cause.gt90day',display.outcomes)
                 Y.toplot  <-  c('death.other.cause',display.outcomes)
            }
            top.plot.height  <- 0.33*length(Y.toplot)
            HDs  <-  HDs[Y.toplot,]
            HDs$y_axis  <-  1:nrow(HDs)
            label_list3  <-  label_list2
            # names(label_list3)  <- gsub( '_unbinned', '', names(label_list3))
            g1.a  <-  make.HD.plot(HDs, label_list3)  + ggtitle(type.titles[[type]])
            plot.list[[type]]  <- g1.a
            ggsave(g1.a, width=6, height=top.plot.height, filename = sprintf('figs/%s.%s.pdf', analysis.name, type))
            # Print the figure
        }
        print('=====================')
        print('=====================')
    }
    g1  <-   wrap_plots(plot.list, nrow =3, ncol = 1, heights = c(1,1,0.9)) + plot_layout(guides = "collect", axes = "collect")
    ggsave(g1, width=6, height=6, filename = sprintf('figs/%s.pdf', analysis.name))

    HDs  <- list()
    for (type in c('raw','adj', 'proximal') ) {
        HDs[[type]]  <- readRDS(sprintf('data/%s.hazard.differences.outcomes.%s.rds', analysis.name, type))
    }
     HD  <- rbind('raw' = HDs[['raw']]['death.lc.specific',],
                  'adj' = HDs[['adj']]['death.lc.specific',] ,
                  'proximal' = HDs[['proximal']]['death.lc.specific',] 
     )
     HD$y_axis   <- 1:nrow(HD)
     HD$'death.lc.specific'  <- rownames(HD)
     label_list_types  <- list('raw'='Unadjusted', 'adj'='Adjusted', 'proximal'='**Proximal**')
     g <-  make.HD.plot(HD, label_list_types)  + ggtitle('Cause-specific mortality')
     ggsave(g, width=6, height=1.2, filename = sprintf('figs/%s.lc.specific.pdf', analysis.name))
}





