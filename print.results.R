source('utilities.R')
library(patchwork)
library('ggplot2')
variable.selection  <- 'automatic'
subset.names  <- list('all.gte.65', 'sens1')
lambda.ss  <- list('lambda.min', 'lambda.min.lowest', 'lambda.min.lower')
for (subset.name in subset.names ){
    for (lambda.s in lambda.ss) {
        analysis.name  <- sprintf('%s.%s', subset.name, lambda.s)
        print(analysis.name)
        for (type in c('raw','adj', 'proximal') ) {
            print(sprintf(">%s - %s", analysis.name, type))
            HDs  <- readRDS(sprintf('data/%s.hazard.differences.outcomes.%s.rds', analysis.name, type))
            ORs  <- readRDS(sprintf('data/%s.odds.ratios.nocs.%s.rds', analysis.name, type))
            toprint  <- rbind(HDs, ORs)
            for (i in 1:nrow(toprint)) {
                cat( sprintf( '%s| %.3f (%.3f, %.3f)\n', rownames(toprint)[i], toprint[i,1], toprint[i,2], toprint[i,3]) ) 
            }
            print('=====================')
            type.titles  <-  list('raw'='Unadjusted', 'adj'='Adjusted', 'proximal'='Proximal')
            height  <-  1.8
            if (type == 'proximal') {
                Y.toplot  <-  c( 'death.cause.specific')
                top.plot.height  <- 0.25
            }else {
                Y.toplot  <-  c('death', 'death.90.day', 'death.other.cause.gt90day', 'death.cause.specific')
                top.plot.height  <- 1
            }
            HDs  <-  HDs[Y.toplot,]
            HDs$y_axis  <-  1:nrow(HDs)
            label_list3  <-  label_list2
            names(label_list3)  <- gsub( '_unbinned', '', names(label_list3))
            g1.a  <-  make.HD.plot(HDs, label_list3) 
            g1.b  <-  make.OR.plot(ORs, label_list3, number.of.spaces = 30) 
            x_start <- 0.9
            x_end <- 1
            y_position <- 10
            g1  <-  g1.a / g1.b+ plot_layout(heights = (c(top.plot.height,height))) + plot_annotation(title=type.titles[[type]])
            g1
            ggsave(g1, width=6, height=2.8, filename = sprintf('figs/%s.%s.pdf', analysis.name, type))
            # Print the figure
        }
        print('=====================')
        print('=====================')
    }
}

