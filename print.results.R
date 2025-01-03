
source('utilities.R')
library(patchwork)
library('ggplot2')
variable.selection  <- 'automatic'
# subset.names  <- list('all.gte.65', 'sens1', 'sens.pet', 'sens.ebus')
 subset.names  <- list('all.gte.65', 'sens1')
   # subset.names  <- list('all.gte.65')
# lambda.ss  <- list('lambda.min', 'lambda.min.lowest', 'lambda.min.lower')
display.outcomes  <- list(
                          'death.noncopd' = c('death.copd'),
                          'death.nonheart' = c('death.heart'),
                          'death.noncopd.nonheart' = c('death.heart', 'death.copd'),
                          'death.stroke' = c('death.heart', 'death.copd'),
                          'death.copd' = c('death.heart', 'death.noncopd', 'death.stroke'),
                          'death.nonstroke' = c('death.heart',  'death.stroke'),
                          'death.other.cause.gt90day' = c()
                          )
display.outcomes  <- c('death.other.cause.gt90day',  'death.copd','death.heart',   'death.other', 'death.stroke')
 lambda.ss  <- list('lambda.min')
for (subset.name in subset.names ){
    for (lambda.s in lambda.ss) {
        analysis.name  <- sprintf('v2.tte.%s', subset.name)
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
            height  <-  2.0
            if (type == 'proximal') {
                Y.toplot  <-  c( 'death.cause.specific', display.outcomes) 
            }else {
                # Y.toplot  <-  c('death', 'death.90.day',  'death.cause.specific',  'death.heart', 'death.copd', 'death.noncopd.nonheart')
                 Y.toplot  <-  c('death', 'death.cause.specific', display.outcomes)
            }
            top.plot.height  <- 0.25*length(Y.toplot)
            HDs  <-  HDs[Y.toplot,]
            HDs$y_axis  <-  1:nrow(HDs)
            label_list3  <-  label_list2
            # names(label_list3)  <- gsub( '_unbinned', '', names(label_list3))
            g1.a  <-  make.HD.plot(HDs, label_list3) 
            g1.b  <-  make.OR.plot(ORs, label_list3, number.of.spaces = 30) 
            x_start <- 0.9
            x_end <- 1
            y_position <- 10
            g1  <-  g1.a / g1.b+ plot_layout(heights = (c(top.plot.height,height))) + plot_annotation(title=type.titles[[type]])
            ggsave(g1, width=6, height=top.plot.height + height, filename = sprintf('figs/%s.%s.pdf', analysis.name, type))
            # Print the figure
        }
        print('=====================')
        print('=====================')
    }
}

# this makes sense
# death.cause.specific| 0.052 (0.043, 0.061)

# sens1 = raw
# death.cause.specific| 0.050 (0.041, 0.059)


# subset.name  <- 'all.gte.65'
# # lambda.ss  <- 0:5
# type  <- 'proximal'
# toplot  <- data.frame()
# for (lambda.s in lambda.ss) {
#     analysis.name  <- sprintf('%s.%s', subset.name, lambda.s)
#     print(sprintf(">%s - %s", analysis.name, type))
#     HDs  <- readRDS(sprintf('data/%s.hazard.differences.outcomes.%s.rds', analysis.name, 'proximal'))
#     toplot  <- rbind(toplot, HDs['death.cause.specific',])
# }
# toplot$lambdas  <- lambda.ss
# library(ggplot2)
# lambda.ss.labels  <-  sprintf('λ + %dδ', lambda.ss)
# lambda.ss.labels[1]  <-  'λ'
# lambda.ss.labels[length(lambda.ss.labels)]  <-  'lambda.min.mse'
# g  <-  ggplot(toplot, aes(x = lambdas, y = estimate, ymin = low_ci, ymax = high_ci)) +  geom_pointrange() + scale_x_continuous (labels = lambda.ss.labels)  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab('Hazard Difference') + xlab('Regularization strength') + ylim(-0.01, 0.05)

