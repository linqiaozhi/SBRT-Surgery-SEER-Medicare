source('utilities.R')
library(patchwork)
library(ggplot2)





label_list3  <-  label_list2
names(label_list3)  <- gsub( '_unbinned', '', names(label_list3))
ll  <- label_list3



for ( analysis.name in c('raw', 'adj', 'proximal' ) ) {


bump  <- -0.1
bump.2  <- -0.2
error.height.1  <- 0.1
error.height.2  <- 0.3
hazard.differences.main  <-  readRDS(sprintf('data/%s.hazard.differences.outcomes.%s.rds', 'lambda.min.lower', analysis.name))[c('death.cause.specific'),]
hazard.differences.above  <- readRDS(sprintf('data/%s.hazard.differences.outcomes.%s.rds', 'lambda.min.lowest',analysis.name))[c('death.cause.specific'),]
hazard.differences.main$color   <-  'N0, less regularization'
hazard.differences.above$color  <-  'N0, least regularization'
hazard.differences.main$y_axis  <-  1:nrow(hazard.differences.main)
hazard.differences.above$y_axis  <- 1:nrow(hazard.differences.main) - bump
hazard.differences  <- rbind(hazard.differences.main, hazard.differences.above)
hazard.differences$color  <- factor(hazard.differences$color, levels = c(c('N0, less regularization', 'N0, least regularization')))
xlims = c(-0.05,0.15)
breaks=c(-0.05,0.025,0, 0.05,0.1, 0.15)
tt  <-  'Hazard Difference'
g.a <- ggplot(hazard.differences, aes(x = estimate, y=y_axis, linetype = color)) + 
geom_vline(aes(xintercept = 0), size = 0.25, linetype = "dashed") +
geom_errorbarh(aes( xmax = high_ci, xmin = low_ci), size = 0.50, height = error.height)+
# geom_point(size=1.5) +
theme_bw() +
theme(panel.grid.minor = element_blank()) +
scale_y_continuous(breaks = 1:max(hazard.differences$y_axis), labels = ll[row.names(hazard.differences.main)], trans='reverse') +
scale_x_continuous(limits = xlims, breaks = breaks ) +
#coord_trans(x = "log10") +
xlab(tt) +
ylab("") +
scale_linetype_manual(values=c("solid","dashed"))+
scale_shape_manual(values=c(15,17))+
theme( panel.grid.major.x = element_blank() ,
      panel.grid.major.y = element_line( size=.05, color="grey", linetype = 'dashed' ),
legend.position = 'bottom',
#  legend.title = element_blank(),
legend.key.size=grid::unit(2,"lines"),
 axis.text.y = ggtext::element_markdown())+
guides(linetype="none")

odds.ratios.main  <-  readRDS(sprintf('data/%s.odds.ratios.nocs.%s.rds', 'lambda.min.lower',analysis.name))
odds.ratios.above  <- readRDS(sprintf('data/%s.odds.ratios.nocs.%s.rds', 'lambda.min.lowest',analysis.name))
odds.ratios.main$color   <-  'N0, less regularization'
odds.ratios.above$color  <-  'N0, least regularization'
#reorder the color
odds.ratios.main$y_axis  <- 1:nrow(odds.ratios.main) 
odds.ratios.above$y_axis  <- 1:nrow(odds.ratios.main) - bump.2
odds.ratios  <- rbind(odds.ratios.main, odds.ratios.above)
odds.ratios$color  <- factor(odds.ratios$color, levels = c('N0, less regularization', 'N0, least regularization'))

xlims <- c(0.3, 3)
tt  <-  'Risk ratio (log scale)'
row.names(odds.ratios) <- gsub("_..._count", "",  row.names(odds.ratios))

g.b <- ggplot(odds.ratios, aes(x = estimate, y=y_axis, linetype=color)) + 
geom_vline(aes(xintercept = 1), size = 0.25, linetype = "dashed") +
geom_errorbarh(aes( xmax = high_ci, xmin = low_ci), size = 0.50, height = error.height.2)+
theme_bw() +
theme(panel.grid.minor = element_blank()) +
scale_y_continuous(breaks = 1:max(odds.ratios.main$y_axis), labels = ll[row.names(odds.ratios.main)], trans='reverse') +
scale_x_continuous(limits = xlims ) +
coord_trans(x = "log10") +
xlab(tt) +
ylab("") +
scale_linetype_manual(values=c("solid","dashed"))+
scale_shape_manual(values=c(15,17))+
theme( panel.grid.major.x = element_blank() ,
      panel.grid.major.y = element_line( size=.05, color="grey", linetype = 'dashed' ),
legend.position = 'bottom',
  legend.title = element_blank(),
legend.key.size=grid::unit(2,"lines"),
 axis.text.y = ggtext::element_markdown()) 






g  <-  g.a / g.b + plot_layout(heights = (c(0.25,1.5))) + plot_annotation(title="Proximal")
g
ggsave(g, width=6.5, height=3.5, filename = sprintf('figs/supplemental.%s.pdf', analysis.name))

}






for ( analysis.name in c('raw', 'adj', 'proximal' ) ) {
    print(analysis.name)
    print(readRDS(sprintf('data/%s.hazard.differences.outcomes.%s.rds', 'lambda.min.lower', analysis.name)))
    print(readRDS(sprintf('data/%s.odds.ratios.nocs.%s.rds', 'lambda.min.lower',analysis.name)))
}

for ( analysis.name in c('raw', 'adj', 'proximal' ) ) {
    print(analysis.name)
    print(readRDS(sprintf('data/%s.hazard.differences.outcomes.%s.rds', 'lambda.min.lowest', analysis.name)))
    print(readRDS(sprintf('data/%s.odds.ratios.nocs.%s.rds', 'lambda.min.lowest',analysis.name)))
}

