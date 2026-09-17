library(ggplot2)
library(tidyverse)
library(patchwork)
source('utilities.R')

analysis.name  <- 'data49.nctime90_all.gte.65'
outcome.name  <- 'death.lc.specific'
res  <- readRDS(sprintf('data/%s.%s.martingale.residuals.rds', analysis.name, outcome.name))
# use ggplot

plot_list <- list()
label_list2  <- c(label_list,
                  'sexMale_bool' = 'Male sex',
                  'A_bool' = 'SBRT',
                  'race2White_bool' = 'White race',
                  list(
                       'treatment.year22012_bool' = 'Treatment Year 2012',
                       'treatment.year22013_bool' = 'Treatment Year 2013',
                       'treatment.year22014_bool' = 'Treatment Year 2014',
                       'treatment.year22015_bool' = 'Treatment Year 2015',
                       'treatment.year22016_bool' = 'Treatment Year 2016',
                       'treatment.year22017_bool' = 'Treatment Year 2017',
                       'treatment.year22018_bool' = 'Treatment Year 2018',
                       'treatment.year22019_2020_bool' = 'Treatment Year 2019-2020'))


toplot  <- data.frame(x_ = res$hd_primary, martingale_residual_primary = res$martingale_residual_primary)
#Exclude the x values greater than 99.5th percentile and 0.5th percentile
xlm0  <- quantile(toplot$x_, 0.01)
xlm1  <- quantile(toplot$x_, 0.98)
g1 <- ggplot(toplot, aes(x = x_, y = martingale_residual_primary)) +
    # geom_point(alpha = 0.4, size = 1) +
   geom_smooth(method = "loess", se = TRUE, color = "blue", level = 0.95) +
    # geom_smooth(method='loess') +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
     labs(x = 'Fitted values', y = "Lung-cancer mortality residuals") +
     xlim(xlm0, xlm1) +
    theme_minimal()
toplot  <- data.frame(x_ = res$hd_cprisk, martingale_residual_primary = res$martingale_residual_cprisk)
xlm0  <- quantile(toplot$x_, 0.01)
xlm1  <- quantile(toplot$x_, 0.98)
g2 <- ggplot(toplot, aes(x = x_, y = martingale_residual_primary)) +
    # geom_point(alpha = 0.4, size = 1) +
   geom_smooth(method = "loess", se = TRUE, color = "blue", level = 0.95) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
     labs(x = 'Fitted values', y = "Other-cause mortality residuals") +
     xlim(xlm0, xlm1) +
    theme_minimal()
gg  <-  g1 + g2 + plot_layout(ncol = 1)
gg
ggsave(gg, filename = sprintf('figs/martingales3.pdf'), width = 6, height = 10)
