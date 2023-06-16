library(ggtext)
# Codes

find.rows.idx  <- function( haystack, needles) {
        found.rows <- ((sapply(haystack, `%in%`, needles)))
        found.rows  <- found.rows  %>% as_tibble %>% mutate(found.idx = ifelse (rowSums(.) ==0, NA,    max.col(., ties.method='first'))) %>% select(found.idx) %>% unlist
    return(found.rows)
}

find.rows  <- function( haystack, needles) {
        found.rows <- rowSums(sapply(haystack, `%in%`, needles)) > 0
    return(found.rows)
}

find.rows.icdsmart  <- function( haystack, needles, icd9or10) {
        found.rows.icd9 <- rowSums(sapply(haystack, `%in%`, needles[['icd9']])) > 0
        found.rows.icd10 <- rowSums(sapply(haystack, `%in%`, needles[['icd10']])) > 0
        found.rows = (found.rows.icd9 & icd9or10 == 'icd9')  | (found.rows.icd10 & icd9or10 == 'icd10')
    return(found.rows)
}


get.dates.of.procedure  <-  function( A, proc.codes ) {
    date.cols  <-  A %>% select( SRGCL_PRCDR_PRFRM_1_DT:SRGCL_PRCDR_PRFRM_25_DT) %>% colnames
    # Date of each procedure 
    # First, obtain the index of each column that has the procedure of interest
    A$proc.colidx  <- find.rows.idx( A %>% select( SRGCL_PRCDR_1_CD:SRGCL_PRCDR_25_CD ), proc.codes )
    # Now that the column index is obtained, need to grab the actual date
    A <-  A %>% 
        rowwise() %>% 
        mutate( proc.date = 
               ifelse ( is.na(proc.colidx), NA, unlist(cur_data()[ date.cols[proc.colidx]] ) )
        )
    return(unlist(A$proc.date))
}


get.dates.of.dx  <-  function( A, proc.codes ) {
    date.cols  <-  A %>% select( SRGCL_PRCDR_PRFRM_1_DT:SRGCL_PRCDR_PRFRM_25_DT) %>% colnames
    # Date of each procedure 
    # First, obtain the index of each column that has the procedure of interest
    A$proc.colidx  <- find.rows.idx( A %>% select( SRGCL_PRCDR_1_CD:SRGCL_PRCDR_25_CD ), proc.codes )
    # Now that the column index is obtained, need to grab the actual date
    A <-  A %>% 
        rowwise() %>% 
        mutate( proc.date = 
               ifelse ( is.na(proc.colidx), NA, unlist(cur_data()[ date.cols[proc.colidx]] ) )
        )
    return(unlist(A$proc.date))
}


make.OR.plot  <-  function (odds.ratios_, label_list2, hazard =F) {
    xlims <- c(0.3, 4)
    tt  <-  ifelse (hazard, 'Hazard Ratio (log scale)', 'Hazard Ratio (log scale)')
    g <- ggplot(odds.ratios_, aes(x = estimate, y=y_axis)) + 
        geom_vline(aes(xintercept = 1), linewidth = 0.25, linetype = "dashed") +
        geom_errorbarh(aes( xmax = high_ci, xmin = low_ci), linewidth = 0.20, height = 0.3)+
        geom_point(size=1.5) +
        theme_bw() +
        theme(panel.grid.minor = element_blank()) +
        #scale_y_continuous(breaks = 1:max(odds.ratio$y_axis), labels = unlist(label_list[rownames(odds.ratios__propensity)]), trans='reverse') +
        scale_y_continuous(breaks = 1:max(odds.ratios_$y_axis), labels = label_list2[row.names(odds.ratios_)], trans='reverse') +
        #scale_x_continuous(breaks = seq(0,1.4,0.2), limits = xlims ) +
        scale_x_continuous(limits = xlims ) +
        coord_trans(x = "log10") +
        xlab(tt) +
        ylab("") +
        scale_linetype_manual(values=c("solid","dashed"))+
        scale_shape_manual(values=c(15,17))+
        theme( panel.grid.major.x = element_blank() ,
              panel.grid.major.y = element_line( linewidth=.05, color="grey", linetype = 'dashed' ),
        legend.position = 'bottom',
      #  legend.title = element_blank(),
        legend.key.size=grid::unit(2,"lines"),
         axis.text.y = ggtext::element_markdown()) 
 g
}

make.HD.plot  <-  function (odds.ratios_, label_list2, hazard =F, xlims = c(-0.05,0.25)) {
    tt  <-  'Hazard Difference'
    g <- ggplot(odds.ratios_, aes(x = estimate, y=y_axis)) + 
        geom_vline(aes(xintercept = 0), linewidth = 0.25, linetype = "dashed") +
        geom_errorbarh(aes( xmax = high_ci, xmin = low_ci), linewidth = 0.20, height = 0.3)+
        geom_point(size=1.5) +
        theme_bw() +
        theme(panel.grid.minor = element_blank()) +
        scale_y_continuous(breaks = 1:max(odds.ratios_$y_axis), labels = label_list2[row.names(odds.ratios_)], trans='reverse') +
        scale_x_continuous(limits = xlims ) +
        #coord_trans(x = "log10") +
        xlab(tt) +
        ylab("") +
        scale_linetype_manual(values=c("solid","dashed"))+
        scale_shape_manual(values=c(15,17))+
        theme( panel.grid.major.x = element_blank() ,
              panel.grid.major.y = element_line( linewidth=.05, color="grey", linetype = 'dashed' ),
        legend.position = 'bottom',
      #  legend.title = element_blank(),
        legend.key.size=grid::unit(2,"lines"),
         axis.text.y = ggtext::element_markdown())
 g
}

make.odds.ratio.df  <-  function(outcome.names ) {
    odds.ratios <- as.data.frame(matrix(NA, ncol=3,nrow=length(outcome.names)))
    rownames(odds.ratios) <- c(outcome.names)
    colnames(odds.ratios) <- c('estimate','low_ci', 'high_ci')
    odds.ratios$y_axis <- 1:nrow(odds.ratios)
    odds.ratios$outcome <- rownames(odds.ratios)
    odds.ratios
}

explain_icd9_10  <-  function (dx) {
    dx.string = explain_table(as.icd10(dx), condense=F)$short_desc
    dx.string = ifelse (is.na(dx.string) ,  explain_table(as.icd9(dx), condense=F)$short_desc, dx.string)
    return(dx.string)
}

nna  <-  function ( x) !is.na(x)
################################
# Backup 
################################
# system.time ( MBSF_2013 <- read.SAScii('../SEER-Medicare-data/data/SEER_Medicare/mbsf.abcd.summary.2013.txt', 
#             '../SEER-Medicare-data/data/SEER_Medicare/2020 Input Statements/MBSF.ABCD.gcl.txt', n=100000) ) 

# lung.SEER <- read.SAScii('../SEER-Medicare-data/data/SEER_Medicare/SEER.lung.cancer.n1000.txt',  
                         # '../SEER-Medicare-data/data/SEER_Medicare/2020 Input Statements/SEER.input.gcl.txt', n=100) 
