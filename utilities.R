library(ggtext)
# Codes
sbrt.icds  <-  c('9230', '9231', '9232', '9233', '9239',
                    'DB22DZ', 'DB22HZZ', 'DB22JZZ')
sublobar.icds  <-  c(  '3230', '3239', '3220', '3229', '0BBC4ZX', '0BBC4ZZ', '0BBC0ZX', '0BBC0ZZ', '0BBD4ZX', '0BBD4ZZ', '0BBD0ZX', '0BBD0ZZ', '0BBF4ZX', '0BBF4ZZ', '0BBF0ZX', '0BBF0ZZ', '0BBG4ZX', '0BBG4ZZ', '0BBG0ZX', '0BBG0ZZ', '0BBH4ZX', '0BBH4ZZ', '0BBH0ZX', '0BBH0ZZ', '0BBJ4ZX', '0BBJ4ZZ', '0BBJ0ZX', '0BBJ0ZZ', '0BBK4ZX', '0BBK4ZZ', '0BBK0ZX', '0BBK0ZZ', '0BBL4ZX', '0BBL4ZZ', '0BBL0ZX', '0BBL0ZZ', '0BBM4ZX', '0BBM4ZZ', '0BBM0ZX', '0BBM0ZZ')
other.resection.icds  <-  c( '3240', '3241', '3249', '3260', '3250', '3259', '0BTC4ZX', '0BTC4ZZ', '0BTC0ZX', '0BTC0ZZ', '0BTD4ZX', '0BTD4ZZ', '0BTD0ZX', '0BTD0ZZ', '0BTF4ZX', '0BTF4ZZ', '0BTF0ZX', '0BTF0ZZ', '0BTG4ZX', '0BTG4ZZ', '0BTG0ZX', '0BTG0ZZ', '0BTH4ZX', '0BTH4ZZ', '0BTH0ZX', '0BTH0ZZ', '0BTJ4ZX', '0BTJ4ZZ', '0BTJ0ZX', '0BTJ0ZZ', '0BTK4ZX', '0BTK4ZZ', '0BTK0ZX', '0BTK0ZZ', '0BTL4ZX', '0BTL4ZZ', '0BTL0ZX', '0BTL0ZZ', '0BTM4ZX', '0BTM4ZZ', '0BTM0ZX', '0BTM0ZZ')

#https://cdn.jamanetwork.com/ama/content_public/journal/surg/931810/soi140036supp1_prod.pdf?Expires=1684345221&Signature=OSPftfHA9rDQPKLJrRnmTuVgeN80y9G1jEx1SFmWPu8xdjKv-GfJLKrub~yKR12CYlUhXH6JMBap3PqTII4naX9ZndV10ahKSFYdWMciw4diKk9daxaK0V~4m7pTmaMw~a7H6H-Xqcke5-g~XeFlSz18EX08a-d5MQsC5Os0Um9v8Ng~y7-xB8t2tJygyztL2UqJqhdejzcXdUWCERoQ~fWXJJh640Qer8NojGa8Va2YT8mEXhu3q3GEsVUJ~dEooFq5X17r5w4~qzSfIylot3LJVVGxo9eslbLjfnumtjUIlbdNlLx8zoguvAjP8VoHUEpSQI88n~wOwPpKht68Bg <- &Key-Pair-Id=APKAIE5G5CRDK6RD3PGA


sbrt.cpts  <-  c('77373', 'G0173', 'G0251', 'G0339', 'G0340', '61793',  '0082T' )


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
    xlims <- c(0.5, 4)
    tt  <-  ifelse (hazard, 'Hazard Ratio (log scale)', 'Incidence Rate Ratio (log scale)')
    g <- ggplot(odds.ratios_, aes(x = estimate, y=y_axis)) + 
        geom_vline(aes(xintercept = 1), size = 0.25, linetype = "dashed") +
        geom_errorbarh(aes( xmax = high_ci, xmin = low_ci), size = 0.20, height = 0.3)+
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
              panel.grid.major.y = element_line( size=.05, color="grey", linetype = 'dashed' ),
        legend.position = 'bottom',
      #  legend.title = element_blank(),
        legend.key.size=grid::unit(2,"lines"),
         axis.text.y = ggtext::element_markdown()) 
 g
}

make.HD.plot  <-  function (odds.ratios_, label_list2, hazard =F) {
    xlims <- c(-0.05,0.25)
    tt  <-  'Hazard Difference'
    g <- ggplot(odds.ratios_, aes(x = estimate, y=y_axis)) + 
        geom_vline(aes(xintercept = 0), size = 0.25, linetype = "dashed") +
        geom_errorbarh(aes( xmax = high_ci, xmin = low_ci), size = 0.20, height = 0.3)+
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
              panel.grid.major.y = element_line( size=.05, color="grey", linetype = 'dashed' ),
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
