
# Codes
sbrt.icds  <-  c('9230', '9231', '9232', '9233', '9239',
                    'DB22DZ', 'DB22HZZ', 'DB22JZZ')
sublobar.icds  <-  c(  '3230', '3239', '3220', '3229', '0BBC4ZX', '0BBC4ZZ', '0BBC0ZX', '0BBC0ZZ', '0BBD4ZX', '0BBD4ZZ', '0BBD0ZX', '0BBD0ZZ', '0BBF4ZX', '0BBF4ZZ', '0BBF0ZX', '0BBF0ZZ', '0BBG4ZX', '0BBG4ZZ', '0BBG0ZX', '0BBG0ZZ', '0BBH4ZX', '0BBH4ZZ', '0BBH0ZX', '0BBH0ZZ', '0BBJ4ZX', '0BBJ4ZZ', '0BBJ0ZX', '0BBJ0ZZ', '0BBK4ZX', '0BBK4ZZ', '0BBK0ZX', '0BBK0ZZ', '0BBL4ZX', '0BBL4ZZ', '0BBL0ZX', '0BBL0ZZ', '0BBM4ZX', '0BBM4ZZ', '0BBM0ZX', '0BBM0ZZ')

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





################################
# Backup 
################################
# system.time ( MBSF_2013 <- read.SAScii('../SEER-Medicare-data/data/SEER_Medicare/mbsf.abcd.summary.2013.txt', 
#             '../SEER-Medicare-data/data/SEER_Medicare/2020 Input Statements/MBSF.ABCD.gcl.txt', n=100000) ) 

# lung.SEER <- read.SAScii('../SEER-Medicare-data/data/SEER_Medicare/SEER.lung.cancer.n1000.txt',  
                         # '../SEER-Medicare-data/data/SEER_Medicare/2020 Input Statements/SEER.input.gcl.txt', n=100) 
