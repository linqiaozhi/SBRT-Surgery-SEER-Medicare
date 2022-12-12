install.packages('https://cran.r-project.org/src/contrib/Archive/SAScii/SAScii_1.0.tar.gz', repos=NULL, type="source")
library('SAScii')
library(lubridate)
library('tidyverse')
library(haven)
library(icd)#devtools::install_github("jackwasey/icd")
source('utilities.R')
data.path  <- '../SEER-Medicare-data/data/SEER_Medicare'


################################
# Plan: 
# Medpar: provides the surgery information. Included SBRT in case the patient
# received this as an inpatient .
#  Carrier lines: One line for each procedure. This is where we obtain the SBRT, which are all done as outpatient procedures. 
# Carrier base, outpatient, and DME: for diagnoses
################################


################################
# Load SEER 
################################

lung.SEER <- read_dta('../SEER-Medicare-data/data/SEER_Medicare/SEER.lung.cancer.dta')
lung.SEER.pids  <-  lung.SEER %>% select(PATIENT_ID) %>% distinct(PATIENT_ID) 

################################
# Process the Medpar files
################################
# unlink(fn.RDS) # to start from scratch
fn.RDS  <- sprintf('%s/medpar.RDS', data.path)
if ( ! file.exists (fn.RDS) ) {
    medpars  <-  list()
    years  <-  as.character(2009:2019)
    for (yeari in 1:length(years)) {
        year  <-  years[yeari]
        print(year)
        medpari  <-   read_dta(sprintf('%s/medpar%s.dta', data.path, year))
        # inner join with the SEER patients to reduce size
        medpars[[year]]  <-  medpari %>% inner_join(lung.SEER.pids) %>% select( PATIENT_ID, ADMSN_DT,  DSCHRG_DT, SRGCL_PRCDR_IND_SW, DGNS_1_CD:DGNS_25_CD, SRGCL_PRCDR_1_CD:SRGCL_PRCDR_25_CD, SRGCL_PRCDR_PRFRM_1_DT:SRGCL_PRCDR_PRFRM_25_DT)
        medpars[[year]]$sbrt.date  <-  get.dates.of.procedure( medpars[[year]], sbrt.icds  )
        medpars[[year]]$sublobar.date  <-  get.dates.of.procedure( medpars[[year]], sublobar.icds  )
        medpars[[year]]  <-   medpars[[year]]  %>% filter( !is.na(sbrt.date) | !is.na(sublobar.date)  )
    }
    medpar  <-  bind_rows(medpars ,  .id='dataset.year')
    saveRDS(object = medpar, file = fn.RDS) 
}else{
    medpar  <-  readRDS(fn.RDS)
}
medpar <- medpar %>% mutate( sbrt.date = ymd(sbrt.date), sublobar.date = ymd(sublobar.date),)


################################
# Process the Carrier line files 
################################
years  <-  as.character(2010:2019)
carriers  <-  list()
for (yeari in 1:length(years)) {
    year  <-  years[yeari]
    fn  <-  sprintf('../SEER-Medicare-data/data/SEER_Medicare/nch%s.line.RDS', year)
    if (!file.exists( fn )) {
        dta.fn  <-  sprintf('../SEER-Medicare-data/data/SEER_Medicare/nch%s.line.dta', year )
        print(sprintf('Reading in %s', dta.fn))
        carrieri  <-   read_dta(dta.fn)
        carrieri.small  <- carrieri %>% inner_join(lung.SEER.pids)  %>% select( PATIENT_ID, CLM_THRU_DT, HCPCS_CD, LINE_ICD_DGNS_CD)
        saveRDS(object = carrieri.small,file = fn )
        unlink(dta.fn)
    }else {
        print(sprintf('Reading in %s', fn))
        carrieri.small  <- readRDS(fn)
    }
    carriers[[year]]  <-  carrieri.small
}
carrier  <-  bind_rows(carriers,  .id='dataset.year')
rm(carriers); gc();

carrier$sbrt  <-  find.rows( carrier %>% select( HCPCS_CD) , sbrt.cpts )
carrier$sbrt.date  <-  ifelse ( carrier$sbrt, carrier$CLM_THRU_DT,NA_character_) %>% ymd 
table( carrier$sbrt, useNA="ifany") # 33,000 SBRTs

#TODO: Check the diganosis codes--make sure it's SBRT for lung cancer

# fofo  <-  find.rows( carrier %>% select( HCPCS_CD) , sublobar.icds )
# table( fofo, useNA="ifany") # 0
# No resections in carrier files

################################
# Process the DME line files 
################################
# year  <- 2016
# dmei  <-   read_dta(sprintf('%s/dme%s.line.dta', data.path, year))
# dmei$sbrt  <-  find.rows( dmei %>% select( HCPCS_CD) , sbrt.cpts )
# table(dmei$sbrt , useNA="ifany")
#dmei$sublobar  <-  find.rows( dmei %>% select( HCPCS_CD) , sublobar.icds )
#table( dmei$sublobar, useNA="ifany") # 0
# No SBRT or sublobar resections in in DME files, will not load

################################
# Process the  Outpatient files
################################
outpati  <-   read_dta('../SEER-Medicare-data/data/SEER_Medicare/outpat2014.base.dta')
outpati.small  <- outpati %>% inner_join(lung.SEER.pids) 
outpati.small[  outpati.small == "" ]  <-  NA



object.size(outpati.small)/1E9
fofo





fn.RDS  <- sprintf('%s/outpat.base.RDS', data.path)
if ( ! file.exists (fn.RDS) ) {
    outpats  <-  list()
    years  <-  as.character(2010:2019)
    for (yeari in 1:length(years)) {
        year  <-  years[yeari]
        print(year)
        outpati  <-   read_dta(sprintf('%s/outpat%s.base.dta', data.path, year))
        # inner join with the SEER patients to reduce size
        outpats[[year]]  <-  outpati %>% 
            inner_join(lung.SEER.pids) %>% 
            select( PATIENT_ID, CLM_FROM_DT, CLM_THRU_DT, PRNCPAL_DGNS_CD:PRCDR_DT25  ) %>% select( ! contains( "PRCDR_DT"))
    }
    outpat  <-  bind_rows(outpats,  .id='dataset.year')
    saveRDS(object = outpat, file = fn.RDS) 
}else{
    outpat  <-  readRDS(fn.RDS)
}
outpat  <-  outpat %>% mutate( CLM_FROM_DT = ymd(CLM_FROM_DT), CLM_THRU_DT = ymd(CLM_THRU_DT))



################################
# SEER variables 
################################

nna  <-  function ( x) !is.na(x)

lung.SEER

label_list  <-  list(  AGERECODEWITHSINGLEAGES_AND_100 = 'Age')

table( lung.SEER$RACE_RECODE_WHITE_BLACK_OTHER, useNA="ifany")
table( lung.SEER$YEAR_OF_DIAGNOSIS, useNA="ifany")

# SEQUENCE_NUMBER, extent of disease width, SURVIVAL_MONTHS, TNM
A  <-  lung.SEER %>% select( PATIENT_ID, AGERECODEWITHSINGLEAGES_AND_100, 
                            MARITAL_STATUS_AT_DIAGNOSIS, 
                            RACE_RECODE_WHITE_BLACK_OTHER, 
                            SEX,
                            MONTH_OF_DIAGNOSIS, YEAR_OF_DIAGNOSIS,
                            SEER_DATEOFDEATH_MONTH, SEER_DATEOFDEATH_YEAR,
                            HISTOLOGIC_TYPE_ICD_O_3,
                            SEERCAUSESPECIFICDEATHCLASSIFIC,
                            SEEROTHERCAUSEOFDEATHCLASSIFICA) %>% 
                    mutate ( 
                            dx.date = ymd( ifelse ( nna(YEAR_OF_DIAGNOSIS) & nna(MONTH_OF_DIAGNOSIS) , sprintf('%d%02d15', YEAR_OF_DIAGNOSIS, MONTH_OF_DIAGNOSIS), NA_character_ ) )  ,
                            death.date = ymd( ifelse ( ""!=(SEER_DATEOFDEATH_YEAR) & ""!=(SEER_DATEOFDEATH_MONTH) , sprintf('%s%s15', SEER_DATEOFDEATH_YEAR, SEER_DATEOFDEATH_MONTH), NA_character_ ) )  
                            ) %>% arrange(dx.date) %>% distinct(PATIENT_ID, .keep_all =T)


                            A %>% print(n=100, width=Inf)



#TODO: Figure out why there are duplicated entries and be smarter about picking one

################################
#  Merge SEER with the treatment status
################################
# Determine the treatment. Sublobar date comes from Medpar, and SBRT date comes
# from both Medpar and carrier. We bind them into one first, and then merge.
medpar.tx  <-   medpar %>% select ( PATIENT_ID, sbrt.date, sublobar.date) %>% filter(!is.na(sbrt.date) | !is.na(sublobar.date) )
carrier.tx  <-  carrier %>% select(PATIENT_ID, sbrt.date) %>% filter(!is.na(sbrt.date))
#medpar.carrier.tx  <-  bind_rows ( medpar.tx, carrier.tx)

# Running these summarize statements on all hundreds of thousands of SEER patients does not make sense, but we 
# do need the dates from medpar. Will restrict to medpar.carrier.tx patients then join back to SEER
medpar.carrier.tx  <- bind_rows ( medpar.tx, carrier.tx)  %>% 
    left_join(A )  %>%
    group_by( PATIENT_ID ) %>% 
    summarise ( 
               dx.date = first(dx.date), 
               tx = factor( case_when ( 
                # any( nna( sbrt.date) ) & ! any( nna(sublobar.date))  ~ first(na.omit(sbrt.date)),
                # ! any( nna( sbrt.date) ) &  any( nna(sublobar.date)) ~ first(na.omit(sublobar.date)),
                any( nna( sbrt.date) ) & ! any( nna(sublobar.date))  ~ 'sbrt',
                ! any( nna( sbrt.date) ) &  any( nna(sublobar.date)) ~ 'sublobar',
                T ~ (NA_character_)
                ), levels = c('sublobar', 'sbrt')),
               tx.date = case_when (
                                    tx == 'sbrt' ~ first(sbrt.date),
                                    tx == 'sublobar' ~ first(sublobar.date),
                                    T ~ ymd(NA_character_)
                                    ) ) %>%
               mutate( tx.after.dx = tx.date > dx.date)
table( medpar.carrier.tx$tx.after.dx, useNA="ifany")
table( medpar.carrier.tx$tx, medpar.carrier.tx$tx.after.dx, useNA="ifany")
medpar.carrier.tx.2 <- medpar.carrier.tx %>% filter( tx.after.dx & nna(tx.date) ) 

A2 <- A %>% right_join(medpar.carrier.tx.2)  %>% mutate (
                                 death = death.date,
                                 tt = as.numeric( if_else ( nna(death.date), death.date, ymd('20191231')  ) - tx.date, units = 'days')
                                 ) %>% 
                     filter( tt >0 )
table( A2$death, useNA="ifany")

                 #TODO: Include tt = 0
sum(A2$tt == 0 )


 # max(ymd(carrier$CLM_THRU_DT)) # 2019-12-32

################################
# Merge with negative outcomes 
################################
# Create negative outcomes
# Get common diagnoses
outpati.long  <-  outpats[['2017']] %>% 
    select( PATIENT_ID, ICD_DGNS_CD1:ICD_DGNS_E_CD12) %>% 
    pivot_longer(!PATIENT_ID, 
                 names_to = 'diagnosis.field', values_to = 'diagnosis.code', 
                 values_drop_na = T) 
conditions  <-  outpati.long %>% 
    filter ( diagnosis.code != "") %>% 
    count(diagnosis.code) %>% arrange(-n)%>% top_n ( 1000) %>% 
    mutate( explanation= explain_table(as.icd10( diagnosis.code),condense=F)$short_desc)
conditions %>% print(n=Inf)

# Create some common lists
negative.outcomes  <-  list(
        'fall' = c( expand_range('E880','E888'),expand_range('W00', 'W19' )),
        'acute_bronchitis' = c('4660', '4661', sprintf('J20%d',0:9))
)


# There are two kinds of negative outcomes: 1) Any outcome that occured after the treatment and 2) Any outcome that occured for the first time after the treatment.  Let's focus on the first for now. 

A3  <-  A2
for (i in 1:length(negative.outcomes) ) {
    noc.name  <- names(negative.outcomes)[i]
    noc.codes  <- negative.outcomes[[noc.name]]
    outpat.noc <- outpat %>% right_join(A2)  %>% 
         mutate( 
                noc.temp = find.rows( across(ICD_DGNS_CD1:ICD_DGNS_E_CD12), noc.codes),
                noc.temp =  if_else(noc.temp & (CLM_FROM_DT > tx.date)& (CLM_FROM_DT < death.date), CLM_FROM_DT, ymd(NA_character_)) ) %>% 
         group_by(PATIENT_ID) %>% 
         summarise( !!noc.name := first(na.omit(noc.temp))) %>% 
         filter(nna(!!rlang::sym(noc.name)))
     A3 <- A3 %>% left_join(outpat.noc)
}

#table( as.numeric(A3$fall - A3$tx.date, units = 'days') > A3$tt, useNA="ifany")


outcome.name  <-  'fall'
A4  <-  A3 %>% mutate( 
                      # Contributing person time until 1) time of outcome, 2) death, or 3) end of follow up. The latter two are in tt
                      outcome.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                      outcome.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F))
#A4  <-  within(A4, tx  <-  relevel(tx, ref = 'sublobar'))
m  <- glm( (outcome.bool) ~ tx + AGERECODEWITHSINGLEAGES_AND_100+   offset( log(outcome.time) ) , data = A4, family = poisson(link=log))
summary(m)
exp(c( coef(m)['txsbrt'], confint(m,'txsbrt'))) # http://www.philender.com/courses/categorical/notes1/pois1.html





#TODO: Include the inclusion year
m  <- glm( (SEERCAUSESPECIFICDEATHCLASSIFIC) ~ tx + offset( log(tt) ) + AGERECODEWITHSINGLEAGES_AND_100 , data = A4, family = poisson(link=log))
summary(m)
exp(coef(m))



A3 %>% 
    mutate(across(fall:acute_bronchitis, nna))%>%
    group_by(tx) %>% 
    summarise( across( 
                  c('AGERECODEWITHSINGLEAGES_AND_100', 'SEERCAUSESPECIFICDEATHCLASSIFIC', 'fall', 'acute_bronchitis') , ~ mean( .x, na.rm = T))) %>% glimpse

summary(m)


# at this point any event is valid, and we just want if it occurred
outpat.noc.out <- outpat.noc %>% outpat.noc.out


table( outpat.noc$fall_, useNA="ifany")

with( outpat.noc, table( nna(fall), (fall_), useNA="ifany"))
glimpse(outpat.noc %>% filter( nna(fall)))


outpat.noc <- outpat.noc %>% mutate( fall = ifelse( 

outpat.noc$fall_date
table(outpat.noc$fall_date , useNA="ifany")

%>% select( ICD_DGNS_CD1:ICD_DGNS_E_CD12) %>% find.rows (negative.outcomes[['fall']])

print(negative.outcomes)

expand_range('J200', 'J209')

table( outpat$dataset.year, useNA="ifany")


outpat.noc %>% mutate( fall)
table( outpat$fall, useNA="ifany")

sum(outpat$CLM_FROM_DT == "")

#TODO: Do I need to do survival analysis for each negative outcome? The length of follow up is so different between the two groups...
# Need to use person-time
# Proportional hazards does not apply here. You have ahug















    table( (fofo$tx.date  > fofo$dx.date), useNA="ifany")
    summary( interval(fofo$tx.date  , fofo$dx.date)%/%months(1), useNA="ifany")
    table( , useNA="ifany")

glimpse(fofo)
    

A2  <- A  %>% left_join(medpar.carrier.tx )  %>% mutate( tx = case_when ( nna( sbrt.date) & is.na( sublobar.date) &  


# 







#TODO: Is there an issue where the bulk of SBRT patients were treated more
#recently so have less time to die? 

colnames(lung.SEER)
l
                            


)


SEERcausespecificdeathclassific
SEERothercauseofdeathclassifica


# Complications








# SBRT 92.3, 92.30-92.39       77373, G0173, G0251, G0339, G0340, 61793,  0082T
#3230 Thorac seg lung resect
#3239 Oth seg lung resect NOS
#3220 Thorac exc lung lesion
#3229 Destroy loc lung les NEC





medpar2016.small$lobectomy.date  <-  get.dates.of.procedure( medpar2016.small, c('3241', '3249' )  )
# https://www.jtcvsopen.org/article/S2666-2736(22)00316-3/pdf

#table(is.na(medpar2012.small$sublobar.date) , useNA="ifany")
#1745 Thoraco robotic ast proc
medpar2012.small$rats.date  <-  get.dates.of.procedure( medpar2012.small, c('1745')  )
#3220 Thorac exc lung lesion
#3421 Transpleura thoracoscopy
#3241 Thorac lobectomy lung
#3230 Thorac seg lung resect
#3250 Thoracospc pneumonectomy
medpar2012.small$vats.date  <-  get.dates.of.procedure( medpar2012.small, c('3220','3421','3241', '3230', '3250')  )
#3229 Destroy loc lung les NEC - Not including, as it has ablation
#3402 Exploratory thoracotomy
#3249 Lobectomy of lung NEC
#3239 Oth seg lung resect NOS
#3259 Other pneumonectomy NOS
#329 Other excision of lung
medpar2012.small$open.date  <-  get.dates.of.procedure( medpar2012.small, c( '3402',  '3249' , '3239', '3259', '329')  )
medpar2012.small <- medpar2012.small %>% mutate( 
                        approach = case_when( 
                                             !is.na(vats.date) & is.na(rats.date) ~ 'VATS',
                                             !is.na(rats.date) & is.na(vats.date) ~ 'RATS',
                                             !is.na(open.date) & is.na (rats.date) & is.na(vats.date) ~ 'open',
                                             T ~ 'exclude' ),
                        surgery.date = case_when( 
                                                 approach == 'VATS' ~ vats.date,
                                                 approach == 'RATS' ~ rats.date,
                                                 approach == 'open' ~ open.date,
                                                 T ~ NA_character_
                                                 ) %>% ymd
                        ) 

fofo  <-  medpar2012.small %>% filter ( !is.na(surgery.date)) %>% group_by( PATIENT_ID) %>% arrange( surgery.date) %>% filter(row_number() == 1)
table( fofo$approach, useNA="ifany")
table( !is.na(fofo$lobectomy.date), useNA="ifany")

# KEep only the first episode
fofo  <-  medpar2012.small %>% group_by( PATIENT_ID) %>% arrange( surgery.date) %>% mutate(nn = n())
fifi  <-  fofo %>% filter (nn > 1)
sum(fifi[1,] != fifi[2,], na.rm = T)
lung.SEER$Histology_ICD_O_2
table( lung.SEER$HISTOLOGIC_TYPE_ICD_O_3, useNA="ifany")

lung.SEER %>% left_join( medpar2012.small %>% select( PATIENT_ID, approach, surgery.date)) 

################################
#  
################################

sum(medpar2012.small$exclude)

with( medpar2012.small, table( !is.na(lobectomy.date), !is.na(segment.date), useNA="ifany"))
with( medpar2012.small, table( !is.na(vats.date), !is.na(open.date), useNA="ifany"))

table(medpar2012.small$lobectomy.date , useNA="ifany")


# Date of each lobectomy 
# First, obtain the index of each column that has a lobectomy
medpar2012.small$lobectomy.colidx  <- find.rows.idx( medpar2012.small %>% select( SRGCL_PRCDR_1_CD:SRGCL_PRCDR_25_CD ), c('3249', '3241') )
# Now that the column index is obtained, need to grab the actual date
medpar2012.small <-  medpar2012.small %>% 
    select(lobectomy.colidx,  SRGCL_PRCDR_PRFRM_1_DT:SRGCL_PRCDR_PRFRM_25_DT ) %>% 
    rowwise() %>% 
    mutate( lobectomy.date = 
           ifelse ( is.na(lobectomy.colidx), NA, unlist(cur_data()[ date.cols[lobectomy.colidx]] ) )
    )






first(which(medpar2012.small$lobectomy.colidx > 1))

glimpse(medpar2012.small[348,])





table( medpar2012.small$lobectomy.colidx, useNA="ifany")



fofo <-  medpar2012.small %>% select(lobectomy.colidx,  SRGCL_PRCDR_PRFRM_1_DT:SRGCL_PRCDR_PRFRM_25_DT ) %>% rowwise() %>% mutate( date = cur_data()[match(date.cols[lobectomy.colidx], names(.))    ])


medpar2012.small
fofo %>%  %>% select(test)

head(fofo$max.col2)

medpar2012.small

#medpar2012.small %>% head(10) %>% mutate(code_description= explain_code(SRGCL_PRCDR_1_CD,condense=F )) %>% select(code_description)
?explain_code
medpar2012.small %>% print(width=Inf)
table( medpar2012$SRGCL_PRCDR_IND_SW, useNA="ifany")
class(medpar2012$SRGCL_PRCDR_IND_SW)
# Check if surgery was present
explain_code('0393')

fofo  <- find.rows( medpar2012 %>% select( SRGCL_PRCDR_1_CD:SRGCL_PRCDR_25_CD ), c('3249', '3241') )
table( fofo, useNA="ifany")
head(fofo)




negative.outcomes.map  <- list ( 
                                fracture = er('80000', '8291'),
                                fecal.incontinence = er('78760', '78763' ) )


seer.outpat  <-  lung.SEER[1:10000,] %>% inner_join( outpat2014.base) %>%  select ( PATIENT_ID, CLM_FROM_DT, ICD_DGNS_CD1:ICD_DGNS_CD24)
seer.outpat <- seer.outpat %>% unite("ID_DATE", PATIENT_ID:CLM_FROM_DT)

colnames(lung.SEER)
table( lung.SEER$SITESPECIFICSURGERY19731997VARY , useNA="ifany")

seer.outpat[2,1] = seer.outpat[1,1]
seer.outpat[1,3] = '78763'


seer.outpat.no  <-  comorbid ( seer.outpat, negative.outcomes.map, restore_id_order = T, return_df = T) %>% as_tibble 

seer.outpat[duplicated(seer.outpat$ID_DATE),]
seer.outpat.no[duplicated(seer.outpat.no$ID_DATE),]

%>% separate ("ID_DATE", c("PATIENT_ID", "CLM_FROM_DT"))

which(!(seer.outpat$ID_DATE %in% seer.outpat.no$ID_DATE) ) 


table( seer.outpat.no$fracture, useNA="ifany")
table( seer.outpat.no$fecal.incontinence, useNA="ifany")
colnames(outpat2014.base)
colnames(lung.SEER)

outpat2014.base %>% print(n = 5, width = Inf)




diagnoses  <-  outpat2014.base[1:500000,] %>% select ( ICD_DGNS_CD1:ICD_DGNS_CD24) %>% unlist

er  <-  icd::expand_range






fofo

outpat2014.base <- read_dta('../SEER-Medicare-data/data/SEER_Medicare/outpat2014.base.dta')
colnames(outpat2014.base)
diagnoses  <-  outpat2014.base[1:500000,] %>% select ( ICD_DGNS_CD1:ICD_DGNS_CD24) %>% unlist
diagnoses.frequencies   <-  tibble( dx = diagnoses )  %>% count(dx) %>% mutate( dx.string = explain_code(dx, condense=F))  %>% arrange(-n)
diagnoses.frequencies %>% filter ( n > 1 ) %>% print (n=Inf)
fifi

length(fifi$dx[1:25] )
explain_code(fifi$dx[22] )
length((unlist(fifi$dx[21:22] )))
length(explain_code(unlist(fifi$dx )[21:22]))

unlist(fifi$dx[21:22] )
explain_code(unlist(fifi$dx )[21:22])


unlist(fifi$dx )[1:22]
fifi %>% print (n=100)

for (i in 1:32 ) {
    cat(sprintf('\n\n %d\n', i))
    print(fifi$dx[i])
    print(explain_code(fifi$dx[i]))
}







comorbid_quan_deyo( fofo[1:10,])
data(vermont_dx)
head(vermont_dx)

icd9_map_elix$CHF


system.time ( outpat_2013 <- read.SAScii('../SEER-Medicare-data/data/SEER_Medicare/outpat2013.base.txt', 
            '../SEER-Medicare-data/data/SEER_Medicare/2020 Input Statements/outpat.gcl.txt', n=100) ) 


head(lung.SEER)

colnames(lung)
lung.SEER[1,] %>% t


View(lung.SEER)
table( lung.SEER$YEAR_OF_DIAGNOSIS)

MBSF_2019
?read.SAScii

df <- tibble(x=1:10,y=1:10)

df

0BBC4ZX Upper lung lobe, right
0BBC4ZZ Upper lung lobe, right
0BBC0ZX Upper lung lobe, right
0BBC0ZZ Upper lung lobe, right

0BBD4ZX Middle lung lobe, right
0BBD4ZZ Middle lung lobe, right
0BBD0ZX Middle lung lobe, right
0BBD0ZZ Middle lung lobe, right

0BBF4ZX Lower lung lobe, right
0BBF4ZZ Lower lung lobe, right
0BBF0ZX Lower lung lobe, right
0BBF0ZZ Lower lung lobe, right

0BBG4ZX Upper lung lobe, left
0BBG4ZZ Upper lung lobe, left
0BBG0ZX Upper lung lobe, left
0BBG0ZZ Upper lung lobe, left

0BBH4ZX Lung lingula
0BBH4ZZ Lung lingula
0BBH0ZX Lung lingula
0BBH0ZZ Lung lingula

0BBJ4ZX Lower lung lobe, left
0BBJ4ZZ Lower lung lobe, left
0BBJ0ZX Lower lung lobe, left
0BBJ0ZZ Lower lung lobe, left

0BBK4ZX Lung, right
0BBK4ZZ Lung, right
0BBK0ZX Lung, right
0BBK0ZZ Lung, right

0BBL4ZX Lung, left
0BBL4ZZ Lung, left
0BBL0ZX Lung, left
0BBL0ZZ Lung, left

0BBM4ZX Lungs, bilateral
0BBM4ZZ Lungs, bilateral
0BBM0ZX Lungs, bilateral
0BBM0ZZ Lungs, bilateral

