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
medpars


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
carrier$sbrt.date  <-  ifelse (carrier$sbrt, carriers$CLM_THRU_DT,NA_character_) %>% ymd 
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
table(carrieri.small$sbrt , useNA="ifany")

outpati  <-   read_dta('../SEER-Medicare-data/data/SEER_Medicare/outpat2014.base.dta')
outpati.small  <- outpati %>% inner_join(lung.SEER.pids) 
fofo  <-  outpati %>% filter (ICD_PRCDR_CD1 !='')

    medpars[[i]]  <-  medpari %>% inner_join(lung.SEER.pids) %>% select( PATIENT_ID, ADMSN_DT,  DSCHRG_DT, SRGCL_PRCDR_IND_SW, DGNS_1_CD:DGNS_25_CD, SRGCL_PRCDR_1_CD:SRGCL_PRCDR_25_CD, SRGCL_PRCDR_PRFRM_1_DT:SRGCL_PRCDR_PRFRM_25_DT)




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

