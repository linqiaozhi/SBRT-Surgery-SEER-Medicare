#install.packages('https://cran.r-project.org/src/contrib/Archive/SAScii/SAScii_1.0.tar.gz', repos=NULL, type="source")

library(arsenal)
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
valid.dxs  <- c( expand_range('1622','1629'), expand_range(as.icd10('C34'), as.icd10('C349')))

################################
# Load SEER 
################################

lung.SEER <- read_dta('../SEER-Medicare-data/data/SEER_Medicare/SEER.lung.cancer.dta')
lung.SEER.valid.dx  <-  lung.SEER %>% filter(PRIMARY_SITE %in% valid.dxs) #exclude rows in the dataset corresponding to NON-lung cancer diagnoses (i.e., other cancers)
lung.SEER.ordered <-  lung.SEER.valid.dx[order(lung.SEER.valid.dx$SEQUENCE_NUMBER, decreasing=FALSE),] #sort by sequence number in ascending order
lung.SEER.first.lc.dx <- lung.SEER.ordered %>% distinct(PATIENT_ID, .keep_all = TRUE) #keep the lung cancer diagnosis corresponding to the LOWEST sequence number (i.e., their first LC diagnosis)
lung.SEER.pids <- lung.SEER.ordered %>% filter(YEAR_OF_DIAGNOSIS>=2010 & YEAR_OF_DIAGNOSIS<=2017) %>% select(PATIENT_ID) #restrict to patients diagnosed from 2010-2017; final SEER patient list; gives you a list of 415,741 patients (includes patients with a first primary LC diagnosis)


################################
# Process the Medpar files
################################
year = "2015"
# unlink(fn.RDS) # to start from scratch
fn.RDS  <- sprintf('%s/medpar.RDS', data.path)
if ( ! file.exists (fn.RDS) ) {
    medpars  <-  list()
    years  <-  as.character(2009:2019)
    for (yeari in 1:length(years)) {
        year  <-  years[yeari]
        print(year)
        medpari  <-   read_dta(sprintf('%s/medpar%s.dta', data.path, year), col_select=c(PATIENT_ID, ADMSN_DT,  DSCHRG_DT, SRGCL_PRCDR_IND_SW, DGNS_1_CD:DGNS_25_CD, SRGCL_PRCDR_1_CD:SRGCL_PRCDR_25_CD, SRGCL_PRCDR_PRFRM_1_DT:SRGCL_PRCDR_PRFRM_25_DT))
        # inner join with the SEER patients to reduce size
        medpars[[year]]  <-  medpari %>% inner_join(lung.SEER.pids) 
        medpars[[year]]$sbrt.date  <-  get.dates.of.procedure( medpars[[year]], sbrt.icds  )
        medpars[[year]]$sublobar.date  <-  get.dates.of.procedure( medpars[[year]], sublobar.icds  )
        #medpars[[year]] <-   medpars[[year]]  %>% filter( !is.na(sbrt.date) | !is.na(sublobar.date)  )
    }
    medpar  <-  bind_rows(medpars ,  .id='dataset.year')
    saveRDS(object = medpar, file = fn.RDS)
}else{
    medpar  <-  readRDS(fn.RDS)
}

medpar %>% group_by(dataset.year) %>% tally() #check the number of observations per year

# The location of the SBRT is not specified, so need to filter to only patients with lung cancer
medpar <- medpar %>% mutate(actually.lung.cancer = find.rows( across(DGNS_1_CD:DGNS_25_CD), valid.dxs),  sbrt.date = ymd(sbrt.date), sublobar.date = ymd(sublobar.date) ) 
medpar$sbrt.date[ ! medpar$actually.lung.cancer ]  <- as.Date(NA_character_)

valid.dxs  <- c( expand_range('1622','1629'), expand_range(as.icd10('C34'), as.icd10('C349')))



#fofo <- medpar %>% mutate( sbrt.date = ymd(sbrt.date), sublobar.date = ymd(sublobar.date),)
#fofo  <-  medpar %>% filter (nna(sbrt.date) & PATIENT_ID %in% A2$PATIENT_ID) %>% mutate(actually.lung.cancer = find.rows( across(DGNS_1_CD:DGNS_25_CD), valid.dxs)) 
#fofo$actually.lung.cancer


################################
# Process the Carrier line files 
################################

#TODO: Rerun the SAS
years  <-  as.character(2010:2019)
carriers  <-  list()
for (yeari in 1:length(years)) {
    year  <-  years[yeari]
    fn  <-  sprintf('../SEER-Medicare-data/data/SEER_Medicare/nch%s.line.RDS', year)
    if (!file.exists( fn )) {
        dta.fn  <-  sprintf('../SEER-Medicare-data/data/SEER_Medicare/nch%s.line.dta', year )
        print(sprintf('Reading in %s', dta.fn))
        carrieri  <-   read_dta(dta.fn, col_select=c('PATIENT_ID', 'CLM_THRU_DT', 'HCPCS_CD', 'LINE_ICD_DGNS_CD'))
        carrieri.small  <- carrieri %>% inner_join(lung.SEER.pids)  
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

carrier$valid.dx  <- carrier$LINE_ICD_DGNS_CD %in% valid.dxs 
carrier$sbrt  <-  find.rows( carrier %>% select( HCPCS_CD) , sbrt.cpts ) & carrier$valid.dx
carrier$sbrt.date  <-  ifelse ( carrier$sbrt, carrier$CLM_THRU_DT,NA_character_) %>% ymd 

# table( (carrier$sbrt), useNA="ifany") # 33,000 SBRTs
# 
#     FALSE      TRUE 
# 131015967     24277 
# fofo  <-  carrier %>% filter (sbrt) 
# fofo %>% count(LINE_ICD_DGNS_CD) %>% filter( LINE_ICD_DGNS_CD %in% valid.dxs)  %>%  mutate(dx.string = explain_icd9_10 (LINE_ICD_DGNS_CD))  %>% arrange(-n) %>% print(n=Inf)





# fofo  <-  find.rows( carrier %>% select( HCPCS_CD) , sublobar.icds )
# table( fofo, useNA="ifany") # 0
# No resections in carrier files

################################
# Process the Carrier base files 
################################

years  <-  as.character(2010:2019)
carrierbases  <-  list()
for (yeari in 1:length(years)) {
  year  <-  years[yeari]
  fn  <-  sprintf('../SEER-Medicare-data/data/SEER_Medicare/nch%s.base.RDS', year)
  if (!file.exists( fn )) {
    dta.fn  <-  sprintf('../SEER-Medicare-data/data/SEER_Medicare/nch%s.base.dta', year )
    print(sprintf('Reading in %s', dta.fn))
    carrierbasei  <-   read_dta(dta.fn, col_select=c('PATIENT_ID', 'CLM_FROM_DT', 'CLM_THRU_DT', 'PRNCPAL_DGNS_CD', ICD_DGNS_CD1:ICD_DGNS_CD12)) #
    carrierbasei.small  <- carrierbasei %>% inner_join(lung.SEER.pids)  
    saveRDS(object = carrierbasei.small,file = fn )
    unlink(dta.fn)
  }else {
    print(sprintf('Reading in %s', fn))
    carrierbasei.small  <- readRDS(fn)
  }
  carrierbases[[year]]  <-  carrierbasei.small
}
carrierbase  <-  bind_rows(carrierbases,  .id='dataset.year')
rm(carrierbases); gc();


# explain_code(as.icd9(carrierbase$ICD_DGNS_CD1[1:50]))

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
# SEER variables 
################################



topography  <-  read_csv(file= './ICDO3topography.csv') %>% rename(site.topography = description) %>% mutate(PRIMARY_SITE = str_remove_all( icdo3_code, fixed(".")))

#A.gt2010 %>% count(sex)

# SEQUENCE_NUMBER, extent of disease width, SURVIVAL_MONTHS, TNM

A.gt2010  <-  lung.SEER.first.lc.dx %>% #Use lung.SEER.first.lc.dx
    filter(YEAR_OF_DIAGNOSIS >=2010) %>%  
    rename ( 
            age                    = AGERECODEWITHSINGLEAGES_AND_100,
            sex                    = SEX,
            race                    = RACE_RECODE_WHITE_BLACK_OTHER,
            marital.status          = MARITAL_STATUS_AT_DIAGNOSIS,
            cause.specific.mortality = SEERCAUSESPECIFICDEATHCLASSIFIC,
            other.cause.mortality = SEEROTHERCAUSEOFDEATHCLASSIFICA,
            seer.surgery = RX_SUMM_SURG_PRIM_SITE_1998, # https://seer.cancer.gov/archive/manuals/2021/AppendixC/Surgery_Codes_Lung_2021.pdf
            ) %>%
    left_join( topography)   %>%
    mutate( 
           sex = case_when ( 
                            sex == 1 ~ 'Male',
                            sex == 2 ~ 'Female', 
                            sex == 9 ~ 'Unknown' ),
           race = case_when ( 
                             race ==1 ~ 'White',
                             race ==2 ~ 'Black',
                             race ==3 ~ 'Other',
                             T ~ 'Unknown' ),
           marital.status = case_when ( 
                                       marital.status ==1 ~ 'Never married',
                                       marital.status ==2 ~ 'Married',
                                       marital.status %>% between( 3,6) ~ 'Other',
                                       T ~ 'Unknown'),
           cause.specific.mortality = case_when ( 
                                                 cause.specific.mortality == 0 ~ 'Alive or other death',
                                                 cause.specific.mortality ==1 ~ 'Death',
                                                 T ~ 'Unknown' ),
           other.cause.mortality = case_when ( 
                                              other.cause.mortality == 0 ~ 'Alive or cancer death',
                                              other.cause.mortality ==1 ~ 'Death',
                                              T ~ 'Unknown' ),
           histology.code = sprintf( '%d/%d',  HISTOLOGIC_TYPE_ICD_O_3, BEHAVIOR_CODE_ICD_O_3 ),
           histology = case_when ( 
                                  histology.code  == '8140/3' ~ 'Adenocarcinoma, NOS',
                                  # INvasive non-mucinous adenocarcinomas
                                  histology.code  == '8250/3' ~ 'Adenocarcinoma, lepidic',
                                  histology.code  == '8550/3' ~ 'Adenocarcinoma, acinar',
                                  histology.code  == '8250/3' ~ 'Adenocarcinoma, papillary',
                                  #histology.code  == '8265/3' ~ 'Adenocarcinoma, micropapillary', # there were not any
                                  #histology.code  == '8230/3' ~ 'Adenocarcinoma, solid',
                                  #Invasive mucinous adenocarcinoma
                                  histology.code  == '8253/3' |histology.code  == '8480/3'  ~ 'Adenocarcinoma, mucinous',
                                  histology.code  == '8255/3' ~ 'Adenocarcinoma, mixed',
                                  histology.code  %in%  c('8070/3', '8071/3', '8073/3', '8075/3', '8076/3', '8078/3') ~ 'Squamous cell carcinoma',
                                  histology.code  == '8041/3' ~ 'Small cell carcinoma',
                                  histology.code  == '8000/3' ~ 'Malignancy, unspecified',
                                  histology.code  == '8046/3' ~ 'Non-small cell carcinoma',
                                  histology.code  == '8010/3' ~ 'Carcinoma, unspecified',
                                  T ~ 'Other' ),
           primary.site = case_when ( 
                                     str_detect(icdo3_code,'^C34*') ~ 'Lung',
                                     T ~ 'Other'
                                     ),
           #Histology Categorical (based on information from -- INSERT CITATION)
           histology.cat = case_when (
             histology.code %in% c('8140/3', '8550/3',  '8551/3',  '8260/3', '8230/3', '8333/3', '8144/3', '8480/3', '8253/3', '8254/3') ~ 'Adenocarcinoma',
             histology.code %in%  c('8070/3', '8071/3', '8072/3', '8073/3', '8074/3', '8075/3', '8076/3', '8078/3', '8083/3') ~ 'Squamous Cell Carcinoma',
             histology.code %in% c('8012/3', '8013/3', '8014/3') ~ 'Large Cell Carcinoma',
             histology.code == '8560/3' ~ 'Adenosquamous Cell Carcincoma',
             histology.code %in% c('8240/3',	'8241/3',	'8242/3',	'8243/3',	'8244/3',	'8245/3', '8246/3',	'8249/3') ~ 'Carcinoid Tumor',
             histology.code %in% c('8250/3', '8251/3', '8252/3', '8255/3') ~ 'Tumors Formerly Classifed as Bronchioloalveolar Carcinoma',
             histology.code == '8046/3' ~ 'Non-small Cell Carcinoma, NOS',
             histology.code == '8041/3' ~ 'Small Cell Carcinoma',
             T ~ 'Other/Unknown'),
           
           histology.simple=case_when(
             histology.cat=='Adenocarcinoma' ~ 'Adenocarcinoma',
             histology.cat=='Squamous Cell Carcinoma' ~ 'Squamous Cell Carcinoma',
             histology.cat=='Tumors Formerly Classifed as Bronchioloalveolar Carcinoma' ~ 'Tumors Formerly Classifed as Bronchioloalveolar Carcinoma', 
             histology.cat=='Non-small Cell Carcinoma, NOS' ~ 'Non-small Cell Carcinoma, NOS',
             histology.cat=='Adenosquamous Cell Carcinoma' | histology.cat=='Large Cell Carcinoma' | histology.cat=='Carcinoid' ~ 'Other/Unknown',
             T ~ 'Other/Unknown'
           ),
           
           tnm.t = case_when ( 
                              str_detect(DERIVED_SEER_COMBINED_T_2016, '^[cp]1') | DERIVED_AJCC_T_7TH_ED_2010 %>% between (100,190) | DERIVED_AJCC_T_7TH_ED_2010 %>% between(800, 810) ~ '1',
                              str_detect(DERIVED_SEER_COMBINED_T_2016, '^[cp]2') | DERIVED_AJCC_T_7TH_ED_2010  %>% between (200,290)~ '2',
                              str_detect(DERIVED_SEER_COMBINED_T_2016, '^[cp]3')| DERIVED_AJCC_T_7TH_ED_2010  %>% between (300,390) ~ '3',
                              str_detect(DERIVED_SEER_COMBINED_T_2016, '^[cp]4')| DERIVED_AJCC_T_7TH_ED_2010  %>% between (400,499) ~ '4',
                              str_detect(DERIVED_SEER_COMBINED_T_2016, '^[cp]X') | DERIVED_AJCC_T_7TH_ED_2010  == 888 ~ 'X',
                              T ~ NA_character_ ),
           tnm.n = case_when ( 
                              str_detect(DERIVED_SEER_COMBINED_N_2016, '^[cp]0') | DERIVED_AJCC_N_7TH_ED_2010  %>% between (0,40) ~ '0',
                              str_detect(DERIVED_SEER_COMBINED_N_2016, '^[cp]1') | DERIVED_AJCC_N_7TH_ED_2010  %>% between (100,199) ~ '1',
                              str_detect(DERIVED_SEER_COMBINED_N_2016, '^[cp]2') | DERIVED_AJCC_N_7TH_ED_2010  %>% between (200,299)~ '2',
                              str_detect(DERIVED_SEER_COMBINED_N_2016, '^[cp]3') | DERIVED_AJCC_N_7TH_ED_2010  %>% between (300,399) ~ '3',
                              str_detect(DERIVED_SEER_COMBINED_N_2016, '^[cp]X') | DERIVED_AJCC_N_7TH_ED_2010  == 99 ~ 'X',
                              T ~ NA_character_ ),
           tnm.m = case_when ( 
                              str_detect(DERIVED_SEER_COMBINED_M_2016, '^[cp]0') | DERIVED_AJCC_M_7TH_ED_2010 %>% between (0, 10) ~ '0',
                              str_detect(DERIVED_SEER_COMBINED_M_2016, '^[cp]1') | DERIVED_AJCC_M_7TH_ED_2010   %>% between (100,199) ~ '1',
                              str_detect(DERIVED_SEER_COMBINED_M_2016, '^[cp]X') | DERIVED_AJCC_M_7TH_ED_2010   == 99 ~ 'X',
                              T ~ NA_character_ ),
           size.lt2015 = case_when ( 
                                    CS_TUMOR_SIZE_2004_2015  >= 0 & CS_TUMOR_SIZE_2004_2015 <= 989 ~  CS_TUMOR_SIZE_2004_2015 / 10, 
                                    CS_TUMOR_SIZE_2004_2015  == 0 ~  0, 
                                    CS_TUMOR_SIZE_2004_2015  == 991 ~  0.99, 
                                    CS_TUMOR_SIZE_2004_2015  == 992 ~  1.99, 
                                    CS_TUMOR_SIZE_2004_2015  == 993 ~  2.99, 
                                    CS_TUMOR_SIZE_2004_2015  == 994 ~  3.99, 
                                    CS_TUMOR_SIZE_2004_2015  == 995 ~  4.99, 
                                    T ~ NA_real_
                                    ),
           size.gt2015 = case_when ( 
                                    TUMOR_SIZE_SUMMARY_2016  > 0 & TUMOR_SIZE_SUMMARY_2016  <= 989 ~ TUMOR_SIZE_SUMMARY_2016 / 10,
                                    TUMOR_SIZE_SUMMARY_2016 == 990 ~ 0, 
                                    T  ~ NA_real_ 
           ),
           size = ifelse ( nna(size.lt2015), size.lt2015, size.gt2015),
           
           t_stage_8 = case_when (
             size>0 & size<1 & tnm.t==1  ~ 'T1a',
             size>=1 & size<2 & tnm.t==1 ~ 'T1b',
             size>=2 & size<3 & tnm.t==1 ~ 'T1c',
             (size>=3 & size<4) | (tnm.t=='2' & size<4) ~ 'T2a',
             size>=4 & size<5 ~ 'T2b',
             (size>=5 & size<7) | (tnm.t=='3' & size<7) ~ 'T3',
             (size >=7 & size<100) | tnm.t=='4'  ~ 'T4',
           T ~ NA_character_)
           ) %>% mutate ( 
           dx.date = ymd( ifelse ( nna(YEAR_OF_DIAGNOSIS) & nna(MONTH_OF_DIAGNOSIS) , sprintf('%d%02d15', YEAR_OF_DIAGNOSIS, MONTH_OF_DIAGNOSIS), NA_character_ ) )  ,
           death.date = ymd( ifelse ( ""!=(SEER_DATEOFDEATH_YEAR) & ""!=(SEER_DATEOFDEATH_MONTH) , sprintf('%s%s15', SEER_DATEOFDEATH_YEAR, SEER_DATEOFDEATH_MONTH), NA_character_ ) )  
           ) %>% arrange(dx.date) 

    #TODO: Is it even possible to get 30 day mortality if we don't have the date of death?
    table( A.gt2010$other.cause.mortality, A.gt2010$cause.specific.mortality, useNA="ifany")
A.gt2010 %>% count(histology.simple ,histology.cat)%>%arrange(-n) %>% print (n=Inf)
#A.gt2010 %>% group_by(histology.simple) %>% reframe((n()/415741)*100)
A.gt2010 %>% count() #415,741 Observations

#Double checking new variables that were added  
table(A.gt2010$t_stage_8, A.gt2010$tnm.t, useNA=c("ifany"))
A.gt2010 %>% filter(is.na(size)=="TRUE") %>% tally() #181,470 people have missing tumor size
A.gt2010 %>% filter(is.na(TUMOR_SIZE_SUMMARY_2016)=="TRUE", tnm.t==1, YEAR_OF_DIAGNOSIS>2015) %>% tally() #730 people have tnm.t==1 but missing tumor size information
A.gt2010 %>% filter(is.na(size)=="TRUE", tnm.t==1) %>% group_by(YEAR_OF_DIAGNOSIS) %>% tally() #730 people have tnm.t==1 but missing tumor size information
miss.t<-A.gt2010 %>% filter(is.na(size)=="TRUE", tnm.t==1, YEAR_OF_DIAGNOSIS>2015) #create dataset of patients with tnm.t=1 but who have missing size information
table(miss.t$YEAR_OF_DIAGNOSIS, miss.t$seer.surgery) 

#Checking the distribution of TNM T, N, and M staging variables across the years 
table(A.gt2010$tnm.t, A.gt2010$YEAR_OF_DIAGNOSIS, useNA=c("ifany"))
table(A.gt2010$t_stage_8, A.gt2010$YEAR_OF_DIAGNOSIS, useNA=c("ifany"))
table(A.gt2010$tnm.n, A.gt2010$YEAR_OF_DIAGNOSIS, useNA=c("ifany"))
table(A.gt2010$tnm.m, A.gt2010$YEAR_OF_DIAGNOSIS, useNA=c("ifany"))


label_list  <-  list(  
                     age = 'Age',  
                     sex = 'Sex',  
                     race = 'Race',  
                     marital.status = 'Marital status',  
                     other.cause.mortality = 'Other cause mortality',
                     cause.specific.mortality = 'Cause specific mortality',
                     primary.site = 'Primary site',
                     histology = 'Histology',
                     tnm.t = 'T stage',
                     tnm.n = 'N stage',
                     tnm.m = 'M stage',
                     t_stage_8 = 'T stage 8th edition',
                     histology.cat = 'Histology Categorical',
                     histology.simple = 'Histology Simple',
                     YEAR_OF_DIAGNOSIS = 'Year of Diagnosis',
                     BEHAVIOR_CODE_ICD_O_3 = 'Behavior'
                     
)


A   <-  A.gt2010 %>% select(PATIENT_ID, names(label_list) , dx.date, death.date, seer.surgery ) %>% distinct(PATIENT_ID, .keep_all =T) #this no longer changes anything. However, I kept it because everything downstream references 'A'



#Take a patient who has multiple observations for example
A.gt2010 %>% filter(PATIENT_ID=='lnK2020w0045894') %>% group_by(SEQUENCE_NUMBER, PRIMARY_SITE) %>% tally() %>% spread(PRIMARY_SITE, n) #patient has three primary cancers
A.gt2010 %>% filter(PATIENT_ID=='lnK2020w0196859') %>% group_by(SEQUENCE_NUMBER, PRIMARY_SITE) %>% tally() %>% spread(PRIMARY_SITE, n) #patient has two primary cancers

A.gt2010_lc <- A.gt2010 %>% filter(PRIMARY_SITE %in% valid.dxs) #goes from 495376 diagnoses to 440247 by removing all non-lung cancer diagnoses
A.gt2010_lc_ordered <-A.gt2010_lc[order(A.gt2010_lc$SEQUENCE_NUMBER, decreasing=FALSE),] #sort by sequence number 
A.gt2010_lc_firstonly <- A.gt2010_lc_ordered %>% distinct(PATIENT_ID, .keep_all = T) #goes from 440247 to 426195 by removing additional lung cancers diagnosed after the first lung cancer
A.gt2010_lc_firstonly %>% group_by(SEQUENCE_NUMBER) %>% tally()


################################
#  Inclusion Criteria 
################################
A2 <-A %>% filter(histology.cat!="Small Cell Carcinoma" & histology.cat!="Other/Unknown"& (t_stage_8=="T1a" | t_stage_8=="T1b" | t_stage_8=="T1c")) #Restricts to 50,613 Observations

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
medpar.carrier.tx  <- bind_rows ( medpar.tx, carrier.tx) %>% 
    left_join(A2)  %>%
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

#table( medpar.carrier.tx$tx.after.dx, useNA="ifany")
#table( medpar.carrier.tx$tx, medpar.carrier.tx$tx.after.dx, useNA="ifany")
medpar.carrier.tx.2 <- medpar.carrier.tx %>% filter( tx.after.dx & nna(tx.date) ) 

medpar.carrier.tx.2 %>% group_by(tx) %>% tally()
################################
#  Filter data and add in Thirty and Ninety Day Mortality Based on "tt"
################################
#TODO: Did Alex figure out the censoring?
A3 <- A2 %>% right_join(medpar.carrier.tx.2)  %>% mutate (
                                                          death = death.date,
                                                          tt = as.numeric( if_else ( nna(death.date), death.date, ymd('20191231')  ) - tx.date, units = 'days')
                                                          ) %>% 
                                filter( tt >0 ) %>% mutate( 
                                                           thirty.day.mortality = ifelse ( nna(death.date) & tt < 30, T, F ) ,
                                                           ninety.day.mortality = ifelse ( nna(death.date) & tt < 90, T, F ) 
                                    )

A3 %>% count(tx, thirty.day.mortality)
A3 %>% count(tx, ninety.day.mortality)

filename.out  <-  'data/A.final.all.3.RDS' #This is the final SEER file, combined with medpar and carrier data, with the inclusion criteria applied
 
 
#table( A3$tnm.n, useNA="ifany")
A3  <-  A3 %>% filter(  tnm.n== "0")
#A2  <-  A2 %>% filter( tnm.t == "1" &  primary.site == "Lung" & tnm.n== "0" & tnm.m =="0")
# filename.out  <-  'data/A.final.age.gte.80.RDS'
# A2  <-  A2 %>% filter( tnm.t == "1" &  primary.site == "Lung" & tnm.n== "0" & tnm.m =="0" & age >= 80 )

#filename.out  <-  'data/A.final.age.65.79.RDS'
 
 
#A2<-A.final 
#A2  <-  A2 %>% filter( tnm.t == "1" &  primary.site == "Lung" & tnm.n== "0" & tnm.m =="0" & age >= 65 & age < 80 )
A3 %>% count(tx)


################################
# Process the  Outpatient files
################################
#outpati  <-   read_dta('../SEER-Medicare-data/data/SEER_Medicare/outpat2014.base.dta')

fn.RDS  <- sprintf('%s/outpat.base.RDS', data.path)
if ( ! file.exists (fn.RDS) ) {
    outpats  <-  list()
    years  <-  as.character(2010:2019)
    for (yeari in 1:length(years)) {
        year  <-  years[yeari]
        print(year)
        outpati  <-   read_dta(sprintf('%s/outpat%s.base.dta', data.path, year), col_select=c('PATIENT_ID', 'CLM_FROM_DT', 'CLM_THRU_DT', PRNCPAL_DGNS_CD:PRCDR_DT25))
        # inner join with the SEER patients to reduce size
        outpats[[year]]  <-  outpati %>% 
            inner_join(lung.SEER.pids) %>%select( ! contains( "PRCDR_DT"))
    }
    outpat  <-  bind_rows(outpats,  .id='dataset.year')
    outpat  <-  outpat %>% mutate( CLM_FROM_DT = ymd(CLM_FROM_DT), CLM_THRU_DT = ymd(CLM_THRU_DT))
    saveRDS(object = outpat, file = fn.RDS) 
}else{
    outpat  <-  readRDS(fn.RDS)
}


# Combine outpat and medpar
medpar.dx <- medpar %>% select( PATIENT_ID, ADMSN_DT,DSCHRG_DT,  DGNS_1_CD:DGNS_25_CD) %>% set_names ( ~ str_replace_all(.,"DGNS_", "ICD_DGNS_CD") %>%  str_replace_all(.,"_CD$", "")) %>%
     mutate( CLM_FROM_DT = ymd(ADMSN_DT),
             CLM_THRU_DT = ymd(DSCHRG_DT))


outpat.medpar  <-  bind_rows ( outpat, medpar.dx) 
outpat.medpar$CLM_THRU_DT[  is.na( outpat.medpar$CLM_THRU_DT) ]  =  outpat.medpar$CLM_FROM_DT[  is.na( outpat.medpar$CLM_THRU_DT) ] 
outpat.medpar <- outpat.medpar %>% mutate( icd9or10 = ifelse( CLM_THRU_DT >= ymd('20151001'), 'icd10', 'icd9'  ))


################################
# Identify comorbidities 
################################

# First for ICD9
outpat.medpar.long  <- outpat.medpar  %>% right_join(A3) %>% filter( CLM_FROM_DT < tx.date ) %>% unite("ID_DATE", PATIENT_ID:CLM_FROM_DT, remove = F) %>%  select( ID_DATE, CLM_FROM_DT, ICD_DGNS_CD1:ICD_DGNS_E_CD12 )
outpat.medpar.long[ outpat.medpar.long == ""] = NA_character_



outpat.medpar.long.icd9  <-  outpat.medpar.long %>% filter ( CLM_FROM_DT < ymd('20151001') ) %>%  select(-CLM_FROM_DT) %>% pivot_longer( !ID_DATE , names_to= NULL, values_to = 'DX', values_drop_na = T) 
outpat.medpar.quan.deyo.icd9  <-  icd9_comorbid_quan_deyo(outpat.medpar.long.icd9, return_df = T)
smoking.hx.icd9  <-  outpat.medpar.long.icd9 %>% filter(DX %in% c( 'V1582', '3051') ) %>% mutate( Smoking = T)
o2.hx.icd9  <-  outpat.medpar.long.icd9 %>% filter(DX %in% c('V462') ) %>% mutate( o2 = T)
outpat.medpar.quan.deyo.icd9  <-  outpat.medpar.quan.deyo.icd9 %>% 
                                    left_join( smoking.hx.icd9 %>% select( - DX )) %>% replace_na( list(Smoking = F)) %>% 
                                    left_join( o2.hx.icd9 %>% select( - DX )) %>% replace_na( list(o2 = F))


outpat.medpar.long.icd10  <-  outpat.medpar.long %>% filter ( CLM_FROM_DT >= ymd('20151001') ) %>%  select(-CLM_FROM_DT) %>% pivot_longer( !ID_DATE , names_to= NULL, values_to = 'DX', values_drop_na = T) 
outpat.medpar.quan.deyo.icd10  <-  icd10_comorbid_quan_deyo(outpat.medpar.long.icd10, return_df = T)
smoking.hx.icd10  <-  outpat.medpar.long.icd10 %>% filter(DX %in% c( 'Z87891', expand_range('F17', 'F17299')) ) %>% mutate( Smoking = T)
o2.hx.icd10  <-  outpat.medpar.long.icd10 %>% filter(DX %in% c('Z9981') ) %>% mutate( o2 = T)
outpat.medpar.quan.deyo.icd10  <-  outpat.medpar.quan.deyo.icd10 %>% 
                                        left_join( smoking.hx.icd10 %>% select( - DX )) %>% replace_na( list(Smoking = F)) %>%
                                        left_join( o2.hx.icd10 %>% select( - DX )) %>% replace_na( list(o2 = F))

outpat.medpar.quan.deyo  <-  rbind(outpat.medpar.quan.deyo.icd9,outpat.medpar.quan.deyo.icd10) %>% as_tibble %>% separate (ID_DATE, c("PATIENT_ID", "CLM_FROM_DT"), sep = '_') %>% mutate( CLM_FROM_DT = ymd(CLM_FROM_DT)) %>% arrange( PATIENT_ID, CLM_FROM_DT)

outpat.medpar.quan.deyo.long  <-  outpat.medpar.quan.deyo  %>% replace(. == F, NA) %>% pivot_longer(-c(PATIENT_ID, CLM_FROM_DT), names_to = 'comorbidity', values_to = 'comorbidity.present', values_drop_na = T)

quan.deyo.final  <-  outpat.medpar.quan.deyo.long %>% 
                group_by(PATIENT_ID, comorbidity) %>% 
                mutate( time.from.last =  CLM_FROM_DT - first(CLM_FROM_DT)) %>% arrange(PATIENT_ID, comorbidity) %>% 
# Use this to require at >1 visits at certain time separation
                #summarise( meets.criteria = max( as.numeric(time.from.last, units='days') ) >= 30 ) %>% 
                summarise( meets.criteria = T) %>% 
                filter(meets.criteria) %>%
                pivot_wider( names_from =comorbidity, values_from = meets.criteria, values_fill = F )  


################################
# Merge with negative outcomes 
################################
# Create negative outcomes

# Create some common lists

expand_range('V01', 'V011')
# There are two kinds of negative outcomes: 1) Any outcome that occured after the treatment and 2) Any outcome that occured for the first time after the treatment.  Let's focus on the first for now. 

# Get common diagnoses
outpat.medpar.long.temp  <-  outpat.medpar %>%  filter (  CLM_THRU_DT >= ymd('20151001')) %>% right_join( A3)  %>% 
    select( PATIENT_ID, ICD_DGNS_CD1:ICD_DGNS_E_CD12) %>% 
    pivot_longer(!PATIENT_ID, 
                 names_to = 'diagnosis.field', values_to = 'diagnosis.code', 
                 values_drop_na = T) 

conditions  <-  outpat.medpar.long.temp %>% 
    filter ( diagnosis.code != "") %>% 
    mutate( explanation= explain_table(as.icd10( diagnosis.code),condense=F)$short_desc) %>%
    group_by(diagnosis.code, explanation) %>%
    summarise(n = n_distinct(PATIENT_ID)) %>%
    select(diagnosis.code, n, explanation) %>%
    arrange(-n)%>%
    top_n(1000)

conditions %>% write_tsv('tbls/potential.neg.outcomes.tsv')
conditions %>% print(n=Inf)






#fofo <- outpat.medpar %>% right_join(A2)   %>% mutate( noc.temp = find.rows.icdsmart( across(ICD_DGNS_CD1:ICD_DGNS_E_CD12), negative.outcomes[['fall']], icd9or10))
#testing  <- fofo %>% filter( ICD_DGNS_CD2 == 'E8889') 
#table( testing$noc.temp, testing$dataset.year, useNA="ifany")
#fofo %>% filter( ICD_DGNS_CD1 == 'E8889') %>% glimpse # other specified metabolic disorder
#explain_code(as.icd10('E8889'))
#explain_code(as.icd9('E8889'))
#fofo %>% filter( ICD_DGNS_CD1 == 'E889') %>% glimpse




#table( outpat.medpar$icd9or10, useNA="ifany")
#
#table( nna(outpat.medpar$CLM_THRU_DT), nna(outpat.medpar$CLM_FROM_DT), useNA="ifany")
##outpat.medpar <- outpat.medpar %>% mutate( icd9or10 = ifelse( CLM_THRU_DT >= ymd(20151001), 'icd10', 'icd9'  ))
##outpat.medpar  <-  outpat
#
#
#fifi  <-  medpar %>% filter( dataset.year == '2016') %>%  select(  DGNS_1_CD:DGNS_25_CD) %>% 
#    mutate( fifi.bool = find.rows( across(  DGNS_1_CD:DGNS_25_CD), c('E889')))
#fifi %>% filter(fifi.bool) %>%glimpse
#




#fofo  <- medpar %>% right_join(A2) %>%   filter( DGNS_1_CD == '57400')
#fofo %>% filter(ymd(ADMSN_DT) >ymd('20150101' )) %>% glimpse
#medpar %>% filter(PATIENT_ID == 'lnK2020w5173227') %>% glimpse


negative.outcomes  <-  list(
    'fall' = list( 
                  'icd9' = expand_range('E880','E888'),
                  'icd10' = expand_range('W00', 'W19' )),
    'other_injury' = list( 
                  'icd9' = expand_range('800', '999' ),
                  'icd10' = expand_range('S00','T79')),
    'acute_bronchitis' = list(
                              'icd9'=c('4660', '4661'),
                              'icd10' =  sprintf('J20%d',0:9)),
    'cholelithiasis' = list(
                            'icd9' = expand_range('5740', '57491'),
                            'icd10' =  expand_range('K80', 'K8081')),
    'oral' = list(
                            'icd9' = expand_range('520','5299'),
                            'icd10' =  expand_range('K00', 'K149')),
    'hpb' = list(
                            'icd9' = expand_range('570','577'),
                            'icd10' =  expand_range('K70', 'K87')),
    'gout' = list(
                            'icd9' = expand_range('2740','2749'),
                            'icd10' =  expand_range('M100', 'M109')),
    'arthropathy' = list(
                            'icd9' = c(expand_range('711', '715')),
                            'icd10' = c(  expand_range('M00', 'M19'))),
    'GU_sx' = list(
                            'icd9' = c(expand_range('590', '599'), expand_range('788','78899')),
                            'icd10' = c(  expand_range('R30', 'R39')), expand_range('N30', 'N39')),
    'diverticular_disease' = list(
                            'icd9' = expand_range('562','56213' ) ,
                            'icd10' =  expand_range('K57', 'K5793') ),
    'hernia' = list(
                            'icd9' = c( expand_range('550', '5539') ) ,
                            'icd10' =  expand_range('K40', 'K469') ),
    'hemorrhoids' = list(
                            'icd9' = c( expand_range('4550', '4559') ) ,
                            'icd10' =  expand_range('K640', 'K649') ),
    'optho' = list(
                            'icd9' = c( expand_range('360', '379') ) ,
                            'icd10' =  expand_range('H00', 'H59') )
)
expand_range('360', '379')

sink('tbls/negative.outcomes.txt'); 
for (namei in (names(negative.outcomes))) { 
    cat(sprintf('Variable Name: %s', namei))
    print('ICD9')
    icd9.codes  <-  negative.outcomes[[namei]][['icd9']]
    cat(sprintf('%s\n', paste(icd9.codes, explain_code(condense=F, as.icd9(icd9.codes)))))
    print('ICD10')
    icd10.codes  <-  negative.outcomes[[namei]][['icd10']]
    cat(sprintf('%s\n', paste(icd10.codes, explain_code(condense=F, as.icd10(icd10.codes)))))
    cat('------------------------------\n\n\n\n')
} 
sink()



#A3  <-  A2
#for (i in 1:length(negative.outcomes) ) {
#    noc.name  <- names(negative.outcomes)[i]
#    noc.codes  <- negative.outcomes[[noc.name]]
#    print(noc.name)
#    outpat.noc <- outpat.medpar %>% right_join(A2)  %>% 
#         mutate( 
#                noc.temp = find.rows.icdsmart( across(ICD_DGNS_CD1:ICD_DGNS_E_CD12), noc.codes, icd9or10),
#                noc.temp =  if_else(noc.temp & (CLM_FROM_DT > tx.date)& (CLM_FROM_DT < death.date), CLM_FROM_DT, ymd(NA_character_)) ) %>% 
#         group_by(PATIENT_ID) %>% 
#         summarise( !!noc.name := first(na.omit(noc.temp))) %>% 
#         filter(nna(!!rlang::sym(noc.name)))
#     A3 <- A3 %>% left_join(outpat.noc)
#}

A4  <-  A3
for (i in 1:length(negative.outcomes) ) {
    noc.name  <- names(negative.outcomes)[i]
    noc.codes  <- negative.outcomes[[noc.name]]
    print(noc.name)
    outpat.noc <- outpat.medpar %>% right_join(A3)  %>% 
         mutate( 
                noc.temp = find.rows.icdsmart( across(ICD_DGNS_CD1:ICD_DGNS_E_CD12), noc.codes, icd9or10),
                noc.temp.pre =  if_else(noc.temp & (CLM_FROM_DT < tx.date)& (CLM_FROM_DT < death.date), CLM_FROM_DT, ymd(NA_character_)), 
                noc.temp.post =  if_else(noc.temp & (CLM_FROM_DT > tx.date)& (CLM_FROM_DT < death.date), CLM_FROM_DT, ymd(NA_character_))  ,
                noc.temp.any =  if_else(noc.temp & (CLM_FROM_DT < death.date), CLM_FROM_DT, ymd(NA_character_))  
                ) %>% 
         group_by(PATIENT_ID) %>% 
         summarise( 
                      !!noc.name := first(na.omit(noc.temp.post)), 
                      !!sprintf('%s_pre', noc.name ) := first(na.omit(noc.temp.pre)),
                      !!sprintf('%s_any', noc.name ) := first(na.omit(noc.temp.any))
         )
     A4 <- A4 %>% left_join(outpat.noc)
}



# Combine with the comorbidities
comorbidities  <-  c('DM','DMcx', 'LiverMild', 'Pulmonary', 'PVD', 'CHF', 'MI', 'Renal', 'Stroke',  'PUD', 'Rheumatic', 'Dementia', 'LiverSevere', 'Paralysis', 'HIV', 'Smoking', 'o2')
A.final  <- A4 %>% left_join(quan.deyo.final)
A.final[,comorbidities] <- A.final[,comorbidities] %>% mutate_all( coalesce, F)
#sum(is.na(A.final))

tblcontrol <- tableby.control(numeric.stats = c('Nmiss', 'meansd'), numeric.simplify = T, cat.simplify =T, digits = 1,total = T,test = F)
f  <-  sprintf( 'tx ~ %s', paste( c(names(label_list),comorbidities), collapse = "+") )
labels(A.final)  <-  label_list
tt <- tableby(as.formula(f), data=A.final, control = tblcontrol)
summary(tt) %>% write2html('/PHShome/gcl20/Research_Local/SEER-Medicare/tbls/all_vars.htm')

#summary(tt, text=T) %>% as.data.frame %>% write_csv('output/table1.csv')
getwd()
write_rds( A.final,filename.out)
write_rds( label_list,'data/label.list.RDS')

