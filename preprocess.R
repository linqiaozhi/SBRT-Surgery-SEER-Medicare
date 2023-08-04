#install.packages('https://cran.r-project.org/src/contrib/Archive/SAScii/SAScii_1.0.tar.gz', repos=NULL, type="source")

library(arsenal)
library('SAScii')
library(lubridate)
library('tidyverse')
library(haven)
library(icd)#devtools::install_github("jackwasey/icd")
source('utilities.R')
source('codes.R')
source('DrugCodes.R')
source('file.paths.R')

################################
# New Plan: 
# I. Load the SEER file. Will filter the following files using patients with lung cancer in SEER, to speed things up a bit. Output: patient.seer, lung.SEER.pids
# II. Identify the treatment received by patients. Output: patient.tx
#       1. Inpatients: Load the Medpar file, and check for SBRT or Resection for inpatient
#       2. Outpatients: The outpatient file contains bills from institutions. The carrier file contains bills from providers. 
#             1. Carrier line file contains both procedure and the daignosis for which it was performed.        
#             2.  Outpatient Revenue file contains procedures, link it with Outpatient Line file to confirm diagnosis for each procedure is lung cancer.
#       3. Combine everything into a single dataframe patient.tx, with three two columns: PATIENT_ID, tx.date, tx
#       4. Create a center volume variable - 1) for SBRT and 2) for sublobar resection 
# III. Identify diagnoses (comorbidities and negative control outcomes). Output: patient.dx
#       1. Combine long versions of medpar, outpatient line, and carrier base. Filter each by patient.tx. Results in dataframes into all.long.icd9.dx, all.long.icd10.dx
#       2. Identify comorbidities of interest, patient.comorbidities
#       3. Identify negative outcome diagnoses of interest, patient.noc
# IV. Identify procedures (for now, just PET). Output: patient.proc
#       1. 
# V.  Part D files 
# VI. Identify death using MBSF, resulting in patient.mbsf.death
# VII. Combine patient.tx,  patient.dx, patient.proc, patient.mbsf.death, patient.seer
# VII. Massage variables
# IV. Filter





################################
#  SECTION I: Load SEER file
################################

lung.SEER <- read_dta(sprintf('%s/SEER.lung.cancer.dta', dta.path))
lung.SEER.valid.dx  <-  lung.SEER %>% filter(PRIMARY_SITE %in% valid.dxs) #exclude rows in the dataset corresponding to NON-lung cancer diagnoses (i.e., other cancers)
lung.SEER.ordered <-  lung.SEER.valid.dx[order(lung.SEER.valid.dx$SEQUENCE_NUMBER, decreasing=FALSE),] #sort by sequence number in ascending order
patient.seer <- lung.SEER.ordered %>% distinct(PATIENT_ID, .keep_all = TRUE) #keep the lung cancer diagnosis corresponding to the LOWEST sequence number (i.e., their first LC diagnosis)
lung.SEER.pids <- patient.seer %>% filter(YEAR_OF_DIAGNOSIS>=2010 & YEAR_OF_DIAGNOSIS<=2017) %>% select(PATIENT_ID) #restrict to patients diagnosed from 2010-2017; final SEER patient list; gives you a list of 415,741 patients (includes patients with a first primary LC diagnosis)



################################
# SECTION II: Identify treatment received by patients.
################################

################################
#  II.1 Inpatients 
################################
year = "2015"
fn.RDS  <- sprintf('%s/medpar.RDS', rds.path)
if ( ! file.exists (fn.RDS) ) {
    medpars  <-  list()
    years  <-  as.character(2009:2019)
    for (yeari in 1:length(years)) {
        year  <-  years[yeari]
        print(year)
        medpari  <-   read_dta(sprintf('%s/medpar%s.dta', data.path, year), col_select=c(PATIENT_ID, ADMSN_DT,  DSCHRG_DT, SRGCL_PRCDR_IND_SW, DGNS_1_CD:DGNS_25_CD, SRGCL_PRCDR_1_CD:SRGCL_PRCDR_25_CD, SRGCL_PRCDR_PRFRM_1_DT:SRGCL_PRCDR_PRFRM_25_DT, 'ORG_NPI_NUM'))
        # inner join with the SEER patients to reduce size
        medpars[[year]]  <-  medpari %>% inner_join(lung.SEER.pids) 
        medpars[[year]]$sbrt.date  <-  get.dates.of.procedure( medpars[[year]], sbrt.icds  )
        medpars[[year]]$sublobar.date  <-  get.dates.of.procedure( medpars[[year]], sublobar.icds  )
        medpars[[year]]$other.resection.date  <-  get.dates.of.procedure( medpars[[year]], other.resection.icds  )
        #medpars[[year]] <-   medpars[[year]]  %>% filter( !is.na(sbrt.date) | !is.na(sublobar.date)  )
    }
    medpar  <-  bind_rows(medpars ,  .id='dataset.year')
    rm(medpars); gc()
    saveRDS(object = medpar, file = fn.RDS)
}else{
    medpar  <-  readRDS(fn.RDS)
}

medpar %>% group_by(dataset.year) %>% tally() #check the number of observations per year

medpar <- medpar %>% mutate(
        actually.lung.cancer = find.rows( across(DGNS_1_CD:DGNS_25_CD), valid.dxs),  
        sbrt.date = ymd(sbrt.date), 
        sublobar.date = ymd(sublobar.date), 
        other.resection.date = ymd(other.resection.date) ) 
medpar$sbrt.date[ ! medpar$actually.lung.cancer ]  <- as.Date(NA_Date_)

################################
# II.2.1 Outpatients: Carrier File
################################
# We won't use the carrier base file for identifying procedures, but we will need it later.
years  <-  as.character(2010:2019)
carrierbases  <-  list()
for (yeari in 1:length(years)) {
  year  <-  years[yeari]
  fn  <-  sprintf('%s/nch%s.base.RDS',rds.path, year)
  if (!file.exists( fn )) {
    dta.fn  <-  sprintf('%s/nch%s.base.dta.gz', dta.path, year )
    print(sprintf('Reading in %s', dta.fn))
    carrierbasei  <-   read_dta(dta.fn, col_select=c('PATIENT_ID', 'CLM_FROM_DT', 'CLM_THRU_DT', 'PRNCPAL_DGNS_CD', ICD_DGNS_CD1:ICD_DGNS_CD12)) #
    carrierbasei.small  <- carrierbasei %>% inner_join(lung.SEER.pids)  
    saveRDS(object = carrierbasei.small,file = fn )
    # unlink(dta.fn)
  }else {
    print(sprintf('Reading in %s', fn))
    carrierbasei.small  <- readRDS(fn)
  }
  carrierbases[[year]]  <-  carrierbasei.small
}
carrierbase  <-  bind_rows(carrierbases,  .id='dataset.year')
rm(carrierbases); gc();

# Carrier line file contains procedures and their diagnosis codes
fn  <-  sprintf("%s/nch.lines.RDS", rds.path)
if (!file.exists( fn )) {
    years  <-  as.character(2010:2019)
    carriers  <-  list()
    for (yeari in 1:length(years)) {
        year  <-  years[yeari]
        dta.fn  <-  sprintf('../SEER-Medicare-data/data/SEER_Medicare/nch%s.line.dta.gz', year )
        print(sprintf('Reading in %s', dta.fn))
       carrieri  <-   read_dta(dta.fn, col_select=c('PATIENT_ID', 'CLM_THRU_DT', 'HCPCS_CD', 'LINE_ICD_DGNS_CD', 'ORG_NPI_NUM'))
        carriers[[year]]  <- carrieri %>% inner_join(lung.SEER.pids)  
    }
    carrier  <-  bind_rows(carriers,  .id='dataset.year')
    rm(carriers); gc()
    saveRDS(object = carrier,file = fn )
}else {
        carrier  <- readRDS(fn)
}

carrier <- carrier %>% mutate(
              valid.dx = LINE_ICD_DGNS_CD %in% valid.dxs ,
              sbrt =   HCPCS_CD %in% sbrt.cpts  & valid.dx,
              sbrt.date  = if_else ( sbrt, CLM_THRU_DT, NA_character_) %>% ymd 
)

table( nna(carrier$sbrt.date), useNA="ifany")



#####################################
# II.2.2 Outpatients: Outpatient files
#####################################
fn.RDS  <- sprintf("%s/outpat.revenue.RDS", rds.path )
if ( ! file.exists (fn.RDS) ) {
  revenue.outpats  <-  list()
  years  <-  as.character(2010:2019)
  for (yeari in 1:length(years)) {
    year  <-  years[yeari]
    print(year)
    revenue.outpati  <-   read_dta(sprintf('%s/outpat%s.revenue.dta.gz', data.path, year), col_select=c('PATIENT_ID','CLM_ID', 'CLM_THRU_DT', 'HCPCS_CD'))
    # inner join with the SEER patients to reduce size
    revenue.outpats[[year]]  <-  revenue.outpati %>% 
      inner_join(lung.SEER.pids) 
  }
  outpat.revenue  <-  bind_rows(revenue.outpats,  .id='dataset.year')
  rm(revenue.outpats); gc()
  outpat.revenue  <-  outpat.revenue %>% mutate(CLM_THRU_DT = ymd(CLM_THRU_DT))
  saveRDS(object = outpat.revenue, file = fn.RDS) 
}else{
  outpat.revenue  <-  readRDS(fn.RDS)
}

fn.RDS  <- sprintf("%s/outpat.base.RDS", rds.path)
if ( ! file.exists (fn.RDS) ) {
    outpats  <-  list()
    years  <-  as.character(2010:2019)
    for (yeari in 1:length(years)) {
        year  <-  years[yeari]
        print(year)
        outpati  <-   read_dta(sprintf('%s/outpat%s.base.dta', data.path, year), col_select=c('PATIENT_ID', 'CLM_ID', 'CLM_FROM_DT', 'CLM_THRU_DT', PRNCPAL_DGNS_CD:PRCDR_DT25, 'ORG_NPI_NUM'))
        # inner join with the SEER patients to reduce size
        outpats[[year]]  <-  outpati %>% 
            inner_join(lung.SEER.pids) %>%select( ! contains( "PRCDR_DT")) 
    }
    outpat  <-  bind_rows(outpats,  .id='dataset.year')
    outpat  <-  outpat %>% mutate( CLM_FROM_DT = ymd(CLM_FROM_DT), CLM_THRU_DT = ymd(CLM_THRU_DT))
    saveRDS(object = outpat, file = fn.RDS) 
    rm(outpats); gc()
}else{
    outpat  <-  readRDS(fn.RDS)
}
outpat %>% count (dataset.year)

#Check to see if the ICD procedure codes in outpatient base correspond to any SBRTs
#outpat %>% filter(ICD_PRCDR_CD1 %in% sbrt.icds) %>% count() #There are no patients who received SBRT 
#outpat %>% select(ICD_PRCDR_CD1) %>% tally() #All ICD_PRCRD_CD1 codes are empty

outpat.outpat.revenue  <-  outpat  %>% inner_join( outpat.revenue, by = c('PATIENT_ID', 'CLM_ID')) 
outpat.outpat.revenue <- outpat.outpat.revenue %>% mutate(
    valid.dx  =   PRNCPAL_DGNS_CD %in%  valid.dxs,
    sbrt  =   HCPCS_CD  %in%  sbrt.cpts  & valid.dx,
    sbrt.date  =  if_else ( sbrt, CLM_THRU_DT.x, as.Date(NA_Date_) )
)



################################
# II.3 Combine all three sources of treatment codes
################################
medpar.tx  <-   medpar %>% select ( PATIENT_ID, sbrt.date, sublobar.date, other.resection.date) %>% filter(!is.na(sbrt.date) | !is.na(sublobar.date) )
carrier.tx  <-  carrier %>% select(PATIENT_ID, sbrt.date) %>% filter(!is.na(sbrt.date))
outpat.tx  <-  outpat.outpat.revenue %>% select(PATIENT_ID, sbrt.date) %>% filter(!is.na(sbrt.date))
patient.tx  <- bind_rows ( medpar.tx, carrier.tx, outpat.tx)
patient.tx <- patient.tx %>% group_by(PATIENT_ID) %>% summarise ( 
               tx = factor( case_when ( 
                    any( nna( sbrt.date) ) & ! any( nna(sublobar.date))  ~ 'sbrt',
                    ! any( nna( sbrt.date) ) &  any( nna(sublobar.date)) ~ 'sublobar',
                    T ~ (NA_character_)
                    ), levels = c('sublobar', 'sbrt')),
               tx.date = case_when (
                                    tx == 'sbrt' ~ first(sbrt.date),
                                    tx == 'sublobar' ~ first(sublobar.date),
                                    T ~ ymd(NA_character_)
                                    ),
               other.resection.date = first(other.resection.date),
               ) 
table( patient.tx$tx, useNA="ifany")

################################
# Section III:  Identify diagnoses 
################################

# Create the long diagnosis data frame
outpat.dx  <-  outpat %>% 
    right_join(patient.tx %>% select( PATIENT_ID, tx.date) , by = 'PATIENT_ID')

carrierbase.dx  <- carrierbase  %>%  
    right_join(patient.tx%>% select( PATIENT_ID, tx.date) , by = 'PATIENT_ID') %>% 
    select( PATIENT_ID, tx.date, CLM_FROM_DT,CLM_THRU_DT,  ICD_DGNS_CD1:ICD_DGNS_CD12) %>%
    mutate(CLM_FROM_DT = ymd(CLM_FROM_DT), CLM_THRU_DT = ymd(CLM_THRU_DT))

medpar.dx <- medpar %>% 
    right_join(patient.tx%>% select( PATIENT_ID, tx.date) , by = 'PATIENT_ID') %>% 
    select( PATIENT_ID, tx.date, ADMSN_DT,DSCHRG_DT,  DGNS_1_CD:DGNS_25_CD) %>% 
    set_names ( ~ str_replace_all(.,"DGNS_", "ICD_DGNS_CD") %>%  str_replace_all(.,"_CD$", "")) %>%
     mutate( CLM_FROM_DT = ymd(ADMSN_DT),
             CLM_THRU_DT = ymd(DSCHRG_DT))


dx.wide  <-  bind_rows ( list(outpat=outpat.dx, medpar=medpar.dx, carrierbase=carrierbase.dx), .id ='source' )   %>% 
    mutate(across(where(is.character), ~ na_if(.,"")))
dx.wide$CLM_THRU_DT[  is.na( dx.wide$CLM_THRU_DT) ]  =  dx.wide$CLM_FROM_DT[  is.na( dx.wide$CLM_THRU_DT) ] 

dx.long  <- dx.wide %>% 
    unite("ID_DATE", c(PATIENT_ID,CLM_THRU_DT), remove = F) %>%  
    select( ID_DATE, tx.date, PATIENT_ID, CLM_THRU_DT,  ICD_DGNS_CD1:ICD_DGNS_E_CD12 ) %>%
    pivot_longer( !c(ID_DATE,tx.date,  PATIENT_ID, CLM_THRU_DT) , names_to= NULL, values_to = 'DX', values_drop_na = T)  %>% 
    distinct()

dx.long  <- dx.long %>% 
    mutate( icd9or10 = ifelse( CLM_THRU_DT >= ymd('20151001'), 'icd10', 'icd9'  ))

#dx.long.icd9.pre  <-  dx.long %>% filter( icd9or10 == 'icd9', CLM_THRU_DT < tx.date)
#dx.long.icd10.pre  <-  dx.long %>% filter( icd9or10 == 'icd10', CLM_THRU_DT < tx.date)

## Using the Quan comorbidity scores
#dx.long.quan.icd9  <-  icd9_comorbid_quan_deyo( dx.long.icd9.pre %>% select(ID_DATE, DX),
#                           return_df = T) 
#dx.long.quan.icd10  <-  icd10_comorbid_quan_deyo( dx.long.icd10.pre %>% select(ID_DATE, DX),
#                           return_df = T) 
#dx.long.quan  <-  rbind(dx.long.quan.icd9,dx.long.quan.icd10) %>% 
#    filter (rowSums(select(., MI:HIV)) > 0 ) %>%
#    as_tibble %>% 
#    separate (ID_DATE, c("PATIENT_ID", "CLM_THRU_DT"), sep = '_') %>% 
#    mutate( CLM_THRU_DT = ymd(CLM_THRU_DT)) %>% 
#    arrange( PATIENT_ID, CLM_THRU_DT)

#dx.long.quan.long  <-  dx.long.quan  %>% replace(. == F, NA) %>% 
#    pivot_longer(-c(PATIENT_ID, CLM_THRU_DT), 
#                 names_to = 'comorbidity', 
#                 values_to = 'comorbidity.present', 
#                 values_drop_na = T)

#dx.quan  <-  dx.long.quan.long %>% 
#                group_by(PATIENT_ID, comorbidity) %>% 
#                #mutate( time.from.last =  CLM_FROM_DT - first(CLM_FROM_DT)) %>% 
#                #arrange(PATIENT_ID, comorbidity) %>% 
## Use this to require at >1 visits at certain time separation
#                #summarise( meets.criteria = max( as.numeric(time.from.last, units='days') ) >= 30 ) %>% 
#                summarise( meets.criteria = T) %>% 
#                filter(meets.criteria) %>%
#                pivot_wider( names_from =comorbidity, values_from = meets.criteria, values_fill = F )  


## Using hardcoded codes

dx.hardcodeds  <- patient.tx %>% select(PATIENT_ID)
for (i in 1:length(dx.icd)) {
    print(sprintf('%d/%d hard coded dx', i, length(dx.icd)))
    dxoi  <- dx.icd[[i]]
    dx.name  <-  names(dx.icd)[i]
    dx.hardcoded  <- dx.long %>% 
             mutate( 
                    temp = if_else ( icd9or10 == 'icd9', DX %in% dxoi$icd9,DX %in% dxoi$icd10 ),
                    temp.pre =  if_else(temp & (CLM_THRU_DT < tx.date), CLM_THRU_DT, ymd(NA_character_)), 
                    temp.post =  if_else(temp & (CLM_THRU_DT > tx.date), CLM_THRU_DT, ymd(NA_character_))  ,
                    temp.any =  if_else(temp , CLM_THRU_DT, ymd(NA_character_))  
                    )  %>% 
             group_by(PATIENT_ID) %>% 
             summarise( 
                          !!dx.name := first(na.omit(temp.post)), 
                          !!sprintf('%s_pre', dx.name ) := first(na.omit(temp.pre)),
                          !!sprintf('%s_any', dx.name ) := first(na.omit(temp.any)),
                          !!sprintf('%s_pre_count', dx.name ) := length((na.omit(temp.pre))),
                          !!sprintf('%s_pre_date_count', dx.name ) := length(unique(na.omit(temp.pre))),
                          !!sprintf('%s_any_count', dx.name ) := length((na.omit(temp.any))),
                          !!sprintf('%s_any_date_count', dx.name ) := length(unique(na.omit(temp.any))),
                          !!sprintf('%s_post_count', dx.name ) := length((na.omit(temp.post))),
                          !!sprintf('%s_post_date_count', dx.name ) := length(unique(na.omit(temp.post)))
             )
             dx.hardcodeds  <- dx.hardcodeds %>% left_join(dx.hardcoded, by='PATIENT_ID')
}

# Get first dx of all time
first.dx  <- dx.long %>% arrange(PATIENT_ID, CLM_THRU_DT)
first.dx <- first.dx %>% group_by(PATIENT_ID) %>% summarise( first.dx.date = first(CLM_THRU_DT))

patient.dx   <-  dx.hardcodeds   %>%
    # left_join(dx.quan, by ='PATIENT_ID') %>%   
    # mutate(across(colnames(dx.quan), ~replace(., is.na(.), FALSE))) %>% 
    mutate(across(contains('count'), ~replace(., is.na(.), 0))) %>%
    left_join(first.dx)



################################
# SECTION IV Procedure codes (not for treatment)
################################
# Load DME
fn.RDS  <- sprintf("%s/dme.line.RDS", rds.path)
if ( ! file.exists (fn.RDS) ) {
    dme.lines  <-  list()
    years  <-  as.character(2010:2019)
    for (yeari in 1:length(years)) {
        year  <-  years[yeari]
        print(year)
        dme.linei  <-   read_dta(sprintf('%s/dme%s.line.dta', data.path, year) , 
                               col_select=c('PATIENT_ID', 'CLM_ID',  'CLM_THRU_DT', 'HCPCS_CD' ))
        dme.lines[[year]]  <-  dme.linei %>% inner_join(lung.SEER.pids)
    }
    dme.line  <-  bind_rows(dme.lines,  .id='dataset.year')
    dme.line  <-  dme.line %>% mutate(  CLM_THRU_DT = ymd(CLM_THRU_DT))
    saveRDS(object = dme.line, file = fn.RDS) 
}else{
    dme.line  <-  readRDS(fn.RDS)
}
dme.line %>% count (dataset.year)



# Combine DME with the other sources of HCPCS codes
carrier.proc   <- carrier %>% 
    select( PATIENT_ID, HCPCS_CD, CLM_THRU_DT) %>% 
    filter ( PATIENT_ID %in% patient.tx$PATIENT_ID) %>%
    mutate(CLM_THRU_DT = ymd(CLM_THRU_DT))
outpat.revenue.proc   <- outpat.revenue %>% 
    select( PATIENT_ID, HCPCS_CD, CLM_THRU_DT) %>% 
    filter ( PATIENT_ID %in% patient.tx$PATIENT_ID)
dme.line.proc  <-  dme.line %>% select( PATIENT_ID, HCPCS_CD, CLM_THRU_DT)
procs.long  <- rbind(carrier.proc, outpat.revenue.proc, dme.line.proc) %>% 
    left_join(patient.tx %>% 
    select(PATIENT_ID, tx.date))

# Look up each of the HCPCS codes from procs, defined in codes.R
proc.hardcodeds  <- patient.tx %>% select(PATIENT_ID)
procois  <- c(procs) 
for (i in 1:length(procois)) {
    print(sprintf('%d/%d procedures', i, length(procois)))
    proc.name  <-  names(procois)[i]
    procoi  <- procois[[i]]
    proc.hardcoded  <- procs.long %>% 
             mutate( 
                    temp =   HCPCS_CD %in% procoi,
                    temp.pre =  if_else(temp & (CLM_THRU_DT < tx.date), CLM_THRU_DT, ymd(NA_character_)), 
                    temp.post =  if_else(temp & (CLM_THRU_DT > tx.date), CLM_THRU_DT, ymd(NA_character_))  ,
                    temp.any =  if_else(temp , CLM_THRU_DT, ymd(NA_character_))  
                    )  %>% 
             group_by(PATIENT_ID) %>% 
             summarise( 
                          !!proc.name := first(na.omit(temp.post)), 
                          # !!sprintf('%s_pre', proc.name ) := first(na.omit(temp.pre)),
                          # !!sprintf('%s_any', proc.name ) := first(na.omit(temp.any)),
                          !!sprintf('%s_pre_count', proc.name ) := length((na.omit(temp.pre))),
                          # !!sprintf('%s_pre_date_count', proc.name ) := length(unique(na.omit(temp.pre))),
                          !!sprintf('%s_any_count', proc.name ) := length((na.omit(temp.any))),
                          # !!sprintf('%s_any_date_count', proc.name ) := length(unique(na.omit(temp.any))),
                          !!sprintf('%s_post_count', proc.name ) := length((na.omit(temp.post))),
                          # !!sprintf('%s_post_date_count', proc.name ) := length(unique(na.omit(temp.post)))
             )
             proc.hardcodeds  <- proc.hardcodeds %>% left_join(proc.hardcoded, by='PATIENT_ID')
}
pet.scans <- procs.long   %>%
        mutate(
            pet.scan    = HCPCS_CD %in% pet.scan.cpts,
            pet.scan.date = if_else(pet.scan, CLM_THRU_DT, as.Date(NA_Date_)),
            days.between.pet.and.treatment =  tx.date - pet.scan.date,
        ) %>% 
        filter (pet.scan) %>% 
        mutate(
               pet.scan.within.year = days.between.pet.and.treatment>=0 & days.between.pet.and.treatment<=360
               ) %>% 
        group_by(PATIENT_ID) %>% 
        summarise( valid.pet.scan = any(pet.scan.within.year))

patient.outpatient.procs  <-  proc.hardcodeds %>% left_join( pet.scans, by = 'PATIENT_ID') %>% replace_na(list( valid.pet.scan = F))  %>%
    mutate(across(contains('count'), ~replace(., is.na(.), 0)))



################################
# SECTION V Part D 
################################

fn.RDS  <- sprintf("%s/PartD.RDS", rds.path)
 if ( ! file.exists (fn.RDS) ) {
  PartD.files  <-  list()
  years  <-  as.character(2009:2019)
  for (yeari in 1:length(years)) {
    year  <-  years[yeari]
    print(year)
    PartD.filei  <-   read_dta(sprintf('%s/pdesaf%s.dta', data.path, year) , 
                               col_select=c('PATIENT_ID', 'SRVC_DT',  'PROD_SRVC_ID'))
    PartD.files[[year]]  <-  PartD.filei %>% inner_join(lung.SEER.pids)
  }
  PartD.file  <-  bind_rows(PartD.files,  .id='dataset.year')
  saveRDS(object = PartD.file, file = fn.RDS) 
} else{
  PartD.file  <-  readRDS(fn.RDS)
}


PartD.file.txdate  <-  PartD.file  %>% 
    right_join(patient.tx %>% select( PATIENT_ID, tx.date)) %>% 
    mutate(SRVC_DT = ymd(SRVC_DT))

#Create variables indicating the receipt and date of receipt of each drug included in the drug code list 
patient.drugs<-patient.tx %>% select(PATIENT_ID)
for (i in 1:length(drugs)) {
  print(sprintf('%d/%d Prescribed Drugs', i, length(drugs)))
  drugoi <- drugs[[i]]
  drug.name <- names(drugs)[i]
  drug <- PartD.file.txdate %>% 
    mutate(
      temp = PROD_SRVC_ID %in% drugoi,
      temp.pre = if_else(temp & (SRVC_DT < tx.date), SRVC_DT, ymd(NA_character_)),
      temp.post = if_else(temp & (SRVC_DT > tx.date), SRVC_DT, ymd(NA_character_)),
      temp.any = if_else(temp, SRVC_DT, ymd(NA_character_))
    ) %>% 
    group_by(PATIENT_ID) %>% 
    summarise(
      !!drug.name := first(na.omit(temp.post)),
      !!sprintf('%s_pre_count', drug.name) := length((na.omit(temp.pre))),
      !!sprintf('%s_any_count', drug.name ) := length((na.omit(temp.any))),
      !!sprintf('%s_post_count', drug.name) := length((na.omit(temp.post))),
    )
  patient.drugs <- patient.drugs %>% left_join(drug, by = 'PATIENT_ID')
}



################################
# SECTION VI MBSF death
################################

fn.RDS  <- sprintf("%s/MBSF.RDS", rds.path)
if ( ! file.exists (fn.RDS) ) {
  mbsfs  <-  list()
  years  <-  as.character(2010:2019)
  for (yeari in 1:length(years)) {
    year  <-  years[yeari]
    print(year)
        mbsfi  <-   read_dta(sprintf('%s/mbsf.abcd.summary.%s.dta', data.path, year), col_select=c('PATIENT_ID', 'BENE_DEATH_DT', 'VALID_DEATH_DT_SW', 'BENE_PTA_TRMNTN_CD', 'BENE_PTB_TRMNTN_CD', 'BENE_HI_CVRAGE_TOT_MONS', 'BENE_SMI_CVRAGE_TOT_MONS', 'BENE_ENROLLMT_REF_YR'))
    # inner join with the SEER patients to reduce size
    mbsfs[[year]]  <-  mbsfi %>% 
      inner_join(lung.SEER.pids) 
  }
  mbsf <-  bind_rows(mbsfs,  .id='dataset.year')
  saveRDS(object = mbsf, file = fn.RDS) 
}else{
  mbsf  <-  readRDS(fn.RDS)
}

patient.mbsf <- mbsf %>% mutate( 
                        death.date.mbsf = ifelse(mbsf$BENE_DEATH_DT!= "", mbsf$BENE_DEATH_DT, NA_Date_) %>% ymd )

# If MBSF is used for any other purpose, need to be more savvy with this step
patient.mbsf  <-  patient.mbsf %>% 
    filter (nna(death.date.mbsf)) %>% 
    group_by( PATIENT_ID) %>% 
    summarise( death.date.mbsf = first(death.date.mbsf))


################################
# Section VII Combine 
################################

A  <-  patient.tx %>% 
    left_join( patient.dx, by = 'PATIENT_ID') %>% 
    left_join( patient.outpatient.procs, by = 'PATIENT_ID',) %>% 
    left_join( patient.mbsf, by = 'PATIENT_ID') %>%
    left_join( patient.seer, by = 'PATIENT_ID') %>% 
    left_join( patient.drugs, by = 'PATIENT_ID')  



################################
# Section VIII Massage Variables 
################################

topography  <-  read_csv(file= './ICDO3topography.csv') %>% rename(site.topography = description) %>% mutate(PRIMARY_SITE = str_remove_all( icdo3_code, fixed(".")))

A  <-  A %>% #Use lung.SEER.first.lc.dx
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
                              str_detect(DERIVED_SEER_COMBINED_N_2016, '^[cp]1') | DERIVED_AJCC_N_7TH_ED_2010  %>% between (100,199) ~ '1', str_detect(DERIVED_SEER_COMBINED_N_2016, '^[cp]2') | DERIVED_AJCC_N_7TH_ED_2010  %>% between (200,299)~ '2',
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
                                  T ~ NA_character_),
           dx.date = ymd( ifelse ( nna(YEAR_OF_DIAGNOSIS) & nna(MONTH_OF_DIAGNOSIS) , sprintf('%d%02d15', YEAR_OF_DIAGNOSIS, MONTH_OF_DIAGNOSIS), NA_character_ ) )  ,
           #death.date = ymd( ifelse ( ""!=(SEER_DATEOFDEATH_YEAR) & ""!=(SEER_DATEOFDEATH_MONTH) , sprintf('%s%s15', SEER_DATEOFDEATH_YEAR, SEER_DATEOFDEATH_MONTH), NA_character_ ) ),
           cod.new = case_when(
             grepl("C", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) & !grepl("C34", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Other Cancer',
             grepl("C34", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Lung Cancer',
             grepl("J", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Respiratory Disease',
             grepl("I", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Circulatory Disease',
             grepl("A", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("B", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Infection_Parasite',
             grepl("E", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Endocrine Disorder',
             grepl("F", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Mental Disorder',
             grepl("G", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Nervous System Disease',
             grepl("H", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Eye and ear Diseases',
             grepl("K", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Digestive System Diseases',
             grepl("M", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Diseases of the musculoskeletal system and connective tissue',
             grepl("N", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Diseases of the genitourinary system',
             grepl("V", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("W", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X1", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X2", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X3", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X4", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X5", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X9", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("Y", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'External Causes of Morbidity Except Suicide',
             grepl("X6", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X7", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X8", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Suicide',
             CAUSE_OF_DEATH_ICD_10 == "" ~ 'Alive',
             CAUSE_OF_DEATH_ICD_10 == '7777' | CAUSE_OF_DEATH_ICD_10 == '7797' ~ 'Unknown',
             T ~ 'Other'
           ))

    A  <- A %>% mutate( 
                       death.date.seer = 
                           ymd( ifelse ( ""!=(SEER_DATEOFDEATH_YEAR) & ""!=(SEER_DATEOFDEATH_MONTH) , sprintf('%s%s15', SEER_DATEOFDEATH_YEAR, SEER_DATEOFDEATH_MONTH), NA_character_ ) )  ,
                       tt = as.numeric( if_else ( nna(death.date.mbsf), death.date.mbsf, ymd('20191231')  ) - tx.date, units = 'days'),
                       time.enrolled = as.numeric( if_else ( nna(death.date.mbsf), death.date.mbsf, ymd('20191231')  ) - first.dx.date, units = 'days'),
                       thirty.day.mortality = ifelse ( nna(death.date.mbsf) & tt < 30, T, F ) ,
                       ninety.day.mortality = ifelse ( nna(death.date.mbsf) & tt < 90, T, F ) ,
                       death = death.date.mbsf,
                       valid.death.indicator = case_when(
                                 is.na(death.date.seer) & is.na( death.date.mbsf)  ~ 'valid', # not death in either
                                 nna(death.date.seer) & nna( death.date.mbsf)  ~ 'valid', # dead in both
                                 nna(death.date.seer) & is.na( death.date.mbsf) ~ 'invalid', # Dead in SEER but not in MBSF is invalid
                                 nna(death.date.mbsf) & is.na( death.date.seer)  & year(death.date.mbsf) == 2019 ~ 'valid', # Dead in MBSF but not in SEER is valid if it occured in 2019
                                 nna(death.date.mbsf) & is.na( death.date.seer)  & year(death.date.mbsf) < 2019 ~ 'invalid', # Dead in MBSF but not in SEER is invalid if it occured <2019
                                                         ))

A  <- A %>% mutate(
                   Smoking = nna(smoking_pre),
                   Oxygen = nna(o2_pre) )

# cause.of.death<-case_when(
#   CAUSE_OF_DEATH_ICD_10 = 
# )

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
                     BEHAVIOR_CODE_ICD_O_3 = 'Behavior',
                     death      = 'Death',
                     thirty.day.mortality = '30-day mortality',
                     ninety.day.mortality = '90-day mortality',
                     cod.new = 'Cause of Death Category'
)

#A   <-  A.gt2010 %>% select(PATIENT_ID, names(label_list) , dx.date, death.date, seer.surgery ) %>% distinct(PATIENT_ID, .keep_all =T) #this no longer changes anything. However, I kept it because everything downstream references 'A'



################################
# Section IX  Exclusion
################################

A.final  <-  A %>% filter ( tx.date > dx.date & 
                           ! ( nna(other.resection.date) & other.resection.date < tx.date ) &
                            valid.death.indicator == 'valid')
A.final %>% count(tx)
A.final <- A.final %>% filter(  
                              age >=65,
                  histology.cat!="Small Cell Carcinoma" & histology.cat!="Other/Unknown" &
                  (t_stage_8=="T1a" | t_stage_8=="T1b" | t_stage_8=="T1c") & 
                  race != 'Unknown',
                  tnm.n==0 & tnm.m==0 
              )
A.final %>% count(tx)
A.final  <- A.final %>% filter ((valid.pet.scan & tx=='sbrt') | tx=='sublobar' ) 
A.final %>% count(tx)
A.final %>% count(year(tx.date), tx)

A.final %>% filter (tx =='sbrt') %>% count(year(tx.date))

tblcontrol <- tableby.control(numeric.stats = c('Nmiss', 'meansd'), numeric.simplify = T, cat.simplify =T, digits = 1,total = T,test = F)
f  <-  sprintf( 'tx ~ %s', paste( c(names(label_list), 
                sprintf('%s_any_count', c(names(dx.icd), names(procs), names(drugs) )), 
                sprintf('%s_pre_count', c(names(dx.icd), names(procs), names(drugs) ) )) , collapse = "+") )
labels(A.final)  <-  label_list
tt <- tableby(as.formula(f), data=A.final, control = tblcontrol)
summary(tt) %>% write2html('/PHShome/gcl20/Research_Local/SEER-Medicare/tbls/all_vars3.htm')

filename.out  <-  'data/A.final4.all.gte.65.RDS' 
A.final %>% select ( PATIENT_ID:death.date.mbsf, names(label_list), tt, time.enrolled, Insulin:Anticoags_post_count) %>%  write_rds( filename.out)
write_rds( label_list,'data/label.list.RDS')

table(A.final$DIAB_C_any_count >0, A.final$Insulin_any_date_count >0)
summary(A.final$Insulin_any_date_count)
A.final %>% group_by(tx) %>% summarise( m = mean((tt)))

A.final %>% group_by(year(tx.date)) %>% summarise ( mean( tx == 'sbrt')) 
count(year(tx.date), tx)

################################
# Testing 
################################
A.final %>% filter ( nna(METS))  %>% pull(PATIENT_ID) %>% head
A.final %>% filter ( nna(METS), is.na(cancer_nonlung_any)) %>% pull (PATIENT_ID) %>% head
A.final %>% filter (PATIENT_ID == 'lnK2020w0043216') %>% print(width=Inf)

dx.long %>% filter ( PATIENT_ID == 'lnK2020w0043216' , DX %in% dx.icd[['METS']]$icd9 | DX %in% dx.icd[['METS']]$icd10  ) 

table(nna(A.final$cancer_nonlung_any) , nna(A.final$METS))
#A.final %>% filter (tx == 'sbrt') %>% sample_n(1) %>% glimpse
#
#A.final %>% filter (PATIENT_ID == 'lnK2020y2293566') %>% glimpse
#outpat.outpat.revenue %>% filter (PATIENT_ID == 'lnK2020y2293566' & HCPCS_CD %in% sbrt.cpts) %>% glimpse
#outpat %>% filter (PATIENT_ID == 'lnK2020y2293566') %>% mutate(ex = explain_code(as.icd9(PRNCPAL_DGNS_CD),condense=F)) %>% select (CLM_THRU_DT, ex) %>% print(n=Inf)
#carrier %>% filter (PATIENT_ID == 'lnK2020y2293566') %>% mutate(ex = explain_code(as.icd9(LINE_ICD_DGNS_CD),condense=F)) %>% select (CLM_THRU_DT, ex) %>% print(n=Inf)



A.final %>% filter( other.cause.mortality == 'Death') %>% count(cod.new)
A.final %>% filter( cause.specific.mortality == 'Death') %>% count(cod.new)
