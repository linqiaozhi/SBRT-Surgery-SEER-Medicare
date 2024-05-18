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
#             1. Carrier base and line files contains both procedure and the daignosis for which it was performed.        
#             2.  Outpatient Revenue file contains procedures, link it with Outpatient Line file to confirm diagnosis for each procedure is lung cancer.
#       3. Combine everything into a single dataframe patient.tx, with three two columns: PATIENT_ID, tx.date, tx
# III. Identify diagnoses (comorbidities and negative control outcomes). Output: patient.dx
#       1. Combine long versions of medpar, outpatient line, and carrier base. Filter each by patient.tx. Results in dataframes into all.long.icd9.dx, all.long.icd10.dx
#       2. Identify comorbidities of interest, patient.comorbidities
#       3. Identify negative outcome diagnoses of interest, patient.noc
# IV. Identify procedure and DME codes (for now, just PET). Output: patient.proc
# V.  Identify medications using Part D files 
# VI. Identify death using MBSF, resulting in patient.mbsf.death
# VII. Combine patient.tx,  patient.dx, patient.proc, patient.mbsf.death, patient.seer
# VII. Massage variables
# IV. Filter





################################
#  SECTION I: Load SEER file
################################

lung.SEER <- read_dta(sprintf('%s/seerlung.%s.dta.gz', dta.path, suffix))
patient.seer  <-  lung.SEER %>% filter(PRIMARY_SITE %in% valid.dxs) %>% filter( YEAR_OF_DIAGNOSIS>=2010 & YEAR_OF_DIAGNOSIS<=2019 )
# Sequence number considers the patient's recorded lifetime: if only one
# primary cancer, then it is 0. But if there are multiple, then the first will
# be 1. So, to obtain patients with no prior malignancy history (ie. including
# those that might have developed a separate malignancy later in life), we
# include 0 and 1. 
patient.seer  <- patient.seer %>% filter(SEQUENCE_NUMBER == '00' | SEQUENCE_NUMBER == '01')
lung.SEER.pids <- patient.seer %>% select(PATIENT_ID) 

################################
# SECTION II: Identify treatment received by patients.
################################

################################
#  II.1 Inpatients 
################################
fn.RDS  <- sprintf('%s/medpar.RDS', rds.path)
if ( ! file.exists (fn.RDS) ) {
    medpars  <-  list()
    years  <-  as.character(2009:2020)
    for (yeari in 1:length(years)) {
        year  <-  years[yeari]
        print(year)
        medpari  <-   read_dta(sprintf('%s/medpar%s%s.dta.gz', data.path, year, suffix), col_select=c(PATIENT_ID, ADMSN_DT,  DSCHRG_DT, SRGCL_PRCDR_IND_SW, DGNS_1_CD:DGNS_25_CD, SRGCL_PRCDR_1_CD:SRGCL_PRCDR_25_CD, SRGCL_PRCDR_PRFRM_1_DT:SRGCL_PRCDR_PRFRM_25_DT, 'ORG_NPI_NUM'))
        # inner join with the SEER patients to reduce size
        medpars[[year]]  <-  medpari %>% inner_join(lung.SEER.pids) 
        medpars[[year]]$sbrt.date  <-  get.dates.of.procedure( medpars[[year]], sbrt.icds  )
        medpars[[year]]$sublobar.date  <-  get.dates.of.procedure( medpars[[year]], sublobar.icds  )
        medpars[[year]]$lobar.date  <-  get.dates.of.procedure( medpars[[year]], lobar.icds  )
        medpars[[year]]$other.resection.date  <-  get.dates.of.procedure( medpars[[year]], other.resection.icds  )
    }
    medpar  <-  bind_rows(medpars ,  .id='dataset.year')
    rm(medpars); gc()
    saveRDS(object = medpar, file = fn.RDS)
}else{
    medpar  <-  readRDS(fn.RDS)
}

medpar %>% group_by(dataset.year) %>% tally() #check the number of observations per year
sum(!is.na(medpar$sbrt.date))

medpar <- medpar %>% mutate(
        actually.lung.cancer = find.rows( across(DGNS_1_CD:DGNS_25_CD), valid.dxs),  
        sbrt.date = ymd(sbrt.date), 
        sublobar.date = ymd(sublobar.date), 
        lobar.date = ymd(lobar.date), 
        other.resection.date = ymd(other.resection.date) ) 
medpar$sbrt.date[ ! medpar$actually.lung.cancer ]  <- as.Date(NA_Date_)

################################
# II.2.1 Outpatients: Carrier File
################################
# We won't use the carrier base file for identifying procedures, but we will need it later.
years  <-  as.character(2010:2020)
carrierbases  <-  list()
for (yeari in 1:length(years)) {
  year  <-  years[yeari]
  fn  <-  sprintf('%s/nch%s.base.RDS',rds.path, year)
  if (!file.exists( fn )) {
    dta.fn  <-  sprintf('%s/nch%s%s.base.dta.gz', dta.path, year, suffix )
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
    years  <-  as.character(2010:2020)
    carriers  <-  list()
    for (yeari in 3:length(years)) {
        year  <-  years[yeari]
        dta.fn  <-  sprintf('%s/nch%s%s.line.dta.gz', dta.path,year , suffix)
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

#####################################
# II.2.2 Outpatients: Outpatient files
#####################################
fn.RDS  <- sprintf("%s/outpat.revenue.RDS", rds.path )
if ( ! file.exists (fn.RDS) ) {
  revenue.outpats  <-  list()
  years  <-  as.character(2010:2020)
  for (yeari in 1:length(years)) {
    year  <-  years[yeari]
    print(year)
    revenue.outpati  <-   read_dta(sprintf('%s/outpat%s%s.revenue.dta.gz', data.path, year, suffix), col_select=c('PATIENT_ID','CLM_ID', 'CLM_THRU_DT', 'HCPCS_CD'))
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
    years  <-  as.character(2010:2020)
    for (yeari in 1:length(years)) {
        year  <-  years[yeari]
        print(year)
        outpati  <-   read_dta(sprintf('%s/outpat%s%s.base.dta.gz', data.path, year, suffix), col_select=c('PATIENT_ID', 'CLM_ID', 'CLM_FROM_DT', 'CLM_THRU_DT', PRNCPAL_DGNS_CD:PRCDR_DT25, 'ORG_NPI_NUM'))
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
medpar.tx  <-   medpar %>% select ( PATIENT_ID, sbrt.date, sublobar.date, lobar.date, other.resection.date) %>% filter(!is.na(sbrt.date) | !is.na(sublobar.date) | ! is.na(lobar.date) )
carrier.tx  <-  carrier %>% select(PATIENT_ID, sbrt.date) %>% filter(!is.na(sbrt.date))
outpat.tx  <-  outpat.outpat.revenue %>% select(PATIENT_ID, sbrt.date) %>% filter(!is.na(sbrt.date))
patient.tx  <- bind_rows ( medpar.tx, carrier.tx, outpat.tx)
patient.tx <- patient.tx %>% arrange(PATIENT_ID, sbrt.date, sublobar.date, lobar.date) 
patient.tx <- patient.tx %>% group_by(PATIENT_ID) %>% 
    summarise ( 
               sbrt.date = first(na.omit(sbrt.date)),
               sublobar.date = first(na.omit(sublobar.date)),
               lobar.date = first(na.omit(lobar.date)),
               other.resection.date = first(na.omit(other.resection.date))
    )
patient.tx <- patient.tx %>% left_join(patient.seer %>% select (PATIENT_ID, RADIATION_RECODE, RX_SUMM_SURG_RAD_SEQ, RX_SUMM_SURG_PRIM_SITE_1998))


# Care must be taken when excluding patients who underwent other types of
# resection.  A patient is included in the sublobar group if he/she underwent
# SBRT as the primary treatment. If later he/she received some additional
# resection or even SBRT as rescue therapy, this patient should still be
# included. This is because we are simulating a target trial where patients are
# ``randomized'' at time of treatment and we are doing an intention to treat
# analysis.  Same goes for the SBRT group.
patient.tx <- patient.tx %>% mutate( 
                tx = factor( case_when ( 
                      nna( sbrt.date)  &  
                          ( is.na(sublobar.date) | sbrt.date < sublobar.date) & ( is.na(lobar.date) | sbrt.date < lobar.date) & ( is.na(other.resection.date) | sbrt.date < other.resection.date) ~ 'sbrt',
                      nna( sublobar.date)  &  
                          ( is.na(sbrt.date) | sublobar.date < sbrt.date) & ( is.na(lobar.date) | sublobar.date < lobar.date) & ( is.na(other.resection.date) | sublobar.date < other.resection.date) ~ 'sublobar',
                      nna( lobar.date)  &  
                          ( is.na(sbrt.date) | lobar.date < sbrt.date) & ( is.na(sublobar.date) | lobar.date < sublobar.date) & ( is.na(other.resection.date) | lobar.date < other.resection.date) ~ 'lobar',
                     T ~ (NA_character_)
                     ), levels = c('sublobar', 'sbrt', 'lobar')),
               tx.date = case_when (
                                    tx == 'sbrt' ~ sbrt.date,
                                    tx == 'sublobar' ~ sublobar.date,
                                    tx == 'lobar' ~ lobar.date,
                                    T ~ ymd(NA_character_)
                                    ),
               ) 


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

rm(dx.wide); gc()

year.month  <-  function(x) {
    sprintf('%d-%d', year(x), month(x) )
}

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
                          !!sprintf('%s_pre_count', dx.name ) := length((na.omit(temp.pre))),
                          !!sprintf('%s_pre_month_count', dx.name ) := length(unique(year.month(na.omit(temp.pre)))),
                          !!sprintf('%s_any_count', dx.name ) := length((na.omit(temp.any))),
                          !!sprintf('%s_any_month_count', dx.name ) := length(unique(year.month(na.omit(temp.any)))),
                          !!sprintf('%s_post_count', dx.name ) := length((na.omit(temp.post))),
                          !!sprintf('%s_post_month_count', dx.name ) := length(unique(year.month(na.omit(temp.post)))),
             )
             dx.hardcodeds  <- dx.hardcodeds %>% left_join(dx.hardcoded, by='PATIENT_ID')
}

patient.dx   <-  dx.hardcodeds   %>%
    mutate(across(contains('count'), ~replace(., is.na(.), 0)))

################################
# SECTION IV Procedure and DME codes (not for treatment)
################################
# Load DME
fn.RDS  <- sprintf("%s/dme.line.RDS", rds.path)
if ( ! file.exists (fn.RDS) ) {
    dme.lines  <-  list()
    years  <-  as.character(2010:2020)
    for (yeari in 1:length(years)) {
        year  <-  years[yeari]
        print(year)
        dme.linei  <-   read_dta(sprintf('%s/dme%s%s.line.dta.gz', dta.path, year, suffix) , 
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
                          !!sprintf('%s_pre_count', proc.name ) := length((na.omit(temp.pre))),
                          !!sprintf('%s_pre_month_count', proc.name ) := length(unique(year.month(na.omit(temp.pre)))),
                          !!sprintf('%s_any_count', proc.name ) := length((na.omit(temp.any))),
                          !!sprintf('%s_any_month_count', proc.name ) := length(unique(year.month(na.omit(temp.any)))),
                          !!sprintf('%s_post_count', proc.name ) := length((na.omit(temp.post))),
                          !!sprintf('%s_post_month_count', proc.name ) := length(unique(year.month(na.omit(temp.post)))),
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
  years  <-  as.character(2010:2020)
  for (yeari in 1:length(years)) {
    year  <-  years[yeari]
    print(year)
    PartD.filei  <-   read_dta(sprintf('%s/pdesaf%s%s.dta.gz', data.path, year,suffix) , 
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
      !!sprintf('%s_pre_month_count', drug.name ) := length(unique(year.month(na.omit(temp.pre)))),
      !!sprintf('%s_any_count', drug.name ) := length((na.omit(temp.any))),
      !!sprintf('%s_any_month_count', drug.name ) := length(unique(year.month(na.omit(temp.any)))),
      !!sprintf('%s_post_count', drug.name) := length((na.omit(temp.post))),
      !!sprintf('%s_post_month_count', drug.name ) := length(unique(year.month(na.omit(temp.post)))),
    )
  patient.drugs <- patient.drugs %>% left_join(drug, by = 'PATIENT_ID')
}



################################
# SECTION VI MBSF death
################################

fn.RDS  <- sprintf("%s/MBSF.RDS", rds.path)
if ( ! file.exists (fn.RDS) ) {
    mbsfs  <-  list()
    years  <-  as.character(2010:2020)
    for (yeari in 1:length(years)) {
        year  <-  years[yeari]
        print(year)
        mbsfi  <-   read_dta(sprintf('%s/mbsfabcd%s%s.dta.gz', data.path, year, suffix), col_select=c('PATIENT_ID', 'BENE_DEATH_DT', 'VALID_DEATH_DT_SW', 'BENE_PTA_TRMNTN_CD', 'BENE_PTB_TRMNTN_CD',sprintf('HMO_IND_%02d', 1:12), sprintf('MDCR_STATUS_CODE_%02d', 1:12), 'BENE_HI_CVRAGE_TOT_MONS', 'BENE_SMI_CVRAGE_TOT_MONS', 'BENE_ENROLLMT_REF_YR'))
        # inner join with the SEER patients to reduce size
        mbsfs[[year]]  <-  mbsfi %>% inner_join(lung.SEER.pids) 
    }
    mbsf <-  bind_rows(mbsfs,  .id='dataset.year')
    saveRDS(object = mbsf, file = fn.RDS) 
}else{
    mbsf  <-  readRDS(fn.RDS)
}

# Need to get the number of months pre-tx during which the patient was enrolled. 
# One record per patient per month enrolled
# https://resdac.org/cms-data/variables/medicare-status-code-january
mbsf.long  <- mbsf  %>%  select( PATIENT_ID,BENE_ENROLLMT_REF_YR,  MDCR_STATUS_CODE_01:MDCR_STATUS_CODE_12,HMO_IND_01:HMO_IND_12)  %>% rename_with(
                                     .fn = ~ gsub("MDCR_STATUS_CODE", "MDCRSTATUSCODE", .x),
                                     .cols = starts_with("MDCR_STATUS_CODE")
                                     ) %>%
                            rename_with(
                                        .fn = ~ gsub("HMO_IND", "HMOIND", .x),
                                        .cols = starts_with("HMO_IND")
                            )
      
mbsf.long.ffs  <- mbsf.long %>% pivot_longer(cols = -c(PATIENT_ID, BENE_ENROLLMT_REF_YR), names_to = c(".value", "month" ), names_sep="_" ) %>%
    mutate(year = BENE_ENROLLMT_REF_YR,
           date = ceiling_date(ymd(sprintf('%d-%s-01', year, month)), 'month') - days(1))
mbsf.pre.tx.months <- mbsf.long.ffs %>% right_join(patient.tx %>% select( PATIENT_ID, tx.date)) %>% 
                            filter (HMOIND %in% c(0,4) ) %>% 
                            group_by(PATIENT_ID) %>% summarise( pre.tx.months = sum( (date < tx.date))  )
patient.mbsf <- mbsf %>% mutate( 
                        death.date.mbsf.temp = ifelse(mbsf$BENE_DEATH_DT!= "", mbsf$BENE_DEATH_DT, NA_Date_) %>% ymd ) %>% 
                        group_by( PATIENT_ID) %>% 
                        summarise( 
                                  death.date.mbsf = first(na.omit(death.date.mbsf.temp))
                            ) %>%
                        right_join(mbsf.pre.tx.months)



################################
# Section VII Combine 
################################
summary( A$HISTOLOGIC_TYPE_ICD_O_3)
summary( A$DERIVEDSEERCOMBINED_T_2016_2017)
topography  <-  read_csv(file= './ICDO3topography.csv') %>% rename(site.topography = description) %>% mutate(PRIMARY_SITE = str_remove_all( icdo3_code, fixed(".")))
A  <-  patient.seer %>% 
    left_join( patient.dx, by = 'PATIENT_ID') %>% 
    left_join( patient.outpatient.procs, by = 'PATIENT_ID',) %>% 
    left_join( patient.mbsf, by = 'PATIENT_ID') %>%
    left_join( patient.tx %>% select( -RADIATION_RECODE, -RX_SUMM_SURG_PRIM_SITE_1998, - RX_SUMM_SURG_RAD_SEQ), by = 'PATIENT_ID') %>% 
    left_join( patient.drugs, by = 'PATIENT_ID')  
A  <-  A %>% 
    rename ( 
            age                    = AGERECODEWITHSINGLEAGES_AND_100,
            sex                    = SEX,
            race                    = RACE_RECODE_WHITE_BLACK_OTHER,
            marital.status          = MARITAL_STATUS_AT_DIAGNOSIS,
            cause.specific.mortality = SEERCAUSESPECIFICDEATHCLASSIFIC,
            other.cause.mortality = SEEROTHERCAUSEOFDEATHCLASSIFICA,
            ) %>%
    left_join( topography)   %>%
    mutate( 
           age = as.numeric(age),
           seer.surgery = RX_SUMM_SURG_PRIM_SITE_1998, # https://seer.cancer.gov/archive/manuals/2021/AppendixC/Surgery_Codes_Lung_2021.pdf
           sex = case_when ( 
                            sex == 1 ~ 'Male',
                            sex == 2 ~ 'Female', 
                            sex == 9 ~ 'Unknown' ),
           race = case_when ( 
                             race ==1 ~ 'White',
                             race ==2 ~ 'Black',
                             T ~ 'Other or unknown'),
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
           histology.code = sprintf( '%s/%s',  HISTOLOGIC_TYPE_ICD_O_3, BEHAVIOR_CODE_ICD_O_3 ),
           microscopically_confirmed = DIAGNOSTIC_CONFIRMATION %in% c(1,2,3,4),
           primary.site = case_when ( 
                                     str_detect(icdo3_code,'^C34*') ~ 'Lung',
                                     T ~ 'Other'
                                     ),
histology =  case_when(
                       grepl("^8250/3", histology.code) ~ 'Adenocarcinoma, lepidic',
                       histology.code %in% c('8550/3', '8551/3') ~ 'Adenocarcinoma, acinar',
                       grepl("^826", histology.code) ~ 'Adenocarcinoma, papillary',
                       histology.code %in% 
                                    c('8265/3', # micropapillary
                                    '8230/3', #solid cell
                                    '8253/3', #invasive mucinous adenocarcinoma
                                    '8254/3', #Mixed invasivemucinous adenocarcinoma
                                    '8254/3', #Mixed invasivemucinous adenocarcinoma
                                    '8333/3', #fetal none
                                    '8144/3', #enteric none
                                    '8256/3', #minimallly invasive none
                                    '8230/3', #Solid
                                    '8252/3', #Bronchiolo-alveolar
                                    '8323/3', #Mixed cell
                                    '8257/3') ~ 'Adenocarcinoma, NOS/other', #minimallly invasive none
                       grepl("^814", histology.code) ~ 'Adenocarcinoma, NOS/other',
                       grepl("^807", histology.code) ~ 'Squamous cell',
                       grepl("^8083/3", histology.code) ~ 'Squamous cell',
                       grepl("^8046", histology.code) ~ 'Non-small cell carcinoma',
                       grepl("^804", histology.code) ~ 'Small cell',
                       grepl("^8255/3", histology.code) ~ 'Adenocarcinoma, mixed',
                       grepl("^824", histology.code) ~ 'Carcinoid',
                       grepl("^856", histology.code) ~ 'Adenosquamous',
                       grepl("^848", histology.code) ~ 'Adenocarcinoma, mucinous',
                       histology.code %in% c('8012/3', '8013/3', '8014/3') ~ 'Large cell', # large cell
                       T ~ 'Other'),
           tnm.t = case_when ( 
                              str_detect(DERIVED_EOD_2018_T_2018, '^T1*') |str_detect(DERIVEDSEERCOMBINED_T_2016_2017, '^[cp]1') | DERIVED_AJCC_T_7TH_ED_2010_2015 %>% between (100,190) | DERIVED_AJCC_T_7TH_ED_2010_2015 %>% between(800, 810) ~ '1',
                              str_detect(DERIVED_EOD_2018_T_2018, '^T2*') |str_detect(DERIVEDSEERCOMBINED_T_2016_2017, '^[cp]2') | DERIVED_AJCC_T_7TH_ED_2010_2015  %>% between (200,290)~ '2',
                              str_detect(DERIVED_EOD_2018_T_2018, '^T3') |str_detect(DERIVEDSEERCOMBINED_T_2016_2017, '^[cp]3') | DERIVED_AJCC_T_7TH_ED_2010_2015  %>% between (300,390) ~ '3',
                              str_detect(DERIVED_EOD_2018_T_2018, '^T4') |str_detect(DERIVEDSEERCOMBINED_T_2016_2017, '^[cp]4') | DERIVED_AJCC_T_7TH_ED_2010_2015  %>% between (400,499) ~ '4',
                              str_detect(DERIVED_EOD_2018_T_2018, '^TX') |str_detect(DERIVEDSEERCOMBINED_T_2016_2017, '^[cp]X') | DERIVED_AJCC_T_7TH_ED_2010_2015  == 888 ~ 'X',
                              T ~ NA_character_ ),
           tnm.n = case_when ( 
                             str_detect(DERIVED_EOD_2018_N_2018, '^N0') | str_detect(DERIVEDSEERCOMBINED_N_2016_2017, '^[cp]0') | DERIVED_AJCC_N_7TH_ED_2010_2015  %>% between (0,40) ~ '0',
                             str_detect(DERIVED_EOD_2018_N_2018, '^N1') | str_detect(DERIVEDSEERCOMBINED_N_2016_2017, '^[cp]1') | DERIVED_AJCC_N_7TH_ED_2010_2015  %>% between (100,199) ~ '1', 
                             str_detect(DERIVED_EOD_2018_N_2018, '^N2') | str_detect(DERIVEDSEERCOMBINED_N_2016_2017, '^[cp]2') | DERIVED_AJCC_N_7TH_ED_2010_2015  %>% between (200,299)~ '2',
                             str_detect(DERIVED_EOD_2018_N_2018, '^N3') | str_detect(DERIVEDSEERCOMBINED_N_2016_2017, '^[cp]3') | DERIVED_AJCC_N_7TH_ED_2010_2015  %>% between (300,399) ~ '3',
                             str_detect(DERIVED_EOD_2018_N_2018, '^NX') | str_detect(DERIVEDSEERCOMBINED_N_2016_2017, '^[cp]X') | DERIVED_AJCC_N_7TH_ED_2010_2015  == 99 ~ 'X',
                              T ~ NA_character_ ),
           tnm.m = case_when ( 
                             str_detect(DERIVED_EOD_2018_M_2018, '^M0') |  str_detect(DERIVEDSEERCOMBINED_M_2016_2017, '^[cp]0') | DERIVED_AJCC_M_7TH_ED_2010_2015 %>% between (0, 10) ~ '0',
                             str_detect(DERIVED_EOD_2018_M_2018, '^M1*') |  str_detect(DERIVEDSEERCOMBINED_M_2016_2017, '^[cp]1') | DERIVED_AJCC_M_7TH_ED_2010_2015   %>% between (100,199) ~ '1',
                              T ~ NA_character_ ),
           CS_TUMOR_SIZE_2004_2015_num  = as.numeric(CS_TUMOR_SIZE_2004_2015),
           size.lt2015 = case_when ( 
                                    CS_TUMOR_SIZE_2004_2015_num  >= 0 & CS_TUMOR_SIZE_2004_2015_num <= 989 ~  CS_TUMOR_SIZE_2004_2015_num / 10, 
                                    CS_TUMOR_SIZE_2004_2015_num  == 0 ~  0, 
                                    CS_TUMOR_SIZE_2004_2015_num  == 991 ~  0.99, 
                                    CS_TUMOR_SIZE_2004_2015_num  == 992 ~  1.99, 
                                    CS_TUMOR_SIZE_2004_2015_num  == 993 ~  2.99, 
                                    CS_TUMOR_SIZE_2004_2015_num  == 994 ~  3.99, 
                                    CS_TUMOR_SIZE_2004_2015_num  == 995 ~  4.99, 
                                    T ~ NA_real_
                                    ),
           TUMOR_SIZE_SUMMARY_2016_num  = as.numeric(TUMOR_SIZE_SUMMARY_2016),
           size.gt2015 = case_when ( 
                                    TUMOR_SIZE_SUMMARY_2016_num  > 0 & TUMOR_SIZE_SUMMARY_2016_num  <= 989 ~ TUMOR_SIZE_SUMMARY_2016_num / 10,
                                    TUMOR_SIZE_SUMMARY_2016_num == 990 ~ 0, 
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
           dx.date = ymd( ifelse ( nna(YEAR_OF_DIAGNOSIS) & nna(MONTH_OF_DIAGNOSIS) , sprintf('%s%sd15', YEAR_OF_DIAGNOSIS, MONTH_OF_DIAGNOSIS), NA_character_ ) )  ,
           #cod.new = case_when(
           ## Note that for a SEQUENCE_NUMBER 0 cancer, the cause of death might
           #    # be listed as a different type of cancer, but
           #    # this is stil lconsidered lung cancer. If the
           #    # patient had a different primary, then the
           #    # sequence would be 1. If it's 0, and the cause
           #    # of death is, say, liver cancer, then the
           #    # thought is that this is actually just a met.
           #    # https://seer.cancer.gov/causespecific/
           #  grepl("C", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) & !grepl("C34", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Other Cancer', 
           #  grepl("C34", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Lung Cancer',
           #  grepl("J", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Respiratory Disease',
           #  grepl("I", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Circulatory Disease',
           #  grepl("A", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("B", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Infection_Parasite',
           #  grepl("E", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Endocrine Disorder',
           #  grepl("F", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Mental Disorder',
           #  grepl("G", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Nervous System Disease',
           #  grepl("H", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Eye and ear Diseases',
           #  grepl("K", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Digestive System Diseases',
           #  grepl("M", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Diseases of the musculoskeletal system and connective tissue',
           #  grepl("N", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Diseases of the genitourinary system',
           #  grepl("V", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("W", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X1", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X2", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X3", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X4", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X5", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X9", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("Y", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'External Causes of Morbidity Except Suicide',
           #  grepl("X6", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X7", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X8", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Suicide',
           #  CAUSE_OF_DEATH_ICD_10 == "" ~ 'Alive',
           #  CAUSE_OF_DEATH_ICD_10 == '7777' | CAUSE_OF_DEATH_ICD_10 == '7797' ~ 'Unknown',
#             T ~ 'Other'
#           )
)
    A  <- A %>% mutate( 
                       dx.to.tx  =  as.numeric( tx.date - dx.date, units = 'days'),
                       death.date.seer = 
                           ymd( ifelse ( ""!=(SEER_DATEOFDEATH_YEAR) & ""!=(SEER_DATEOFDEATH_MONTH) , sprintf('%s%s15', SEER_DATEOFDEATH_YEAR, SEER_DATEOFDEATH_MONTH), NA_character_ ) )  ,
                       tt = as.numeric( if_else ( nna(death.date.mbsf), death.date.mbsf, ymd('20191231')  ) - tx.date, units = 'days'),
                       thirty.day.mortality = ifelse ( nna(death.date.mbsf) & tt < 30, T, F ) ,
                       ninety.day.mortality = ifelse ( nna(death.date.mbsf) & tt < 90, T, F ) ,
                       death = death.date.mbsf, valid.death.indicator = case_when(
                                 is.na(death.date.seer) & is.na( death.date.mbsf)  ~ 'valid', # not death in either
                                 nna(death.date.seer) & nna( death.date.mbsf)  ~ 'valid', # dead in both
                                 nna(death.date.seer) & is.na( death.date.mbsf) ~ 'invalid', # Dead in SEER but not in MBSF is invalid
                                 nna(death.date.mbsf) & is.na( death.date.seer)  & year(death.date.mbsf) == 2019 ~ 'valid', # Dead in MBSF but not in SEER is valid if it occured in 2019
                                 nna(death.date.mbsf) & is.na( death.date.seer)  & year(death.date.mbsf) < 2019 ~ 'invalid', # Dead in MBSF but not in SEER is invalid if it occured <2019
                                                         ))
    # A %>% filter(tnm.t==1) %>% count(DERIVED_SEER_COMBINED_T_2016)

table( A$YEAR_OF_DIAGNOSIS,A$RADIATION_RECODE, useNA="ifany")
table( A$dx.date, useNA="ifany")
################################
# Section IX  Exclusion
################################
# A large number of SEER patients underwent surgery based on the seer.surgery variable which are not included, as they were note enrolled in Medicare at the time of treatment. Recall that SEER includes non Medicare patietns as well. 
incex  <-  function ( A.frame ) {
    print(sprintf('%d (SR: %d, SBRT: %d)', nrow(A.frame), sum(A.frame$tx == 'sublobar',na.rm = T), sum(A.frame$tx == 'sbrt', na.rm =T)))
}
print(nrow(A))
A.final  <- A %>% filter ( histology !="Small cell" & histology != 'Other' & histology != 'Carcinoid' & histology != 'Adenosquamous' & histology != 'Large cell' & histology != 'Non-small cell carcinoma' )
incex(A.final)
A.final  <-  A.final %>% filter ( (t_stage_8=="T1a" | t_stage_8=="T1b" | t_stage_8=="T1c") & tnm.n==0 & tnm.m==0 )
incex(A.final)
A.final  <-  A.final %>% filter (tx %in% c('sublobar', 'sbrt') )
A.final$tx  <- droplevels(A.final$tx)
incex(A.final)
A.final  <-  A.final %>% filter (
                                 (tx == 'sbrt' & ( (RADIATION_RECODE == "") | RADIATION_RECODE == "1") & (( RX_SUMM_SURG_PRIM_SITE_1998 != "") | RX_SUMM_SURG_PRIM_SITE_1998 == "0")  )  |
                                 ( tx == 'sublobar' & ( ( RX_SUMM_SURG_PRIM_SITE_1998 == "") |  RX_SUMM_SURG_PRIM_SITE_1998 %in% c("20", "21", "22", "23" ) ) & ((RADIATION_RECODE =="") | RADIATION_RECODE == "0") )
 )
A.final <- A.final %>% filter( RX_SUMM_SYSTEMIC_SURG_SEQ == "0")
incex(A.final)
table( A$age, useNA="ifany")
A.final  <-  A.final %>% filter (age >= 65 )
incex(A.final)
A.final  <-  A.final %>% filter (dx.to.tx <= 135 & dx.to.tx >= -16) # The diagnosis date is chosen to be the 15th of each month, as the date itself is not available
incex(A.final)
A.final  <-  A.final %>% filter (pre.tx.months >= 12 | is.na(tx))
incex(A.final)
A.final  <-  A.final %>% filter ((valid.pet.scan & tx=='sbrt') | tx=='sublobar')
incex(A.final)
A.final  <-  A.final %>% filter (valid.death.indicator == 'valid' )
incex(A.final)
A.final  <-  A.final %>% filter (microscopically_confirmed)
table( A.final$YEAR_OF_DIAGNOSIS, useNA="ifany")
incex(A.final)
#TODO: Cross each treatment variable iwth a year to make sure no errors

# incex(A.final)
# # summary(A.final$dx.to.tx)
# incex(A.final)
# table( A.final$histology, useNA="ifany")

# table( A.final$RX_SUMM_SYSTEMIC_SURG_SEQ,A.final$CHEMOTHERAPY_RECODE_YES_NO_UNK, useNA="ifany")
# table( A.final$histology, useNA="ifany")

# data.frame ( A.final$YEAR_THERAPY_STARTED, A.final$MONTH_OF_DIAGNOSIS, A.final$tx.date)

#TODO: GCL
# RX_Summ_Scope_Reg_LN_Sur_2003

A.final  %>%  write_rds( 'data/A.final.all.gte.65.RDS' )

table (A.final$cod.new)


table( A.final$arthropathy_pre_month_count, useNA="ifany")
table( A.final$pre.tx.months, useNA="ifany")
/ nrow(A.final)

A.final %>% count(histology.code,  histology, histology2) %>% arrange(-n) %>% write_csv('tbls/histology_codes.csv')

table( A.final$histology2, useNA="ifany")


################################
# Testing 
################################

table( A.final$tx, A.final$CHF_pre_date_count, useNA="ifany")
A.final %>% filter (CHF_pre_count > 200) %>% glimpse

outpat %>% filter ( PATIENT_ID == 'lnK2020x9593394')
fofo  <-  carrier %>% filter ( PATIENT_ID == 'lnK2020x9593394')  %>% mutate( desc = explain_code(LINE_ICD_DGNS_CD, condense =F ))
%>% print(n=30, width=Inf)


A.final %>% group_by(tx, DIAGNOSTIC_CONFIRMATION) %>% tally()
table(A.final$histology, A.final$microscopically_confirmed, useNA="ifany")
table( A.final$FIRSTMALIGNANTPRIMARY_INDICATOR, useNA="ifany")

A.final %>% filter (dx.to.tx > 135) %>% select(PATIENT_ID, tx.date, tx, dx.date)


A.final %>% filter ( nna(METS))  %>% pull(PATIENT_ID) %>% head
A.final %>% filter ( nna(METS), is.na(cancer_nonlung_any)) %>% pull (PATIENT_ID) %>% head
A.final %>% filter (PATIENT_ID == 'lnK2020w0043216') %>% print(width=Inf)

dx.long %>% filter ( PATIENT_ID == 'lnK2020w0043216' , DX %in% dx.icd[['METS']]$icd9 | DX %in% dx.icd[['METS']]$icd10  ) 

table( A.final$time.enrolled== 0, useNA="ifany")
table( year(A.final$tx.date), useNA="ifany")

table(nna(A.final$cancer_nonlung_any) , nna(A.final$METS))

################################
# Internal validity checks
################################

# We use the medicare codes to determine the cancer treatment. Make sure these agree with the SEER variables.
# Confirm that 
table( A.final$tx, A.final$RADIATION_RECODE, useNA="ifany")
table( A.final$tx, A.final$seer.surgery, useNA="ifany")
table( A.final$tx, A.final$RX_SUMM_SYSTEMIC_SURG_SEQ, useNA="ifany")
table( A.final$tx, A.final$RX_SUMM_SURG_RAD_SEQ, useNA="ifany")
table( A.final$tx, A.final$REASONNOCANCER_DIRECTED_SURGERY, useNA="ifany")

# Confirm that the month/year therapy started are somewhat consistent
# A.final %>% filter(PATIENT_ID == 

# with(A.final  ,     table( tx, RX_SUMM_SURG_PRIM_SITE_1998 , useNA="ifany"))
# with(A.final %>% filter (nna(first.month.coverage)) ,     table( tx, RX_SUMM_SURG_PRIM_SITE_1998 , useNA="ifany"))

##begin TODELETE
#    table( A.final$tx, A.final$PRIMARY_PAYER_AT_DX , useNA="ifany")
#    table( A.final$tx, A.final$RADIATION_RECODE, useNA="ifany")
#with(A.final %>% filter (nna(first.month.coverage)) ,     table( tx, RX_SUMM_SURG_PRIM_SITE_1998 , useNA="ifany"))
#    A.final %>% filter ( RX_SUMM_SURG_PRIM_SITE_1998 == 21  & is.na(tx))  %>% select(PATIENT_ID, dx.date) %>%  print(width=Inf, n=30)
#    A.final %>% filter ( RX_SUMM_SURG_PRIM_SITE_1998 == 33  & is.na(lobar.date))  %>% select(PATIENT_ID, dx.date) %>%  print(width=Inf, n=30)
#              0   12   13   15   19   20   21   22   23   24   25   30   33   45   46   55   56   65   70   80   90   99
#  sublobar   75    2    0    0    0    8 1653  553    4    0    0  197  900   21    3    2    5    1    1    0    4    0
#  sbrt     2350    3    1    6    0    0    9    2    0    0    0    2    3    0    0    0    0    0    0    0    2    5
#  lobar      58    0    0    0    0    5  152   44    4    0    0  702 3735   83    2    1    7    0    0    1    3    2
#  <NA>     3171   23    4   23    1    4  358   65    4    4    2  256 1380   23    2    2   20    1    0    0    4   20
#    A.final %>%filter ( PATIENT_ID == 'lnK2020w0028319') %>% t
#    medpar %>% filter ( PATIENT_ID == 'lnK2020w0028319') %>% t
#    outpat %>% filter ( PATIENT_ID == 'lnK2020w0266908')
#    mbsf.2016  <-   read_dta(sprintf('%s/mbsf.abcd.summary.%s.dta', data.path, year))
#    mbsf.2016 %>% print(width=Inf)
#    mbsf.2016 %>% filter( PATIENT_ID == 'lnK2020w0266908') %>% select(contains('HMO')) %>%  t
#    mbsf.2016 %>% filter( PATIENT_ID == 'lnK2020w0113668') %>% select(contains('HMO')) %>%  t
#    mbsf %>% filter ( PATIENT_ID == 'lnK2020w0028319') %>% print(width=Inf)



#    A.final %>%filter ( PATIENT_ID == 'lnK2020w0113668') %>% t
#MONTH_OF_DIAGNOSIS                          "7"                        
#YEAR_OF_DIAGNOSIS                           "2016"                     
#RX_SUMM_SURG_PRIM_SITE_1998                 "21"                       
#MONTH_OF_LAST_FOLLOW_UP_RECODE              "01"                       
#YEAR_OF_LAST_FOLLOW_UP_RECODE               "2018"                     
#MONTH_THERAPY_STARTED                       "09"                       
#YEAR_THERAPY_STARTED                        "2016"                     
#    medpar %>% filter ( PATIENT_ID == 'lnK2020w0113668') %>% t
#    carrier %>% filter ( PATIENT_ID == 'lnK2020w0113668') %>% print(n= Inf)
#    mbsf %>% filter ( PATIENT_ID == 'lnK2020w0113668') %>% print(n= Inf)

    #TODO: What is the survival variable of SEER?
    # A.final %>%filter ( PATIENT_ID == 'lnK2020w0067344') %>% t

    # outpat.outpat.revenue %>% filter ( PATIENT_ID == 'lnK2020w0113668') 
# end TODELETE




# A.fofo  <-  A.final2 %>% filter ( seer.surgery == 22 & is.na(tx) & PRIMARY_PAYER_AT_DX  == 60)
# A.fofo  <-  A.final2 %>% filter (dx.to.tx < -30)
# A.fofo  <-  A.final2  %>% filter (METS_pre_count >0 )
# table( A.final2$tx, A.final2$REASONNOCANCER_DIRECTED_SURGERY, useNA="ifany")
# table( A.final2$tx, A.final2$seer.surgery, useNA="ifany")
 # A.fofo %>% filter (PATIENT_ID == 'lnK2020w0153118') %>% glimpse
 # dx.long %>% filter (PATIENT_ID == 'lnK2020w0153118') %>% filter(tx.date > CLM_THRU_DT) %>%  count(DX) %>%arrange(-n) %>% print(n=100)
# A.fofo %>% filter (PATIENT_ID == 'lnK2020w0176506') %>% glimpse
# mbsf.long %>% filter (PATIENT_ID == 'lnK2020w0134942') %>% glimpse


################################
# Sensitivity analysis 1: Using SEER surgery treatment assignments. 
################################

A.final.sens1  <- A %>% filter ( histology !="Small cell" & histology != 'Other' & histology != 'Carcinoid' & histology != 'Adenosquamous' & histology != 'Large cell' & histology != 'Non-small cell carcinoma' )
 # A.final.sens1  <- A %>% filter ( histology !="Small cell" & histology != 'Other' & histology != 'Carcinoid' & histology != 'Large cell' & histology != 'Non-small cell carcinoma' )
A.final.sens1  <-  A.final.sens1 %>% filter ( (t_stage_8=="T1a" | t_stage_8=="T1b" | t_stage_8=="T1c") & 
                                             tnm.m ==0 &
                                             ( (tx == 'sbrt' & tnm.n == 0 ) |
                                               ( tx == 'sublobar' & tnm.n %in% c('0','1','2', '3') )
                                           ) )

A.final.




A.final.sens1  <-  A.final.sens1 %>% filter ( (t_stage_8=="T1a" | t_stage_8=="T1b" | t_stage_8=="T1c") & 
                                             tnm.m ==0 &
                                             ( (tx == 'sbrt' & tnm.n == 0 ) |
                                               ( tx == 'sublobar' & tnm.n %in% c('0','1','2', '3') )
                                           ) )
print(nrow(A.final.sens1))
 A.final.sens1  <-  A.final.sens1 %>% filter (tx %in% c('sublobar', 'sbrt') )
A.final.sens1$tx  <- droplevels(A.final.sens1$tx)
incex(A.final.sens1)
A.final.sens1 %>% filter ( tx  == 'sublobar') %>% count (RX_SUMM_SURG_RAD_SEQ)
A.final.sens1 %>% filter ( tx  == 'sublobar') %>% count (tnm.n,RADIATION_RECODE)
A.final.sens1  <-  A.final.sens1 %>% filter (
                                 (tx == 'sbrt' & ( is.na(RADIATION_RECODE) | RADIATION_RECODE == 1) & (is.na( RX_SUMM_SURG_PRIM_SITE_1998) | RX_SUMM_SURG_PRIM_SITE_1998 == 0)  )  |
                                 # ( tx == 'sublobar' & ( is.na( RX_SUMM_SURG_PRIM_SITE_1998) |  RX_SUMM_SURG_PRIM_SITE_1998 %in% c(20, 21, 22, 23 ) ) & (is.na(RADIATION_RECODE) | RADIATION_RECODE %in% c( 0, 3)) )
                                  ( tx == 'sublobar' & ( is.na( RX_SUMM_SURG_PRIM_SITE_1998) |  RX_SUMM_SURG_PRIM_SITE_1998 %in% c(20, 21, 22, 23 ) )  )
                             )
incex(A.final.sens1)
A.final.sens1 <- A.final.sens1 %>% filter(  ( tx == 'sbrt' & RX_SUMM_SYSTEMIC_SURG_SEQ == 0) | (tx == 'sublobar' & RX_SUMM_SYSTEMIC_SURG_SEQ %in% c(0, 3)))
incex(A.final.sens1)
A.final.sens1  <-  A.final.sens1 %>% filter (age >= 65 ) #TODO: <80
incex(A.final.sens1)
A.final.sens1  <-  A.final.sens1 %>% filter (dx.to.tx <= 135 & dx.to.tx >= -16) # The diagnosis date is chosen to be the 15th of each month, as the date itself is not available
incex(A.final.sens1)
# summary(A.final.sens1$dx.to.tx)
A.final.sens1  <-  A.final.sens1 %>% filter (pre.tx.months >= 12)
incex(A.final.sens1)
A.final.sens1  <-  A.final.sens1 %>% filter ((valid.pet.scan & tx=='sbrt') | tx=='sublobar')
incex(A.final.sens1)
A.final.sens1  <-  A.final.sens1 %>% filter (valid.death.indicator == 'valid' )
incex(A.final.sens1)
A.final.sens1  <-  A.final.sens1 %>% filter (microscopically_confirmed)
sum(A.final.sens1$tx == 'sublobar' & A.final.sens1$tnm.n %in% c('1','2','3')) / sum(A.final.sens1$tx == 'sublobar')
table( A.final.sens1$RX_SUMM_SYSTEMIC_SURG_SEQ, useNA="ifany")
table( A.final$tx, useNA="ifany")
A.final.sens1 %>% filter(tx == 'sublobar' ) %>% count(tnm.n,RX_SUMM_SCOPE_REG_LN_SUR_2003)
sum(A.final.sens1$tx == 'sublobar' & A.final.sens1$RX_SUMM_SCOPE_REG_LN_SUR_2003 ==0) / sum(A.final.sens1$tx == 'sublobar')

table( A.final$tx, useNA="ifany")

A.final.sens1  %>%  write_rds( 'data/A.final33.sens1.RDS' )


with(A.final.sens1%>%, table( tnm.n, nna(death)))






# print(nrow(A))
# A.final.sens1  <- A %>% filter ( histology !="Small cell")
# print(nrow(A.final.sens1))
#  A.final.sens1  <-  A.final.sens1 %>% filter (tx %in% c('sublobar', 'sbrt') )
# incex(A.final.sens1)
# A.final.sens1  <-  A.final.sens1 %>% filter (
#                                  (tx == 'sbrt' & ( is.na(RADIATION_RECODE) | RADIATION_RECODE == 1) & (is.na( RX_SUMM_SURG_PRIM_SITE_1998) | RX_SUMM_SURG_PRIM_SITE_1998 == 0)  )  |
#                                  ( tx == 'sublobar' & ( is.na( RX_SUMM_SURG_PRIM_SITE_1998) |  RX_SUMM_SURG_PRIM_SITE_1998 %in% c(20, 21, 22, 23 ) ) & (is.na(RADIATION_RECODE) | RADIATION_RECODE == 0) )
#                              )
# incex(A.final.sens1)
# A.final.sens1  <-  A.final.sens1 %>% filter ( (t_stage_8=="T1a" | t_stage_8=="T1b" | t_stage_8=="T1c") & 
#                                              tnm.m ==0 &
#                                              ( (tx == 'sbrt' & tnm.n == 0 ) |
#                                                ( tx == 'sublobar' & tnm.n %in% c('0','1','2', '3') )
#                                            ) 
# )
# incex(A.final.sens1)
# A.final.sens1 %>% filter ( tx  == 'sublobar') %>% count (RX_SUMM_SURG_RAD_SEQ)
# incex(A.final.sens1)
# A.final.sens1  <-  A.final.sens1 %>% filter (age >= 65 ) 
# incex(A.final.sens1)
# A.final.sens1  <-  A.final.sens1 %>% filter (dx.to.tx <= 135 & dx.to.tx >= -16) # The diagnosis date is chosen to be the 15th of each month, as the date itself is not available
# incex(A.final.sens1)
# # summary(A.final.sens1$dx.to.tx)
# A.final.sens1  <-  A.final.sens1 %>% filter (pre.tx.months >= 12)
# incex(A.final.sens1)
# A.final.sens1  <-  A.final.sens1 %>% filter ((valid.pet.scan & tx=='sbrt') | tx=='sublobar')
# incex(A.final.sens1)
# A.final.sens1  <-  A.final.sens1 %>% filter (valid.death.indicator == 'valid' )
# incex(A.final.sens1)
# A.final.sens1  <-  A.final.sens1 %>% filter (microscopically_confirmed)
# incex(A.final.sens1)
# sum(A.final.sens1$tx == 'sublobar' & A.final.sens1$tnm.n %in% c('1','2','3')) / sum(A.final.sens1$tx == 'sublobar')



# A.final.sens1  <- A %>% filter ( histology !="Small cell" & histology != 'Other' & histology != 'Carcinoid' & histology != 'Adenosquamous' )
# print(nrow(A.final.sens1))
# A.final.sens1  <-  A.final.sens1 %>% filter ( (t_stage_8=="T1a" | t_stage_8=="T1b" | t_stage_8=="T1c") & 
#                                              tnm.m ==0 &
#                                              ( (tx == 'sbrt' & tnm.n == 0 ) |
#                                                ( tx == 'sublobar' & tnm.n %in% c('0','1','2', '3') )
#                                            ) 
# )
# print(nrow(A.final.sens1))
#  A.final.sens1  <-  A.final.sens1 %>% filter (tx %in% c('sublobar', 'sbrt') )
# incex(A.final.sens1)
# A.final.sens1 %>% filter ( tx  == 'sublobar') %>% count (RX_SUMM_SURG_RAD_SEQ)
# A.final.sens1  <-  A.final.sens1 %>% filter (
#                                  (tx == 'sbrt' & ( is.na(RADIATION_RECODE) | RADIATION_RECODE == 1) & (is.na( RX_SUMM_SURG_PRIM_SITE_1998) | RX_SUMM_SURG_PRIM_SITE_1998 == 0)  )  |
#                                  ( tx == 'sublobar' & ( is.na( RX_SUMM_SURG_PRIM_SITE_1998) |  RX_SUMM_SURG_PRIM_SITE_1998 %in% c(20, 21, 22, 23 ) ) & (is.na(RADIATION_RECODE) | RADIATION_RECODE == 0) )
#                              )
# incex(A.final.sens1)
# A.final.sens1 <- A.final.sens1 %>% filter( RX_SUMM_SYSTEMIC_SURG_SEQ == 0)
# incex(A.final.sens1)
# A.final.sens1  <-  A.final.sens1 %>% filter (age >= 65 ) #TODO: <80
# incex(A.final.sens1)
# A.final.sens1  <-  A.final.sens1 %>% filter (dx.to.tx <= 135 & dx.to.tx >= -16) # The diagnosis date is chosen to be the 15th of each month, as the date itself is not available
# incex(A.final.sens1)
# # summary(A.final.sens1$dx.to.tx)
# A.final.sens1  <-  A.final.sens1 %>% filter (pre.tx.months >= 12)
# incex(A.final.sens1)
# A.final.sens1  <-  A.final.sens1 %>% filter ((valid.pet.scan & tx=='sbrt') | tx=='sublobar')
# incex(A.final.sens1)
# A.final.sens1  <-  A.final.sens1 %>% filter (valid.death.indicator == 'valid' )
# incex(A.final.sens1)
# A.final.sens1  <-  A.final.sens1 %>% filter (microscopically_confirmed)
# incex(A.final.sens1)
# A.final.sens1  %>%  write_rds( 'data/A.final20.sens1.RDS' )

# table( A.final$RX_SUMM_SYSTEMIC_SURG_SEQ,A.final$CHEMOTHERAPY_RECODE_YES_NO_UNK, useNA="ifany")
# table( A.final$histology, useNA="ifany")

# data.frame ( A.final$YEAR_THERAPY_STARTED, A.final$MONTH_OF_DIAGNOSIS, A.final$tx.date)












A.final.sens1  %>%  write_rds( 'data/A.final15.sens1.RDS' )


#A.final %>% filter (tx == 'sbrt') %>% sample_n(1) %>% glimpse
#
#A.final %>% filter (PATIENT_ID == 'lnK2020y2293566') %>% glimpse
#outpat.outpat.revenue %>% filter (PATIENT_ID == 'lnK2020y2293566' & HCPCS_CD %in% sbrt.cpts) %>% glimpse
#outpat %>% filter (PATIENT_ID == 'lnK2020y2293566') %>% mutate(ex = explain_code(as.icd9(PRNCPAL_DGNS_CD),condense=F)) %>% select (CLM_THRU_DT, ex) %>% print(n=Inf)
#carrier %>% filter (PATIENT_ID == 'lnK2020y2293566') %>% mutate(ex = explain_code(as.icd9(LINE_ICD_DGNS_CD),condense=F)) %>% select (CLM_THRU_DT, ex) %>% print(n=Inf)



A.final %>% filter( other.cause.mortality == 'Death') %>% count(cod.new)
A.final %>% filter( cause.specific.mortality == 'Death') %>% count(cod.new)


# label_list  <-  list(  
#                      age = 'Age',  
#                      sex = 'Sex',  
#                      race = 'Race',  
#                      marital.status = 'Marital status',  
#                      other.cause.mortality = 'Other cause mortality',
#                      cause.specific.mortality = 'Cause specific mortality',
#                      primary.site = 'Primary site',
#                      # histology = 'Histology',
#                      tnm.t = 'T stage',
#                      tnm.n = 'N stage',
#                      tnm.m = 'M stage',
#                      t_stage_8 = 'T stage 8th edition',
#                      # histology.cat = 'Histology Categorical',
#                      histology = 'Histology',
#                      YEAR_OF_DIAGNOSIS = 'Year of Diagnosis',
#                      BEHAVIOR_CODE_ICD_O_3 = 'Behavior',
#                      death      = 'Death',
#                      thirty.day.mortality = '30-day mortality',
#                      ninety.day.mortality = '90-day mortality',
#                      cod.new = 'Cause of Death Category'
# )
# write_rds( label_list,'data/label.list.RDS')
# tblcontrol <- tableby.control(numeric.stats = c('Nmiss', 'meansd'), numeric.simplify = T, cat.simplify =T, digits = 1,total = T,test = F)
# # f  <-  sprintf( 'tx ~ %s', paste( c(names(label_list), 
# #                 sprintf('%s_any_count', c(names(dx.icd), names(procs), names(drugs) )), 
# #                 sprintf('%s_pre_count', c(names(dx.icd), names(procs), names(drugs) ) )) , collapse = "+") )
# # labels(A.final)  <-  label_list
# # tt <- tableby(as.formula(f), data=A.final, control = tblcontrol)
# # summary(tt) %>% write2html('/PHShome/gcl20/Research_Local/SEER-Medicare/tbls/all_vars11.htm')


#################################
## TO DELETE 
#################################

#A.old  <-  patient.tx %>% 
#    left_join( patient.dx, by = 'PATIENT_ID') %>% 
#    left_join( patient.outpatient.procs, by = 'PATIENT_ID',) %>% 
#    left_join( patient.mbsf, by = 'PATIENT_ID') %>%
#    left_join( patient.seer, by = 'PATIENT_ID') %>% 
#    left_join( patient.drugs, by = 'PATIENT_ID')  


#topography  <-  read_csv(file= './ICDO3topography.csv') %>% rename(site.topography = description) %>% mutate(PRIMARY_SITE = str_remove_all( icdo3_code, fixed(".")))
#A.old  <-  A.old %>% #Use lung.SEER.first.lc.dx
#    filter(YEAR_OF_DIAGNOSIS >=2010) %>%  
#    rename ( 
#            age                    = AGERECODEWITHSINGLEAGES_AND_100,
#            sex                    = SEX,
#            race                    = RACE_RECODE_WHITE_BLACK_OTHER,
#            marital.status          = MARITAL_STATUS_AT_DIAGNOSIS,
#            cause.specific.mortality = SEERCAUSESPECIFICDEATHCLASSIFIC,
#            other.cause.mortality = SEEROTHERCAUSEOFDEATHCLASSIFICA,
#            seer.surgery = RX_SUMM_SURG_PRIM_SITE_1998, # https://seer.cancer.gov/archive/manuals/2021/AppendixC/Surgery_Codes_Lung_2021.pdf
#            ) %>%
#    left_join( topography)   %>%
#    mutate( 
#           sex = case_when ( 
#                            sex == 1 ~ 'Male',
#                            sex == 2 ~ 'Female', 
#                            sex == 9 ~ 'Unknown' ),
#           race = case_when ( 
#                             race ==1 ~ 'White',
#                             race ==2 ~ 'Black',
#                             T ~ 'Other or unknown'),
#           marital.status = case_when ( 
#                                       marital.status ==1 ~ 'Never married',
#                                       marital.status ==2 ~ 'Married',
#                                       marital.status %>% between( 3,6) ~ 'Other',
#                                       T ~ 'Unknown'),
#           cause.specific.mortality = case_when ( 
#                                                 cause.specific.mortality == 0 ~ 'Alive or other death',
#                                                 cause.specific.mortality ==1 ~ 'Death',
#                                                 T ~ 'Unknown' ),
#           other.cause.mortality = case_when ( 
#                                              other.cause.mortality == 0 ~ 'Alive or cancer death',
#                                              other.cause.mortality ==1 ~ 'Death',
#                                              T ~ 'Unknown' ),
#           histology.code = sprintf( '%d/%d',  HISTOLOGIC_TYPE_ICD_O_3, BEHAVIOR_CODE_ICD_O_3 ),
#           histology = case_when ( 
#                                  histology.code  == '8140/3' ~ 'Adenocarcinoma, NOS',
#                                  # INvasive non-mucinous adenocarcinomas
#                                  histology.code  == '8250/3' ~ 'Adenocarcinoma, lepidic',
#                                  histology.code  == '8550/3' ~ 'Adenocarcinoma, acinar',
#                                  histology.code  == '8250/3' ~ 'Adenocarcinoma, papillary',
#                                  #histology.code  == '8265/3' ~ 'Adenocarcinoma, micropapillary', # there were not any
#                                  #histology.code  == '8230/3' ~ 'Adenocarcinoma, solid',
#                                  #Invasive mucinous adenocarcinoma
#                                  histology.code  == '8253/3' |histology.code  == '8480/3'  ~ 'Adenocarcinoma, mucinous',
#                                  histology.code  == '8255/3' ~ 'Adenocarcinoma, mixed',
#                                  histology.code  %in%  c('8070/3', '8071/3', '8073/3', '8075/3', '8076/3', '8078/3') ~ 'Squamous cell carcinoma',
#                                  histology.code  == '8041/3' ~ 'Small cell carcinoma',
#                                  histology.code  == '8000/3' ~ 'Malignancy, unspecified',
#                                  histology.code  == '8046/3' ~ 'Non-small cell carcinoma',
#                                  histology.code  == '8010/3' ~ 'Carcinoma, unspecified',
#                                  T ~ 'Other' ),
#           primary.site = case_when ( 
#                                     str_detect(icdo3_code,'^C34*') ~ 'Lung',
#                                     T ~ 'Other'
#                                     ),
#           ##Histology Categorical (based on information from -- INSERT CITATION)
#           #histology.cat = case_when (
#           #  histology.code %in% c('8140/3', '8550/3',  '8551/3',  '8260/3', '8230/3', '8333/3', '8144/3', '8480/3', '8253/3', '8254/3') ~ 'Adenocarcinoma',
#           #  histology.code %in%  c('8070/3', '8071/3', '8072/3', '8073/3', '8074/3', '8075/3', '8076/3', '8078/3', '8083/3') ~ 'Squamous Cell Carcinoma',
#           #  histology.code %in% c('8012/3', '8013/3', '8014/3') ~ 'Large Cell Carcinoma',
#           #  histology.code == '8560/3' ~ 'Adenosquamous Cell Carcincoma',
#           #  histology.code %in% c('8240/3',	'8241/3',	'8242/3',	'8243/3',	'8244/3',	'8245/3', '8246/3',	'8249/3') ~ 'Carcinoid Tumor',
#           #  histology.code %in% c('8250/3', '8251/3', '8252/3', '8255/3') ~ 'Tumors Formerly Classifed as Bronchioloalveolar Carcinoma',
#           #  histology.code == '8046/3' ~ 'Non-small Cell Carcinoma, NOS',
#           #  histology.code == '8041/3' ~ 'Small Cell Carcinoma',
#           #  T ~ 'Other/Unknown'),
#           # histology.simple=case_when(
#           #   histology.cat=='Adenocarcinoma' ~ 'Adenocarcinoma',
#           #   histology.cat=='Squamous Cell Carcinoma' ~ 'Squamous Cell Carcinoma',
#           #   histology.cat=='Tumors Formerly Classifed as Bronchioloalveolar Carcinoma' ~ 'Tumors Formerly Classifed as Bronchioloalveolar Carcinoma', 
#           #   histology.cat=='Non-small Cell Carcinoma, NOS' ~ 'Non-small Cell Carcinoma, NOS',
#           #   histology.cat=='Adenosquamous Cell Carcinoma' | histology.cat=='Large Cell Carcinoma' | histology.cat=='Carcinoid' ~ 'Other/Unknown',
#           #   T ~ 'Other/Unknown'
#           # ),
#           histology =  case_when(
#                                   grepl("^814", histology.code) ~ 'Adenocarcinoma, NOS',
#                                   grepl("^807", histology.code) ~ 'Squamous cell',
#                                   grepl("^8046", histology.code) ~ 'Non-small cell carcinoma',
#                                   grepl("^8041", histology.code) ~ 'Small cell',
#                                   grepl("^8255/3", histology.code) ~ 'Adenocarcinoma, mixed',
#                                   grepl("^8550/3", histology.code) ~ 'Adenocarcinoma, acinar',
#                                   grepl("^824", histology.code) ~ 'Carcinoid',
#                                   grepl("^856", histology.code) ~ 'Adenosquamous',
#                                   grepl("^848", histology.code) ~ 'Adenocarcinoma, mucinous',
#                                   grepl("^8250/3", histology.code) ~ 'Adenocarcinoma, lepidic',
#                                   grepl("^8013/34", histology.code) ~ 'Large cell neuroendocrine',
#                                   grepl("^8012/34", histology.code) ~ 'Large cell carcinoma, NOS',
#                                   T ~ 'Other'),
#           tnm.t = case_when ( 
#                              str_detect(DERIVED_SEER_COMBINED_T_2016, '^[cp]1') | DERIVED_AJCC_T_7TH_ED_2010 %>% between (100,190) | DERIVED_AJCC_T_7TH_ED_2010 %>% between(800, 810) ~ '1',
#                              str_detect(DERIVED_SEER_COMBINED_T_2016, '^[cp]2') | DERIVED_AJCC_T_7TH_ED_2010  %>% between (200,290)~ '2',
#                              str_detect(DERIVED_SEER_COMBINED_T_2016, '^[cp]3')| DERIVED_AJCC_T_7TH_ED_2010  %>% between (300,390) ~ '3',
#                              str_detect(DERIVED_SEER_COMBINED_T_2016, '^[cp]4')| DERIVED_AJCC_T_7TH_ED_2010  %>% between (400,499) ~ '4',
#                              str_detect(DERIVED_SEER_COMBINED_T_2016, '^[cp]X') | DERIVED_AJCC_T_7TH_ED_2010  == 888 ~ 'X',
#                              T ~ NA_character_ ),
#           tnm.n = case_when ( 
#                              str_detect(DERIVED_SEER_COMBINED_N_2016, '^[cp]0') | DERIVED_AJCC_N_7TH_ED_2010  %>% between (0,40) ~ '0',
#                              str_detect(DERIVED_SEER_COMBINED_N_2016, '^[cp]1') | DERIVED_AJCC_N_7TH_ED_2010  %>% between (100,199) ~ '1', str_detect(DERIVED_SEER_COMBINED_N_2016, '^[cp]2') | DERIVED_AJCC_N_7TH_ED_2010  %>% between (200,299)~ '2',
#                              str_detect(DERIVED_SEER_COMBINED_N_2016, '^[cp]3') | DERIVED_AJCC_N_7TH_ED_2010  %>% between (300,399) ~ '3',
#                              str_detect(DERIVED_SEER_COMBINED_N_2016, '^[cp]X') | DERIVED_AJCC_N_7TH_ED_2010  == 99 ~ 'X',
#                              T ~ NA_character_ ),
#           tnm.m = case_when ( 
#                              str_detect(DERIVED_SEER_COMBINED_M_2016, '^[cp]0') | DERIVED_AJCC_M_7TH_ED_2010 %>% between (0, 10) ~ '0',
#                              str_detect(DERIVED_SEER_COMBINED_M_2016, '^[cp]1') | DERIVED_AJCC_M_7TH_ED_2010   %>% between (100,199) ~ '1',
#                              str_detect(DERIVED_SEER_COMBINED_M_2016, '^[cp]X') | DERIVED_AJCC_M_7TH_ED_2010   == 99 ~ 'X',
#                              T ~ NA_character_ ),
#           size.lt2015 = case_when ( 
#                                    CS_TUMOR_SIZE_2004_2015  >= 0 & CS_TUMOR_SIZE_2004_2015 <= 989 ~  CS_TUMOR_SIZE_2004_2015 / 10, 
#                                    CS_TUMOR_SIZE_2004_2015  == 0 ~  0, 
#                                    CS_TUMOR_SIZE_2004_2015  == 991 ~  0.99, 
#                                    CS_TUMOR_SIZE_2004_2015  == 992 ~  1.99, 
#                                    CS_TUMOR_SIZE_2004_2015  == 993 ~  2.99, 
#                                    CS_TUMOR_SIZE_2004_2015  == 994 ~  3.99, 
#                                    CS_TUMOR_SIZE_2004_2015  == 995 ~  4.99, 
#                                    T ~ NA_real_
#                                    ),
#           size.gt2015 = case_when ( 
#                                    TUMOR_SIZE_SUMMARY_2016  > 0 & TUMOR_SIZE_SUMMARY_2016  <= 989 ~ TUMOR_SIZE_SUMMARY_2016 / 10,
#                                    TUMOR_SIZE_SUMMARY_2016 == 990 ~ 0, 
#                                    T  ~ NA_real_ 
#           ),
#           size = ifelse ( nna(size.lt2015), size.lt2015, size.gt2015),
#           t_stage_8 = case_when (
#                                  size>0 & size<1 & tnm.t==1  ~ 'T1a',
#                                  size>=1 & size<2 & tnm.t==1 ~ 'T1b',
#                                  size>=2 & size<3 & tnm.t==1 ~ 'T1c',
#                                  (size>=3 & size<4) | (tnm.t=='2' & size<4) ~ 'T2a',
#                                  size>=4 & size<5 ~ 'T2b',
#                                  (size>=5 & size<7) | (tnm.t=='3' & size<7) ~ 'T3',
#                                  (size >=7 & size<100) | tnm.t=='4'  ~ 'T4',
#                                  T ~ NA_character_),
#           dx.date = ymd( ifelse ( nna(YEAR_OF_DIAGNOSIS) & nna(MONTH_OF_DIAGNOSIS) , sprintf('%d%02d15', YEAR_OF_DIAGNOSIS, MONTH_OF_DIAGNOSIS), NA_character_ ) )  ,
#           #death.date = ymd( ifelse ( ""!=(SEER_DATEOFDEATH_YEAR) & ""!=(SEER_DATEOFDEATH_MONTH) , sprintf('%s%s15', SEER_DATEOFDEATH_YEAR, SEER_DATEOFDEATH_MONTH), NA_character_ ) ),
#           cod.new = case_when(
#             grepl("C", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) & !grepl("C34", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Other Cancer',
#             grepl("C34", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Lung Cancer',
#             grepl("J", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Respiratory Disease',
#             grepl("I", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Circulatory Disease',
#             grepl("A", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("B", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Infection_Parasite',
#             grepl("E", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Endocrine Disorder',
#             grepl("F", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Mental Disorder',
#             grepl("G", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Nervous System Disease',
#             grepl("H", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Eye and ear Diseases',
#             grepl("K", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Digestive System Diseases',
#             grepl("M", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Diseases of the musculoskeletal system and connective tissue',
#             grepl("N", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Diseases of the genitourinary system',
#             grepl("V", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("W", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X1", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X2", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X3", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X4", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X5", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X9", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("Y", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'External Causes of Morbidity Except Suicide',
#             grepl("X6", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X7", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) | grepl("X8", CAUSE_OF_DEATH_ICD_10, fixed=TRUE) ~ 'Suicide',
#             CAUSE_OF_DEATH_ICD_10 == "" ~ 'Alive',
#             CAUSE_OF_DEATH_ICD_10 == '7777' | CAUSE_OF_DEATH_ICD_10 == '7797' ~ 'Unknown',
#             T ~ 'Other'
#           ),
#           Smoking = nna(smoking_pre),
#           Oxygen = nna(o2_pre) )
#    A.old  <- A.old %>% mutate( 
#                       dx.to.tx  =  as.numeric( tx.date - dx.date, units = 'days'),
#                       death.date.seer = 
#                           ymd( ifelse ( ""!=(SEER_DATEOFDEATH_YEAR) & ""!=(SEER_DATEOFDEATH_MONTH) , sprintf('%s%s15', SEER_DATEOFDEATH_YEAR, SEER_DATEOFDEATH_MONTH), NA_character_ ) )  ,
#                       tt = as.numeric( if_else ( nna(death.date.mbsf), death.date.mbsf, ymd('20191231')  ) - tx.date, units = 'days'),
#                       # time.enrolled = as.numeric( if_else ( nna(death.date.mbsf), death.date.mbsf, ymd('20191231')  ) - first.dx.date, units = 'days'),
#                       #TODO: The current time.enrolled is based on the first diagnosis. Instead, it should be based on the mbsf
#                       thirty.day.mortality = ifelse ( nna(death.date.mbsf) & tt < 30, T, F ) ,
#                       ninety.day.mortality = ifelse ( nna(death.date.mbsf) & tt < 90, T, F ) ,
#                       death = death.date.mbsf,
#                       valid.death.indicator = case_when(
#                                 is.na(death.date.seer) & is.na( death.date.mbsf)  ~ 'valid', # not death in either
#                                 nna(death.date.seer) & nna( death.date.mbsf)  ~ 'valid', # dead in both
#                                 nna(death.date.seer) & is.na( death.date.mbsf) ~ 'invalid', # Dead in SEER but not in MBSF is invalid
#                                 nna(death.date.mbsf) & is.na( death.date.seer)  & year(death.date.mbsf) == 2019 ~ 'valid', # Dead in MBSF but not in SEER is valid if it occured in 2019
#                                 nna(death.date.mbsf) & is.na( death.date.seer)  & year(death.date.mbsf) < 2019 ~ 'invalid', # Dead in MBSF but not in SEER is invalid if it occured <2019
#                                                         ))
#A.old.final  <- A.old%>% filter ( histology !="Small cell")
#print(nrow(A.old.final))
#A.old.final  <-  A.old.final %>% filter ( (t_stage_8=="T1a" | t_stage_8=="T1b" | t_stage_8=="T1c") & tnm.n==0 & tnm.m==0 )
#print(nrow(A.old.final))
# A.old.final  <-  A.old.final %>% filter (tx %in% c('sublobar', 'sbrt') )
#incex(A.old.final)
#A.old.final  <-  A.old.final %>% filter (age >= 65 ) #TODO: <80
#incex(A.old.final)
## A.old.final  <-  A.old.final %>% filter (dx.to.tx <= 120 & dx.to.tx >= -16) # The diagnosis date is chosen to be the 15th of each month, as the date itself is not available
# A.old.final  <-  A.old.final %>% filter (dx.to.tx <= 365 & dx.to.tx >= 0) # The diagnosis date is chosen to be the 15th of each month, as the date itself is not available
#incex(A.old.final)
#summary(A.old.final$dx.to.tx)
#A.old.final  <-  A.old.final %>% filter (pre.tx.months >= 12)
#incex(A.old.final)
#A.old.final  <-  A.old.final %>% filter ((valid.pet.scan & tx=='sbrt') | tx=='sublobar')
#incex(A.old.final)
#A.old.final  <-  A.old.final %>% filter (valid.death.indicator == 'valid' )
#incex(A.old.final)

#filename.out  <-  'data/A.final10.all.gte.65_old_recreate.RDS' 
# A.old.final  %>%  write_rds( filename.out)
#A %>% filter(PATIENT_ID == 'lnK2020w0096681') %>% glimpse
#A.final %>% filter(PATIENT_ID == 'lnK2020w0096681') %>% glimpse
#A.old %>% filter(PATIENT_ID == 'lnK2020w0096681')
#A.old.final %>% filter(PATIENT_ID == 'lnK2020w0096681')
#incex(A.old.final)
#incex(A.final)

## end




# Begin figuring out
# A.fofo  <-  A.final2 %>% filter ( seer.surgery == 22 & is.na(tx) & PRIMARY_PAYER_AT_DX  == 60)
# A.fofo  <-  A.final2 %>% filter (dx.to.tx < -30)
# A.fofo  <-  A.final2  %>% filter (METS_pre_count >0 )
# table( A.final2$tx, A.final2$REASONNOCANCER_DIRECTED_SURGERY, useNA="ifany")
# table( A.final2$tx, A.final2$seer.surgery, useNA="ifany")
 # A.fofo %>% filter (PATIENT_ID == 'lnK2020w0153118') %>% glimpse
 # dx.long %>% filter (PATIENT_ID == 'lnK2020w0153118') %>% filter(tx.date > CLM_THRU_DT) %>%  count(DX) %>%arrange(-n) %>% print(n=100)
# A.fofo %>% filter (PATIENT_ID == 'lnK2020w0176506') %>% glimpse
# mbsf.long %>% filter (PATIENT_ID == 'lnK2020w0134942') %>% glimpse
# End figuring out

# A.final %>% count(histology) %>% arrange(-n) %>% print(Inf)
for (x in ls()) {
    print( x);
    print(format(object.size(get(x)), unit = 'auto'))
}

sapply(ls(), function(x) print( x); print(format(object.size(get(x)), unit = 'auto'))) 
