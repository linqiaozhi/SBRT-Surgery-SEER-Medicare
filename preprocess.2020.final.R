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
# II. Identify the treatment received by patients using the medicare files. Output: patient.tx
#       1. Inpatients: Load the Medpar file, and check for SBRT or Resection for inpatient
#       2. Outpatients: The outpatient file contains bills from institutions (ie and not inpatient stays). The carrier file contains bills from providers. 
#             1. Carrier line files contains both procedure and the daignosis for which it was performed. We also load the carrier base file for later diagnosis extraction.
#             2.  Outpatient Revenue file contains procedures, link it with Outpatient Line file to confirm diagnosis for each procedure is lung cancer.
#       3. Combine everything into a single dataframe patient.tx, with three two columns: PATIENT_ID, tx.date, tx
# III. MBSF processing. Output: patient.mbsf
# IV. Identify diagnoses (comorbidities and negative control outcomes). Output: patient.dx
#       1. Combine long versions of medpar, outpatient line, and carrier base. Filter each by patient.tx. Results in dataframes into all.long.icd9.dx, all.long.icd10.dx
#       2. Identify comorbidities of interest, patient.comorbidities
#       3. Identify negative outcome diagnoses of interest, patient.noc
# V. Identify procedure and DME codes (for now, just PET). Output: patient.proc
# VI.  Identify medications using Part D files 
# VIII. Combine patient.tx,  patient.dx, patient.proc, patient.mbsf.death, patient.seer
# IX Massage variables
# X. Filter





################################
#  SECTION I: Load SEER file
################################

lung.SEER <- read_dta(sprintf('%s/seerlung.%s.dta.gz', dta.path, suffix))
# This filtering is just to speed things up.
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
# Only valid if the diagnosis for each is actually lung cancer. Note that this includes nodules.
medpar$sbrt.date[ ! medpar$actually.lung.cancer ]  <- as.Date(NA_Date_)
medpar$sublobar.date[ ! medpar$actually.lung.cancer ]  <- as.Date(NA_Date_)
medpar$lobar.date[ ! medpar$actually.lung.cancer ]  <- as.Date(NA_Date_)
medpar$other.resection[ ! medpar$actually.lung.cancer ]  <- as.Date(NA_Date_)

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
  carrierbasei.small  <- carrierbasei.small %>% distinct()
  carrierbases[[year]]  <-  carrierbasei.small
}
carrierbase  <-  bind_rows(carrierbases,  .id='dataset.year')
rm(carrierbases); gc();

# Carrier line file contains procedures and their diagnosis codes
fn  <-  sprintf("%s/nch.lines.RDS", rds.path)
if (!file.exists( fn )) {
    years  <-  as.character(2010:2020)
    carriers  <-  list()
    for (yeari in 1:length(years)) {
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
carrier  <- carrier %>% distinct()

carrier <- carrier %>% mutate(
              valid.dx = LINE_ICD_DGNS_CD %in% valid.dxs ,
              sbrt =   HCPCS_CD %in% sbrt.cpts  & valid.dx,
              sbrt.date  = if_else ( sbrt, CLM_THRU_DT, NA_character_) %>% ymd ,
              ebus =   HCPCS_CD %in% ebus.cpts,
              ebus.date  = if_else ( ebus, CLM_THRU_DT, NA_character_) %>% ymd ,
              med =   HCPCS_CD %in% med.cpts,
              med.date  = if_else ( med, CLM_THRU_DT, NA_character_) %>% ymd ,
              sublobar =   HCPCS_CD %in% sublobar.cpts  & valid.dx,
              sublobar.date  = if_else ( sublobar, CLM_THRU_DT, NA_character_) %>% ymd ,
              lobar =   HCPCS_CD %in% lobar.cpts  & valid.dx,
              lobar.date  = if_else ( lobar, CLM_THRU_DT, NA_character_) %>% ymd ,
              other.resection =   HCPCS_CD %in% other.resection.cpts  & valid.dx,
              other.resection.date  = if_else ( other.resection, CLM_THRU_DT, NA_character_) %>% ymd 
)

table( nna(carrier$ebus.date), useNA="ifany")

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


outpat.outpat.revenue  <-  outpat  %>% inner_join( outpat.revenue, by = c('PATIENT_ID', 'CLM_ID')) 
outpat.outpat.revenue <- outpat.outpat.revenue %>% mutate(
    valid.dx  =   PRNCPAL_DGNS_CD %in%  valid.dxs,
    sbrt  =   HCPCS_CD  %in%  sbrt.cpts  & valid.dx,
    sbrt.date  =  if_else ( sbrt, CLM_THRU_DT.x, as.Date(NA_Date_) ),
      ebus =   HCPCS_CD %in% ebus.cpts,
      ebus.date  = if_else ( ebus, CLM_THRU_DT.x, as.Date(NA_Date_)) %>% ymd ,
)



################################
# II.3 Combine all three sources of treatment codes
################################
# Filter each for efficiency; there are a LOT of rows we do not need here
medpar.tx  <-   medpar %>% select ( PATIENT_ID, sbrt.date, sublobar.date, lobar.date, other.resection.date) %>% filter(!is.na(sbrt.date) | !is.na(sublobar.date) | ! is.na(lobar.date) )
carrier.tx  <-  carrier %>% select( PATIENT_ID, sbrt.date, sublobar.date, lobar.date, other.resection.date, ebus.date, med.date) %>% filter(!is.na(sbrt.date) | !is.na(sublobar.date) | ! is.na(lobar.date) | !is.na(other.resection.date) | !is.na(ebus.date)  )
outpat.tx  <-  outpat.outpat.revenue %>% select(PATIENT_ID, sbrt.date, ebus.date) %>% filter(!is.na(sbrt.date) | !is.na(ebus.date))
# This join filters to only surgeries which match the type and date of surgery between both medpar and carrier files
# surgery.tx  <- medpar.tx %>% inner_join (carrier.tx, by = c('PATIENT_ID', 'sublobar.date', 'lobar.date', 'other.resection.date')) %>% select(PATIENT_ID, sublobar.date, lobar.date, other.resection.date)
surgery.tx  <- medpar.tx %>% full_join (carrier.tx, by = c('PATIENT_ID'), suffix = c('.me', '.ca')) 

sbrt.tx  <- bind_rows ( medpar.tx, carrier.tx, outpat.tx) %>% arrange(PATIENT_ID, sbrt.date, ebus.date) %>% select(PATIENT_ID, sbrt.date, ebus.date) %>% filter (nna(sbrt.date)| nna(ebus.date))
seer.tx  <- patient.seer %>% mutate (    seer.sublobar  = RX_SUMM_SURG_PRIM_SITE_1998 %in% c("20", "21", "22", "23" ),
                                         seer.lobar     = RX_SUMM_SURG_PRIM_SITE_1998 %in% c("30", "32", "33", "45", "46", "47", "48"),
                                         seer.other.resection     = RX_SUMM_SURG_PRIM_SITE_1998 %in% c("55", "65", "66") ) 

# The patient's treatment is specified by the SEER variable, and the
# medicare code tells when it occured. Both of these need to be present: SEER
# code specifying the type of surgery and a Medpar or Carrier code specifying
# the same type of surgery. The Medpar code is for the instution, the carrier
# is for the surgeon.
# Similarly, for SBRT there needs to a corresponding code in the SEER file and
# the Carrier file. Medpar is not relevant as this is not an inpatient
# procedure.
patient.tx  <- seer.tx %>% 
    left_join (surgery.tx, by = 'PATIENT_ID')  %>%
    left_join(sbrt.tx, by = 'PATIENT_ID') %>% 
     filter(nna(sbrt.date) | seer.sublobar + seer.lobar + seer.other.resection > 0) %>% 
    group_by(PATIENT_ID) %>% 
    summarise (  
               tx = factor( case_when ( 
                                       all(seer.sublobar) &  ( any(nna ( sublobar.date.me)) | any(nna ( sublobar.date.ca)) )     ~ 'sublobar',
                                       all(seer.lobar) & ( any(nna ( lobar.date.me)) | any(nna ( lobar.date.ca)) )  ~ 'lobar',
                                       all(seer.other.resection) &( any(nna ( other.resection.date.me)) | any(nna ( other.resection.date.ca)) ) ~ 'other.resection',
                                        any(nna ( sbrt.date))  ~ 'sbrt',
                                       T ~ (NA_character_)
                                       ), levels = c('sublobar', 'sbrt', 'lobar', 'other.resection')),
               sublobar.date = min(min (na.omit(sublobar.date.me)),min (na.omit(sublobar.date.ca)), na.rm = T),
               lobar.date = min(min (na.omit(lobar.date.me)),min (na.omit(lobar.date.ca)), na.rm = T),
               sbrt.date = min(na.omit(sbrt.date)),
               tx.date = case_when (
                                    tx == 'sbrt' ~ sbrt.date,
                                    tx == 'sublobar' ~ sublobar.date,
                                    tx == 'lobar' ~ lobar.date,
                                    T ~ ymd(NA_character_)
                                    )
               )


gc()
################################
# SECTION III MBSF processing
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
# The issue is that only months where patient is enrolled as FFS count, because
# those are the only motnhs where all billing codes are collected. We need two
# things: the number of pre-treatment FFS months (pre.tx.months), and the number of continuous
# months immediately preceeding the treatment where the patient was enrolled
# FFS (ffs.months.counter)
mbsf.wide  <- mbsf  %>%  select( PATIENT_ID,BENE_ENROLLMT_REF_YR,  MDCR_STATUS_CODE_01:MDCR_STATUS_CODE_12,HMO_IND_01:HMO_IND_12)  %>% rename_with(
                                     .fn = ~ gsub("MDCR_STATUS_CODE", "MDCRSTATUSCODE", .x),
                                     .cols = starts_with("DCR_STATUS_CODE")
                                     ) %>%
                            rename_with(
                                        .fn = ~ gsub("HMO_IND", "HMOIND", .x),
                                        .cols = starts_with("HMO_IND")
                            )
sum_reset   <- function(x) accumulate( x, ~(.y!=0)*(.y+.x))
year.month  <-  function(x) {
    sprintf('%d-%d', year(x), month(x) )
}



# colnames(mbsf.wide)
#  mbsf.wide %>% select( -c(PATIENT_ID, BENE_ENROLLMT_REF_YR), -contains('MDCR')) 
#  mbsf.wide %>% select( -) 
# fofo %>% print(width=Inf)
# fofo  <-  mbsf.wide[1:100,] %>% pivot_longer(cols = c(-c(PATIENT_ID, BENE_ENROLLMT_REF_YR), -contains('MDCR')), names_to = c(".value", "month" ), names_sep="_" )
# fofo  <-  mbsf.wide[1:100,] %>% select(-contains('MDCR') ) %>%  pivot_longer(cols = contains('HMO'), names_to = c(".value", "month" ), names_sep="_"  )
# fofo %>% print(width=Inf)

mbsf.long.ffs  <- mbsf.wide %>% 
    select(-contains('MDCR') ) %>%  pivot_longer(cols = contains('HMO'), names_to = c(".value", "month" ), names_sep="_"  ) %>% 
    # pivot_longer(cols = -c(PATIENT_ID, BENE_ENROLLMT_REF_YR), names_to = c(".value", "month" ), names_sep="_" ) %>%
    mutate(year = BENE_ENROLLMT_REF_YR,
           # Last date of the month. Will use this later to count the number of months where tx.date is <= date for FFS months
           date = ceiling_date(ymd(sprintf('%d-%s-01', year, month)), 'month') -1 ,
           ym = year.month(date),
           ffs.month = HMOIND %in% c(0,4) ) 
mbsf.long.ffs.valid  <-  mbsf.long.ffs   %>%  filter (ffs.month)

    
mbsf.long.ffs.grouped  <- mbsf.long.ffs %>% 
    arrange(PATIENT_ID, date) %>% 
    group_by(PATIENT_ID) %>% 
    mutate(ffs.months.counter = sum_reset(1*ffs.month)) %>% 
    ungroup()

month.year  <- function(x)  sprintf('%s-%s', month(x), year(x))
mbsf.long.ffs.grouped <- mbsf.long.ffs.grouped %>% mutate(  month.year = month.year(date) )
# month.year corresponds to tx.date, so it will pull the counter for that month
patient.mbsf.months.counter  <- patient.tx %>% select( PATIENT_ID, tx.date) %>% 
                                    mutate ( month.year = month.year(tx.date)) %>%
                                left_join(mbsf.long.ffs.grouped %>% select(PATIENT_ID, month.year, ffs.months.counter)) %>%
                                select(PATIENT_ID, ffs.months.counter)



patient.mbsf.total.ffs <- mbsf.long.ffs.grouped %>% 
    right_join(patient.tx %>% select( PATIENT_ID, tx.date)) %>% 
    group_by(PATIENT_ID) %>% 
    summarise( 
              pre.tx.months = sum(ffs.month & (date <= tx.date)),
              total.ffs.months = sum(ffs.month))

patient.mbsf <- mbsf %>% mutate( 
                        death.date.mbsf.temp = ifelse(mbsf$BENE_DEATH_DT!= "" & mbsf$BENE_DEATH_DT != "       .", mbsf$BENE_DEATH_DT, NA_Date_) %>% ymd ) %>% 
                        group_by( PATIENT_ID) %>% 
                        summarise( 
                                  death.date.mbsf = first(na.omit(death.date.mbsf.temp))
                            ) %>% 
                        right_join(patient.mbsf.months.counter) %>% 
                        right_join(patient.mbsf.total.ffs)

gc()
################################
# Section IV:  Identify diagnoses 
################################

# Create the long diagnosis data frame
# We join each one with patient.tx for efficiency, since we only are studying a
# small fraction of these patients
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
rm(medpar); rm(outpat);  gc()

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

# Filter to only include diagnoses assigned during FFS months
dx.long  <- dx.long %>%  
    mutate ( ym = year.month(CLM_THRU_DT) ) %>%
    inner_join( mbsf.long.ffs.valid )
rm(dx.wide); gc()


dx.hardcodeds  <- patient.tx %>% select(PATIENT_ID)
for (i in 1:length(dx.icd)) {
    print(sprintf('%d/%d hard coded dx', i, length(dx.icd)))
    dxoi  <- dx.icd[[i]]
    dx.name  <-  names(dx.icd)[i]
    dx.hardcoded  <- dx.long %>% 
             mutate( 
                    temp = if_else ( icd9or10 == 'icd9', DX %in% dxoi$icd9,DX %in% dxoi$icd10 ),
                    temp.pre =  if_else(temp & (CLM_THRU_DT < tx.date), CLM_THRU_DT, ymd(NA_character_)), 
                    temp.pre.12months =  if_else(temp & (CLM_THRU_DT < tx.date & CLM_THRU_DT >= (tx.date - months(12))), CLM_THRU_DT, ymd(NA_character_)), 
                    temp.post =  if_else(temp & (CLM_THRU_DT > tx.date), CLM_THRU_DT, ymd(NA_character_))  ,
                    temp.any =  if_else(temp , CLM_THRU_DT, ymd(NA_character_))  
                    )  %>% 
             group_by(PATIENT_ID) %>% 
             summarise( 
                          !!dx.name := first(na.omit(temp.post)), 
                          !!sprintf('%s_pre_count', dx.name ) := length((na.omit(temp.pre))),
                          !!sprintf('%s_pre_month_count', dx.name ) := length(unique(year.month(na.omit(temp.pre)))),
                          !!sprintf('%s_pre_12months_count', dx.name ) := length((na.omit(temp.pre.12months))),
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
# SECTION V Procedure and DME codes (not for treatment)
################################
print('SECTION V Procedure and DME codes')
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
rm(carrier);  gc()

# Filter to only include diagnoses assigned during FFS months
procs.long  <- procs.long %>%  
    mutate ( ym = year.month(CLM_THRU_DT) ) %>%
    inner_join( mbsf.long.ffs.valid )

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
                     temp.pre.12months =  if_else(temp & (CLM_THRU_DT < tx.date & CLM_THRU_DT >= (tx.date - months(12))), CLM_THRU_DT, ymd(NA_character_)), 
                    temp.post =  if_else(temp & (CLM_THRU_DT > tx.date), CLM_THRU_DT, ymd(NA_character_))  ,
                    temp.any =  if_else(temp , CLM_THRU_DT, ymd(NA_character_))  
                    )  %>% 
             group_by(PATIENT_ID) %>% 
             summarise( 
                          !!proc.name := first(na.omit(temp.post)), 
                          !!sprintf('%s_pre_count', proc.name ) := length((na.omit(temp.pre))),
                           !!sprintf('%s_pre_12months_count', proc.name ) := length((na.omit(temp.pre.12months))),
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
# SECTION VI Part D 
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


# # Filter to only include diagnoses assigned during FFS months

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
      temp.pre.12months =  if_else(temp & (SRVC_DT < tx.date & SRVC_DT >= (tx.date - months(12))), SRVC_DT, ymd(NA_character_)), 
      temp.post = if_else(temp & (SRVC_DT > tx.date), SRVC_DT, ymd(NA_character_)),
      temp.any = if_else(temp, SRVC_DT, ymd(NA_character_))
    ) %>% 
    group_by(PATIENT_ID) %>% 
    summarise(
      !!drug.name := first(na.omit(temp.post)),
      !!sprintf('%s_pre_count', drug.name) := length((na.omit(temp.pre))),
      !!sprintf('%s_pre_12months_count', drug.name ) := length((na.omit(temp.pre.12months))),
      !!sprintf('%s_pre_month_count', drug.name ) := length(unique(year.month(na.omit(temp.pre)))),
      !!sprintf('%s_any_count', drug.name ) := length((na.omit(temp.any))),
      !!sprintf('%s_any_month_count', drug.name ) := length(unique(year.month(na.omit(temp.any)))),
      !!sprintf('%s_post_count', drug.name) := length((na.omit(temp.post))),
      !!sprintf('%s_post_month_count', drug.name ) := length(unique(year.month(na.omit(temp.post)))),
    )
  patient.drugs <- patient.drugs %>% left_join(drug, by = 'PATIENT_ID')
}


################################
# Section VIII Combine 
################################
topography  <-  read_csv(file= './ICDO3topography.csv') %>% rename(site.topography = description) %>% mutate(PRIMARY_SITE = str_remove_all( icdo3_code, fixed(".")))
# Note that each of the patient.* data frames were filtered to only include
# patients who underwent treatment (i.e. filtered by patient.tx). Thus SEER
# includes many more patients than those in the patient.* data frames, so many
# fields will be NA in this join, simply because those patients were not
# preprocessed in the above code.
A  <-  patient.seer %>% 
    left_join( patient.dx, by = 'PATIENT_ID') %>% 
    left_join( patient.outpatient.procs, by = 'PATIENT_ID',) %>% 
    left_join( patient.mbsf, by = 'PATIENT_ID') %>%
     left_join( patient.tx , by = 'PATIENT_ID') %>% 
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
           seer.surgery = case_when (
                                        RX_SUMM_SURG_PRIM_SITE_1998 %in% c("20", "21", "22", "23" ) ~ 'sublobar',
                                        RX_SUMM_SURG_PRIM_SITE_1998 %in% c("30", "32", "33", "45", "46", "47", "48") ~ 'lobar',
                                        RX_SUMM_SURG_PRIM_SITE_1998 %in% c("55", "65", "66") ~ 'other_resection',
                                        # RADIATION_RECODE %in% c("1", "5") ~ 'radiation',
                                        RX_SUMM_SURG_PRIM_SITE_1998 %in% c("00") ~ 'no_surgery',
                                         T ~ NA_character_ ),
           age = as.numeric(age),
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
           dx.date =  ymd(ifelse ( YEAR_OF_DIAGNOSIS != "" & MONTH_OF_DIAGNOSIS != "" , sprintf('%s%s15', YEAR_OF_DIAGNOSIS, MONTH_OF_DIAGNOSIS), NA_character_ ) )  ,
           # tx.date.seer =  ymd(ifelse ( YEAR_THERAPY_STARTED != "" & MONTH_THERAPY_STARTED != "" , sprintf('%s%s15', YEAR_OF_DIAGNOSIS, MONTH_OF_DIAGNOSIS), NA_character_ ) )  ,
           # tx.month.year.seer =  month.year(tx.date.seer)  ,
                       dx.to.tx  =  as.numeric( tx.date - dx.date, units = 'days'),
                       death.date.seer = 
                           ymd( ifelse ( ""!=(SEER_DATEOFDEATH_YEAR) & ""!=(SEER_DATEOFDEATH_MONTH) , sprintf('%s%s15', SEER_DATEOFDEATH_YEAR, SEER_DATEOFDEATH_MONTH), NA_character_ ) )  ,
                       tt = as.numeric( if_else ( nna(death.date.mbsf), death.date.mbsf, ymd('20201231')  ) - tx.date, units = 'days'),
                       thirty.day.mortality = ifelse ( nna(death.date.mbsf) & tt < 30, T, F ) ,
                       ninety.day.mortality = ifelse ( nna(death.date.mbsf) & tt < 90, T, F ) ,
                       death = death.date.mbsf, 
                       valid.death.indicator = is.na (death.date.seer) == is.na(death.date.mbsf) & cause.specific.mortality != 'Unknown' & other.cause.mortality != 'Unkown')

################################
# Section IX  Exclusion
################################
# A large number of SEER patients underwent surgery based on the seer.surgery variable which are not included, as they were note enrolled in Medicare at the time of treatment. Recall that SEER includes non Medicare patietns as well. 
incex  <-  function ( A.frame ) {
    # print(sprintf('%d (SR: %d, SBRT: %d)', nrow(A.frame), sum(A.frame$tx == 'sublobar' & A.frame$seer.surgery == 'sublobar' ,na.rm = T), sum(A.frame$tx == 'sbrt' & A.frame$seer.surgery == 'no_surgery', na.rm =T)))
     # print(sprintf('%d (SR: %d, SBRT: %d)', nrow(A.frame), sum(A.frame$tx == 'sublobar' & A.frame$seer.surgery == 'sublobar' ,na.rm = T), sum(A.frame$tx == 'sbrt' & A.frame$seer.surgery == 'no_surgery', na.rm =T)))
     print(sprintf('%d (SR: %d, SBRT: %d)', nrow(A.frame), sum(A.frame$tx == 'sublobar' ,na.rm = T), sum(A.frame$tx == 'sbrt' , na.rm =T)))
}
print(nrow(A))
A.final  <- A %>% filter ( histology !="Small cell" & histology != 'Other' & histology != 'Carcinoid' & histology != 'Adenosquamous' & histology != 'Large cell' & histology != 'Non-small cell carcinoma' )
incex(A.final)
A.final  <-  A.final %>% filter ( (t_stage_8=="T1a" | t_stage_8=="T1b" | t_stage_8=="T1c") & tnm.n=='0'& tnm.m=='0')
A.final <- A.final %>% filter( RX_SUMM_SYSTEMIC_SURG_SEQ == "0")
incex(A.final)
A.final  <-  A.final %>% filter (age >= 65  & age < 90)
incex(A.final)
A.final  <-  A.final %>% filter ( ffs.months.counter >= 12 )
incex(A.final)
A.final <- A.final %>% filter ( tx %in% c( 'sbrt', 'sublobar') )
A.final$tx  <- droplevels(A.final$tx)
incex(A.final)
A.final  <- A.final %>% filter (tx == 'sublobar' | 
                                ( tx == 'sbrt' &  
                                    ( is.na(lobar.date)   |  (  lobar.date  >  sbrt.date ) ) & 
                                    ( is.na(sublobar.date)   |  (  sublobar.date  >  sbrt.date ) ) 
                                )
                            )
# A.final %>% group_by(tx,year(tx.date)) %>% summarise(n=n()) %>% print (n=Inf)
incex(A.final)
A.final  <-  A.final %>% filter (dx.to.tx <= 135 & dx.to.tx >= -16) # The diagnosis date is chosen to be the 15th of each month, as the date itself is not available
incex(A.final)
A.final  <-  A.final %>% filter (valid.death.indicator )
incex(A.final)
A.final  <-  A.final %>% filter (valid.pet.scan)
incex(A.final)
A.final  <-  A.final %>% filter (microscopically_confirmed)
incex(A.final)
A.final  %>%  write_rds( 'data/A.final13.all.gte.65.RDS' )

A.final %>% filter (tx == 'sublobar' ) %>% count (seer.surgery)
A.final %>% filter (tx == 'sublobar') %>% count(RX_SUMM_SURG_PRIM_SITE_1998)


A.final %>% filter (nna(cause.specific.mortality)) %>% count(  COD_TO_SITE_RECODE) %>% arrange(-n)


################################
# Section IX  sbrt.nodes
################################
# A large number of SEER patients underwent surgery based on the seer.surgery variable which are not included, as they were note enrolled in Medicare at the time of treatment. Recall that SEER includes non Medicare patietns as well. 
print(nrow(A))
A.sens1  <- A %>% filter ( histology !="Small cell" & histology != 'Other' & histology != 'Carcinoid' & histology != 'Adenosquamous' & histology != 'Large cell' & histology != 'Non-small cell carcinoma' )
incex(A.sens1)
A.sens1  <-  A.sens1 %>% filter (( tx == 'sbrt'  & ( (t_stage_8=="T1a" | t_stage_8=="T1b" | t_stage_8=="T1c")  & tnm.m ==0 & tnm.n =='0')) | 
                                 (tx == 'sublobar'  & ( (t_stage_8=="T1a" | t_stage_8=="T1b" | t_stage_8=="T1c")  & tnm.m ==0 ))
                             )
A.sens1$tnm.n[is.na(A.sens1$tnm.n)]  <- 'X'
 incex(A.sens1)
A.sens1 <- A.sens1 %>% filter( tx == 'sublobar' | (tx == 'sbrt' & RX_SUMM_SYSTEMIC_SURG_SEQ == "0"))
A.sens1  <- A.sens1 %>% filter (tx == 'sublobar' | 
                                ( tx == 'sbrt' &  
                                    ( is.na(lobar.date)   |  (  lobar.date  >  sbrt.date ) ) & 
                                    ( is.na(sublobar.date)   |  (  sublobar.date  >  sbrt.date ) ) 
                                )
                            )
incex(A.sens1)
A.sens1  <-  A.sens1 %>% filter (age >= 65  & age < 90)
incex(A.sens1)
A.sens1  <-  A.sens1 %>% filter ( ffs.months.counter >= 12)
incex(A.sens1)
A.sens1  <-  A.sens1 %>% filter (dx.to.tx <= 135 & dx.to.tx >= -16) # The diagnosis date is chosen to be the 15th of each month, as the date itself is not available
A.sens1  <-  A.sens1 %>% filter (valid.death.indicator )
incex(A.sens1)
A.sens1  <-  A.sens1 %>% filter (valid.pet.scan)
incex(A.sens1)
A.sens1  <-  A.sens1 %>% filter (microscopically_confirmed)
incex(A.sens1)
A.sens1  %>%  write_rds( 'data/A.final13.sens1.RDS' )
table( A.sens1$tnm.n, useNA="ifany")


table( A.final$cause.specific.mortality, useNA="ifany")

dx.icd['o2']
procs$O2accessories
