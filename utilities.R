library(ggtext)
label_list  <-  readRDS('data/label.list.RDS')
label_list  <-  list(
                     age = 'Age',
                     sex = 'Male Sex',
                     race = 'Race',
                     size = 'Tumor size (cm)',
                     marital.status = 'Marital status',
                     histology =  'Histology',
                     histology2 =  'Squamous cell histology',
                     t_stage_8 = 'T stage 8th edition',
                     treatment.year = 'Treatment year',
                     treatment.year2 = 'Treatment year',
                     smoking_pre_month_count_unbinned = 'Smoking',
                     o2_pre_month_count_unbinned = 'Oxygen',
                     pneumonia_and_influenza_pre_month_count_unbinned = 'Pneumonia and influenza',
                     other_bacterial_diseases_pre_month_count_unbinned = 'Other bacterial diseases',
                     pressure_ulcer_pre_month_count_unbinned = 'Pressure ulcer',
                     ischemic_heart_disease_pre_month_count_unbinned = 'Ischemic heart disease',
                     CHF_pre_month_count_unbinned = 'Congestive heart failure',
                     other_heart_disease_pre_month_count_unbinned = 'Other heart disease',
                     PVD_pre_month_count_unbinned = 'Peripheral vascular disease',
                     CVD_pre_month_count_unbinned = 'Cerebrovascular disease',
                     dementia_pre_month_count_unbinned = 'Dementia',
                     COPD_pre_month_count_unbinned = 'Chronic obstructive pulmonary disease',
                     asthma_pre_month_count_unbinned = 'Asthma',
                     interstitial_lung_pre_month_count_unbinned = 'Interstitial lung disease',
                     other_lung_pre_month_count_unbinned = 'Other lung disease',
                     PUD_pre_month_count_unbinned = 'Peptic ulcer disease',
                     MILDLD_pre_month_count_unbinned = 'Mild liver disease',
                     MSLD_pre_month_count_unbinned = 'Moderate or severe liver disease',
                     DIAB_UC_pre_month_count_unbinned = 'Uncomplicated diabetes',
                     DIAB_C_pre_month_count_unbinned = 'Complicated diabetes',
                     # PARA_pre_month_count_unbinned = 'Paralysis',
                     RD_pre_month_count_unbinned = 'Renal disease',
                     # cancer_nonlung_pre_month_count_unbinned = 'Cancer (non-lung)',
                     # METS_pre_month_count_unbinned = 'Metastatic disease',
                     mental_disorders_pre_month_count_unbinned = 'Mental disorders',
                     nervous_system_pre_month_count_unbinned = 'Neurological disorders',
                     # veins_lymphatics_other_circulatory_pre_month_count_unbinned = 'Veins, lymphatics, other circulatory',
                     dialysis_pre_month_count_unbinned = 'Dialysis',
                     echo_pre_month_count_unbinned = 'Echocardiogram',
                     # rheum_pre_month_count_unbinned = 'Rheumatologic diseases',
                     Insulin_pre_month_count_unbinned = 'Insulin',
                     Anticoags_pre_month_count_unbinned = 'Anticoagulation',
                     # Z
                     O2accessories_pre_month_count_unbinned = 'Oxygen accessories',
                     walking_aids_pre_month_count_unbinned = 'Walking aids',
                     hospital_beds_and_supplies_pre_month_count_unbinned = 'Hospital beds and supplies',
                     wheelchairs_accessories_pre_month_count_unbinned = 'Wheelchairs and accessories',
                     transportation_services_pre_month_count_unbinned = 'Transportation services',
                     other_supplies_pre_month_count_unbinned = 'Other supplies',
                     diabetic_footwear_pre_month_count_unbinned = 'Diabetic footwear',
                     # W
                     fall_pre_month_count_unbinned = 'Fall',
                     other_injury_pre_month_count_unbinned = 'Other injury',
                     diverticular_disease_pre_month_count_unbinned = 'Diverticular disease',
                     hernia_pre_month_count_unbinned = 'Hernia',
                     arthropathy_pre_month_count_unbinned = 'Arthropathy',
                     GU_sx_pre_month_count_unbinned = 'Genitourinary symptoms',
                     optho2_pre_month_count_unbinned = 'Ophthalmologic disease'
                     # other.cause.mortality = 'Other cause mortality',
                     # cause.specific.mortality = 'Cause specific mortality',
                     # primary.site = 'Primary site',
                     # histology = 'Histology',
                     # YEAR_OF_DIAGNOSIS = 'Year of Diagnosis',
                     # BEHAVIOR_CODE_ICD_O_3 = 'Behavior',
                     # death      = 'Death',
                     # thirty.day.mortality = '30-day mortality',
                     # ninety.day.mortality = '90-day mortality',
                     # cod.new = 'Cause of Death Category'
)
label_list2  <-  c( label_list,
                   death = 'Overall mortality', 
                   death.cause.specific = 'Cancer-specific mortality', 
                   death.other.cause = 'Other-cause mortality', 
                   death.other.cause.gt90day = 'Other mortality (>90 days)', 
                   death.90.day = '90-day mortality', 
                   fall = 'Fall',
                   other_injury = 'Injury',
                   GU_sx = 'GU-related',
                   arthropathy = 'Arthropathy',
                   cholelithiasis = 'Cholelithiasis-related',
                   gout = 'Gout',
                   obstruction = 'Intestinal obstruction',
                   hernia = 'Abdominal hernia',
                   diverticular_disease = 'Diverticular disease',
                   hemorrhoids = 'Hemorrhoids',
                   pancreatic = 'Pancreatic',
                   optho2 = 'Ophthalmic',
                   oral = 'Oral'
)
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
    xlims <- c(0.3, 3)
    tt  <-  ifelse (hazard, 'Hazard Ratio (log scale)', 'Hazard Ratio (log scale)')
    # row.names(odds.ratios_) <- gsub("_pre_month_count", "", row.names(odds.ratios_))
     row.names(odds.ratios_) <- gsub("_..._count", "",  row.names(odds.ratios_))
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

make.HD.plot  <-  function (odds.ratios_, label_list2, hazard =F, xlims = c(-0.05,0.15), breaks=c(-0.05,0.025,0, 0.05,0.1, 0.15)) {
    tt  <-  'Hazard Difference'
    g <- ggplot(odds.ratios_, aes(x = estimate, y=y_axis)) + 
        geom_vline(aes(xintercept = 0), size = 0.25, linetype = "dashed") +
        geom_errorbarh(aes( xmax = high_ci, xmin = low_ci), size = 0.20, height = 0.3)+
        geom_point(size=1.5) +
        theme_bw() +
        theme(panel.grid.minor = element_blank()) +
        scale_y_continuous(breaks = 1:max(odds.ratios_$y_axis), labels = label_list2[row.names(odds.ratios_)], trans='reverse') +
        scale_x_continuous(limits = xlims, breaks = breaks ) +
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
scale_  <-  function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T)


################################
# Backup 
################################
# system.time ( MBSF_2013 <- read.SAScii('../SEER-Medicare-data/data/SEER_Medicare/mbsf.abcd.summary.2013.txt', 
#             '../SEER-Medicare-data/data/SEER_Medicare/2020 Input Statements/MBSF.ABCD.gcl.txt', n=100000) ) 

# lung.SEER <- read.SAScii('../SEER-Medicare-data/data/SEER_Medicare/SEER.lung.cancer.n1000.txt',  
                         # '../SEER-Medicare-data/data/SEER_Medicare/2020 Input Statements/SEER.input.gcl.txt', n=100) 

