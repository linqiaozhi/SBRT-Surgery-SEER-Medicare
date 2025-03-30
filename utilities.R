library(ggtext)
label_list  <-  readRDS('data/label.list.RDS')
cod.df <- tribble(
  ~COD_TO_SITE_RECODE,   ~Name,
  '00000', 'Alive',
  '22030', 'Lung cancer death',
  '50130', 'COPD',
  '50060', 'Cardiac',
  '50080', 'Cerebrovascular',
  '50300', 'Other non-lung ancer',
  # '50120', 'Pneumonia'
)

label_list  <-  list(
                     death.stroke = '*CVD-specific mortality*',
                     death.heart = '*Heart-specific mortality*',
                     death.copd = '*COPD-specific mortality*',
                     death.noncopd.nonheart = '*Non-heart or COPD-cause mortality*',
                     death.stroke = '*Stroke-specific mortality*',
                     death.other = '*Other-specific mortality*',
                     death.noncopd.nonheart.nonstroke = '*Other non-cancer mortality*',
                     age = 'Age',
                     sex = 'Male Sex',
                     race = 'Race',
                     race2 = 'Race',
                     size = 'Tumor size (cm)',
                     marital.status = 'Marital status',
                     histology =  'Histology',
                     histology2 =  'Squamous cell histology',
                     t_stage_8 = 'T stage 8th edition',
                     treatment.year = 'Treatment year',
                     treatment.year2 = 'Treatment year',
                     smoking_pre_12months_count_bool = 'Smoking',
                     o2_pre_12months_count_bool = 'Oxygen',
                     pneumonia_and_influenza_pre_12months_count_bool = 'Pneumonia and influenza',
                     other_bacterial_diseases_pre_12months_count_bool = 'Other bacterial diseases',
                     pressure_ulcer_pre_12months_count_bool = 'Pressure ulcer',
                     ischemic_heart_disease_pre_12months_count_bool = 'Ischemic heart disease',
                     CHF_pre_12months_count_bool = 'Congestive heart failure',
                     other_heart_disease_pre_12months_count_bool = 'Other heart disease',
                     PVD_pre_12months_count_bool = 'Peripheral vascular disease',
                     CVD_pre_12months_count_bool = 'Cerebrovascular disease',
                     dementia_pre_12months_count_bool = 'Dementia',
                     COPD_pre_12months_count_bool = 'Other chronic pulmonary disease, including COPD',
                     asthma_pre_12months_count_bool = 'Asthma',
                     interstitial_lung_pre_12months_count_bool = 'Interstitial lung disease',
                     other_lung_pre_12months_count_bool = 'Other lung disease',
                     PUD_pre_12months_count_bool = 'Peptic ulcer disease',
                     MILDLD_pre_12months_count_bool = 'Mild liver disease',
                     MSLD_pre_12months_count_bool = 'Moderate or severe liver disease',
                     LD_pre_12months_count_bool = 'Liver disease',
                     DIAB_UC_pre_12months_count_bool = 'Uncomplicated diabetes',
                     DIAB_C_pre_12months_count_bool = 'Complicated diabetes',
                     # PARA_pre_12months_count_bool = 'Paralysis',
                     RD_pre_12months_count_bool = 'Renal disease',
                     # cancer_nonlung_pre_12months_count_bool = 'Cancer (non-lung)',
                     # METS_pre_12months_count_bool = 'Metastatic disease',
                     mental_disorders_pre_12months_count_bool = 'Mental disorders',
                     nervous_system_pre_12months_count_bool = 'Neurological disorders',
                     # veins_lymphatics_other_circulatory_pre_12months_count_bool = 'Veins, lymphatics, other circulatory',
                     mobility_aids_pre_12months_count_bool = 'Motility aids',
                     dialysis_pre_12months_count_bool = 'Dialysis',
                     echo_pre_12months_count_bool = 'Echocardiogram',
                     # rheum_pre_12months_count_bool = 'Rheumatologic diseases',
                     Insulin_pre_12months_count_bool = 'Insulin',
                     Anticoags_pre_12months_count_bool = 'Anticoagulation',
                     # Z
                     O2accessories_pre_12months_count_bool = 'Oxygen accessories',
                     walking_aids_pre_12months_count_bool = 'Walking aids',
                     hospital_beds_and_supplies_pre_12months_count_bool = 'Hospital beds and supplies',
                     wheelchairs_accessories_pre_12months_count_bool = 'Wheelchairs and accessories',
                     transportation_services_pre_12months_count_bool = 'Transportation services',
                     other_supplies_pre_12months_count_bool = 'Other supplies',
                     diabetic_footwear_pre_12months_count_bool = 'Diabetic footwear',
                     # W
                     fall_pre_12months_count_bool = '*Fall*',
                     other_injury_pre_12months_count_bool = '*Other injury*',
                     diverticular_disease_pre_12months_count_bool = '*Diverticular disease*',
                     hernia_pre_12months_count_bool = '*Hernia*',
                     arthropathy_pre_12months_count_bool = '*Arthropathy*',
                     GU_sx_pre_12months_count_bool = '*Genitourinary symptoms*',
                     optho2_pre_12months_count_bool = 'Ophthalmologic disease'
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
                   death = '**Overall mortality**', 
                   death.cause.specific = '**Cancer-specific mortality**', 
                   death.other.cause = '*Other-cause mortality*', 
                   death.other.cause = '*Overall other-cause mortality*', 
                   death.other.cause.gt90day = '*Overall other-cause mortality*', 
                   death.90.day = '**90-day mortality**', 
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


make.OR.plot  <-  function (odds.ratios_, label_list2, number.of.spaces = 50) {
    xlims <- c(0.5, 2)
    tt  <- sprintf("Risk Ratio\n  \U2190 favors SBRT %s favors Surgery \U2192", strrep(' ', number.of.spaces))
    # row.names(odds.ratios_) <- gsub("_pre_month_count", "", row.names(odds.ratios_))
     row.names(odds.ratios_) <- gsub("_..._count", "",  row.names(odds.ratios_))
    g <- ggplot(odds.ratios_, aes(x = estimate, y=y_axis)) + 
        geom_vline(aes(xintercept = 1), size = 0.25, linetype = "dashed") +
        geom_errorbarh(aes( xmax = high_ci, xmin = low_ci), size = 0.20, height = 0.3)+
        geom_point(size=1.5) +
        theme_bw() +
        theme(
              panel.grid.minor = element_blank(),
                  axis.text.y = ggtext::element_markdown(),
              ) +
        #scale_y_continuous(breaks = 1:max(odds.ratio$y_axis), labels = unlist(label_list[rownames(odds.ratios__propensity)]), trans='reverse') +
        scale_y_continuous(breaks = 1:max(odds.ratios_$y_axis), labels = label_list2[row.names(odds.ratios_)], trans='reverse') +
        #scale_x_continuous(breaks = seq(0,1.4,0.2), limits = xlims ) +
        scale_x_continuous(limits = xlims,breaks=c(0.5,0.75,1, 1.5,2) ) +
        coord_trans(x = "log10") +
        xlab(tt) +
        ylab("") +
        scale_linetype_manual(values=c("solid","dashed"))+
        scale_shape_manual(values=c(15,17))+
        theme( panel.grid.major.x = element_blank() ,
              panel.grid.major.y = element_line( size=.05, color="grey", linetype = 'dashed' ),
        legend.position = 'bottom',
      #  legend.title = element_blank(),
        legend.key.size=grid::unit(2,"lines")) 
 g
}

make.HD.plot  <-  function (odds.ratios_, label_list2, hazard =F, xlims = c(-0.05,0.10), breaks=c(-0.05,0.025,0, 0.05,0.075, 0.1)) {
    tt  <-  'Hazard Difference'
    tt  <- sprintf("Hazard Difference\n  \U2190 favors SBRT %s favors Surgery \U2192", strrep(' ', 30))
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
        theme( plot.title = element_text(size=12),
              axis.text=element_text(size=10),
              axis.title=element_text(size=10),
              panel.grid.major.x = element_blank() ,
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
# Given a lasso fit, there are several different ways to select the variables
# The default is lambda.min, which is used in the main analysis. Several other
# options are implemented here
################################
get.selected.columns.ahaz  <-  function(fit, s='lambda.min', cn=NULL, verbose=F, min.vars = 1, denom = 3) {
    lambda.min.idx  <-  which.min(fit$tunem)
    # Get the index of the smallest lambda which drops variables
    if (s == 'lambda.min') {
        coefs_  <-  fit$fit$beta[,lambda.min.idx]
    } else {
        last.idx  <-  max(which( diff(fit$fit$df) !=0))
        if (last.idx < lambda.min.idx) {
            if (verbose) print('All variables selected')
        }
        min.to.end.third  <- floor((last.idx - lambda.min.idx )/denom)
        if (s == 'lambda.min.lower' ) {
            lambda.idx  <- lambda.min.idx + min.to.end.third
            coefs_  <-  fit$fit$beta[,lambda.idx]
            if (verbose) print(sprintf('Lambda min is %f (%d), selecting lambda %f (%d)', fit$lambda[lambda.min.idx], fit$fit$df[lambda.min.idx], fit$lambda[lambda.idx], fit$fit$df[lambda.idx]))
        }else if (s == 'lambda.min.lowest' ) {
            lambda.idx  <- lambda.min.idx + 2*min.to.end.third
            coefs_  <-  fit$fit$beta[,lambda.idx]
            if (verbose) print(sprintf('Lambda min is %f (%d), selecting lambda %f (%d)', fit$lambda[lambda.min.idx], fit$fit$df[lambda.min.idx], fit$lambda[lambda.idx], fit$fit$df[lambda.idx]))
        }else if (s== 'lambda.1se') {
            candidates  <-  which(fit$tunem < min(fit$tunem)+ sd(fit$tunem)/ fit$foldsused[[1]]$nfolds)
            lambda.se1.idx  <-  candidates[which.min(fit$tunem[candidates])]
            lambda.se1  <-  fit$lambda[lambda.se1.idx]
            coefs_  <-  fit$fit$beta[,lambda.se1.idx]
            if (length(coefs_) < min.vars) {
                stop('Too few variables selected')
            }
        }else if (is.numeric(s) ) {
            lambda.idx  <- lambda.min.idx + min.to.end.third*s
            coefs_  <-  fit$fit$beta[,lambda.idx]
            if (verbose) print(sprintf('Lambda min is %f (%d), selecting lambda %f (%d)', fit$lambda[lambda.min.idx], fit$fit$df[lambda.min.idx], fit$lambda[lambda.idx], fit$fit$df[lambda.idx]))
        }
    }
    coefs  <-  coefs_ != 0
    if (!is.null(cn)) {
        names(coefs)  <-  cn
    }
    selected.columns  <- names(coefs)[c(-1)][coefs[c(-1)]]
    if (verbose ==T) print(sprintf('Excluded %s', paste(names(coefs)[c(-1)][!coefs[c(-1)]], collapse = ', ')))
    return(selected.columns)
}

################################
# Given a lasso fit, there are several different ways to select the variables
# The default is lambda.min, which is used in the main analysis. Several other
# options are implemented here
################################
get.selected.columns  <-  function(fit, s='lambda.min', cn=NULL, verbose=F, min.vars = 1, denom = 3) {
    lambda.min.idx  <-  which(fit$lambda == fit['lambda.min'])
    last.idx  <-  max(which( diff(fit$glmnet.fit$df) !=0))
    min.to.end.third  <- floor((last.idx - lambda.min.idx )/denom)
    if (s == 'lambda.min.lower' ) {
        lambda.idx  <- lambda.min.idx + min.to.end.third
        if (verbose) print(sprintf('Lambda min is %f (%d), selecting lambda %f (%d)', fit$lambda[lambda.min.idx], fit$fit$df[lambda.min.idx], fit$lambda[lambda.idx], fit$fit$df[lambda.idx]))
    }else if (s == 'lambda.min.lowest' ) {
        lambda.idx  <- lambda.min.idx + 2*min.to.end.third
        if (verbose) print(sprintf('Lambda min is %f (%d), selecting lambda %f (%d)', fit$lambda[lambda.min.idx], fit$glmnet.fit$df[lambda.min.idx], fit$lambda[lambda.idx], fit$glmnet.fit$df[lambda.idx]))
    }else if (s %in% c('lambda.1se')) {
        lambda.idx  <-  which(fit$lambda == fit[s])
        idx_2  <-  lambda.idx
        if (min.vars > fit$nzero[idx_2]) {
            if (verbose) print('Too few variables selected. Increasing lambda')
            idx_2  <-  idx_2 + 1
            while (fit$nzero[idx_2] < min.vars) idx_2  <-  idx_2 + 1
            if (verbose) print(sprintf('Decreased lambda to %f', fit$lambda[idx_2]))
        }
    }else if (s %in% c('lambda.min')) {
        lambda.idx  <-  which(fit$lambda == fit[s])
        idx_2  <-  lambda.idx
        if (min.vars > fit$nzero[idx_2]) {
            stop('Too few variables selected.' )
        }
    }else if (is.numeric(s) ) {
        lambda.idx  <- lambda.min.idx + min.to.end.third*s
        if (verbose) print(sprintf('Lambda min is %f (%d), selecting lambda %f (%d)', fit$lambda[lambda.min.idx], fit$fit$df[lambda.min.idx], fit$lambda[lambda.idx], fit$fit$df[lambda.idx]))
    }
    if (!is.null(cn)) {
        names(coefs)  <-  cn
    }
    coefs  <-  coef(fit, s = fit$lambda[lambda.idx])[,1] != 0
    selected.columns  <- names(coefs)[c(-1,-2)][coefs[c(-1,-2)]]
    if (verbose ==T) print(sprintf('Excluded %s', paste(names(coefs)[c(-1,-2)][!coefs[c(-1,-2)]], collapse = ', ')))
    return(selected.columns)
}
quartile  <- function(x) {
    scaled  <-  as.numeric(x)
    breaks  <- c(0,quantile(scaled[scaled!=0], probs = c(0, 0.25, 0.5, 0.75), na.rm = T, names=F), max(scaled, na.rm = T))
    scaled  <-  cut(x, breaks, include.lowest = T, labels = c('0', '1', '2', '3', '4'), right=F) %>% as.character %>% as.numeric
    return(scaled)
}
check.if.present  <- function(needles,haystack) sapply(needles, function(x) any(grepl(sprintf('^%s',x), haystack)))

################################
# Backup 
################################
# system.time ( MBSF_2013 <- read.SAScii('../SEER-Medicare-data/data/SEER_Medicare/mbsf.abcd.summary.2013.txt', 
#             '../SEER-Medicare-data/data/SEER_Medicare/2020 Input Statements/MBSF.ABCD.gcl.txt', n=100000) ) 

# lung.SEER <- read.SAScii('../SEER-Medicare-data/data/SEER_Medicare/SEER.lung.cancer.n1000.txt',  
                         # '../SEER-Medicare-data/data/SEER_Medicare/2020 Input Statements/SEER.input.gcl.txt', n=100) 

