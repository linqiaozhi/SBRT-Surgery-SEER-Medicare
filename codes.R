library(icd)
library(tidyverse)
CPT_Codes  <-  read_csv('CPT_Codes_2020.csv') %>% rename(HCPC=`HCPC/MOD`, short_desc = `SHORT DESCRIPTION`, long_desc=`LONG DESCRIPTION`)

expand_range_procs   <-  function( from, to, CPT_Codes) {
    # get first character
    if (substr(from, 1, 1) == substr(to, 1, 1)) {
        nums_between = as.numeric(gsub('[^0-9]', '', from)):as.numeric(gsub('[^0-9]', '', to))
        # add first character back
        nums_between = sprintf('%s%04d', substr(from, 1, 1), nums_between)
        codes  <- CPT_Codes %>% filter ( HCPC %in% nums_between) %>% pull(HCPC)
    }else {
        error('First character of from and to must be the same')
    }
    return (codes)
}

# apply expand_range to each
expand.each.code  <-  function(code.list, icd9or10 = 'undefined') {
    out.list  <-  c()
    for (i in 1:length(code.list)) {
        code  <- code.list[i]
        if (icd9or10 == 'undefined') {
            try.out <- try(expand_range(code, code), silent = T)
            if (!inherits(try.out, "try-error")) {
                out.list  <-  c(out.list, try.out)
            }else{
                out.list  <-  c( out.list, code)
                warning(sprintf('Could not find %s', code))
            }
        } else if (icd9or10 == 'icd10') {
            try.out <- try(expand_range(as.icd10cm(code), as.icd10cm(code)), silent = T)
            if (!inherits(try.out, "try-error")) {
                out.list  <-  c(out.list, try.out)
            }else{
                out.list  <-  c( out.list, code)
                warning(sprintf('Could not find %s', code))
            }
        }else{
            out.list  <-  c(out.list, expand_range(as.icd9(code), as.icd9(code)))
        }
    }
    return(out.list)
}


valid.dxs  <- c( expand_range('1622','1629'), expand_range(as.icd10('C34'), as.icd10('C349')))
sbrt.icds  <-  c('9230', '9231', '9232', '9233', '9239',
                    'DB22DZ', 'DB22HZZ', 'DB22JZZ')
sublobar.icds  <-  c(  '3230', '3239', '3220', '3229', '0BBC4ZX', '0BBC4ZZ', '0BBC0ZX', '0BBC0ZZ', '0BBD4ZX', '0BBD4ZZ', '0BBD0ZX', '0BBD0ZZ', '0BBF4ZX', '0BBF4ZZ', '0BBF0ZX', '0BBF0ZZ', '0BBG4ZX', '0BBG4ZZ', '0BBG0ZX', '0BBG0ZZ', '0BBH4ZX', '0BBH4ZZ', '0BBH0ZX', '0BBH0ZZ', '0BBJ4ZX', '0BBJ4ZZ', '0BBJ0ZX', '0BBJ0ZZ', '0BBK4ZX', '0BBK4ZZ', '0BBK0ZX', '0BBK0ZZ', '0BBL4ZX', '0BBL4ZZ', '0BBL0ZX', '0BBL0ZZ', '0BBM4ZX', '0BBM4ZZ', '0BBM0ZX', '0BBM0ZZ')
other.resection.icds  <-  c( '3240', '3241', '3249', '3260', '3250', '3259', '0BTC4ZX', '0BTC4ZZ', '0BTC0ZX', '0BTC0ZZ', '0BTD4ZX', '0BTD4ZZ', '0BTD0ZX', '0BTD0ZZ', '0BTF4ZX', '0BTF4ZZ', '0BTF0ZX', '0BTF0ZZ', '0BTG4ZX', '0BTG4ZZ', '0BTG0ZX', '0BTG0ZZ', '0BTH4ZX', '0BTH4ZZ', '0BTH0ZX', '0BTH0ZZ', '0BTJ4ZX', '0BTJ4ZZ', '0BTJ0ZX', '0BTJ0ZZ', '0BTK4ZX', '0BTK4ZZ', '0BTK0ZX', '0BTK0ZZ', '0BTL4ZX', '0BTL4ZZ', '0BTL0ZX', '0BTL0ZZ', '0BTM4ZX', '0BTM4ZZ', '0BTM0ZX', '0BTM0ZZ')

#https://cdn.jamanetwork.com/ama/content_public/journal/surg/931810/soi140036supp1_prod.pdf?Expires=1684345221&Signature=OSPftfHA9rDQPKLJrRnmTuVgeN80y9G1jEx1SFmWPu8xdjKv-GfJLKrub~yKR12CYlUhXH6JMBap3PqTII4naX9ZndV10ahKSFYdWMciw4diKk9daxaK0V~4m7pTmaMw~a7H6H-Xqcke5-g~XeFlSz18EX08a-d5MQsC5Os0Um9v8Ng~y7-xB8t2tJygyztL2UqJqhdejzcXdUWCERoQ~fWXJJh640Qer8NojGa8Va2YT8mEXhu3q3GEsVUJ~dEooFq5X17r5w4~qzSfIylot3LJVVGxo9eslbLjfnumtjUIlbdNlLx8zoguvAjP8VoHUEpSQI88n~wOwPpKht68Bg <- &Key-Pair-Id=APKAIE5G5CRDK6RD3PGA


sbrt.cpts  <-  c('77373', 'G0173', 'G0251', 'G0339', 'G0340', '61793',  '0082T' )


dx.icd   <- list (
     'smoking'       = list (
        'icd9' = expand.each.code(c( 'V1582', '3051')),
        'icd10' = expand.each.code(c('Z87891'))
     ),
     'o2'       = list (
        'icd9' = expand.each.code(c('V462')),
        'icd10' = expand.each.code(c('Z9981'))
     ),
     # 'cerebrovascular_disease' = list( 
     #    'icd9' = expand_range( '430', '438'),
     #    'icd10' = expand.each.code(c( 'I60', 'I61', 'I62', 'I65', 'I66', 'I67', 'I69', 'G45')) 
     # ),
     # 'copd_and_allied_conditions' = list(
     #    'icd9' = expand_range( '490', '496'),
     #    'icd10' = expand.each.code(c( expand_range('J40', 'J45'), 'J47', 'J67')) 
     # ),
     'other_bacterial_diseases' = list(
        'icd9' = c(expand_range( '030', '041')),
        'icd10' =c(expand.each.code(c( 'A30', 'A31', 'A46', 'A48', 'A49')), expand_range('A35', 'A42')) 
     ),
     'pneumonia_and_influenza' = list( 
        'icd9' = c(expand_range( '480', '487')),
        'icd10' = c( expand_range('J09', 'J13'), expand_range('J15', 'J18')) 
     ),
     'non_diabetes_endocrine' = list(
        'icd9' = setdiff( expand_range( '240', '259'), expand_range( '249', '250'     )),
        'icd10' = expand_range( as.icd10cm('E00'), as.icd10cm('E07')     ), expand_range( as.icd10cm('E15'), as.icd10cm('E35')     )
     ),
     'pressure_ulcer' = list(
        'icd9' = expand_range( '707', '70725'),
        'icd10' = expand_range('L89', 'L89')
     ),
     # 'kidney_diseases' = list(
     #    'icd9' = expand_range( '580', '589'),
     #    'icd10' = expand.each.code(c( 'N00', 'N03', 'N04', 'N05', 'N07', 'N17', 'N18', 'N19', 'N25', 'N269', 'N27')     )
     # ),
     'fall' = list( 
        'icd9' = expand_range('E880','E888'),
        'icd10' = expand_range('W00', 'W19'     )),
     'other_injury' = list( 
        'icd9' = expand_range('800', '999'     ),
        'icd10' = expand_range('S00','T79')),
     'acute_bronchitis' = list(
        'icd9'=c('4660', '4661'),
        'icd10' =  sprintf('J20%d',0:9)),
     'oral' = list(
        'icd9' = expand_range('520','5299'),
        'icd10' =  expand_range('K00', 'K149')),
      'pancreatic' = list(
         'icd9' = expand_range('577','577'),
         'icd10' =  expand_range('K85', 'K86')),
     'gout' = list(
        'icd9' = expand_range('2740','2749'),
        'icd10' =  expand_range('M100', 'M109')),
     'arthropathy' = list(
        'icd9' = setdiff(expand_range('711', '715'), expand_range( '2740', '2749')),
        'icd10' = setdiff ( expand_range('M00', 'M19'), expand_range('M100', 'M109'))),
     'GU_sx' = list(
        'icd9' = c(expand_range('590', '599'), expand_range('788','78899')),
        'icd10' = c(  expand_range('R30', 'R39')), expand_range('N30', 'N39')),
     'diverticular_disease' = list(
        'icd9' = expand_range('562','56213'     ) ,
        'icd10' =  expand_range('K57', 'K5793')     ),
     'hernia' = list(
        'icd9' = c( expand_range('550', '5539')     ) ,
        'icd10' =  expand_range('K40', 'K469')     ),
     'hemorrhoids' = list(
        'icd9' = c( expand_range('4550', '4559')     ) ,
        'icd10' =  expand_range('K640', 'K649')     ),
     'optho' = list(
        'icd9' = c( expand_range('360', '379')     ) ,
        'icd10' =  expand_range('H00', 'H59')     ),
     # 'optho2' = list(
     #    'icd9' = c( 
     #            expand_range('362', '36218'),  # Diabetic, hypertensive, and other retinopathy
     #            expand_range('363', '36335') , # Uveitis
     #            expand_range('364', '3643')  ,
     #            expand_range('3623', '36237'),  # Retinal vascular occlusion                                     
     #            expand_range('37034', '37034'),  # exposure keratitis                                      
     #            expand_range('37741', '37741')  # ischemic optic neuropathy
     # ),
     #    'icd10' =  c( 
     #            expand_range(as.icd10cm('E083'),as.icd10cm('E0839')), # Diabetic retinop
     #            expand_range(as.icd10cm('E093'),as.icd10cm('E0939')),
     #            expand_range(as.icd10cm('E103'),as.icd10cm('E1039')),
     #            expand_range(as.icd10cm('E113'),as.icd10cm('E1139')),
     #            expand_range(as.icd10cm('H35'),as.icd10cm('H3509')), # Other retinal disorders, including hypertensive retinopathy
     #            expand_range(as.icd10cm('H20'),as.icd10cm('H209')), # Uveitis
     #            expand_range(as.icd10cm('H30'),as.icd10cm('H309')), 
     #            expand_range(as.icd10cm('H4411'),as.icd10cm('H44119')),
     #            expand_range(as.icd10cm('H34'),as.icd10cm('H349')), # retinal vascular occlusiosn
     #            expand_range(as.icd10cm('H4701'),as.icd10cm('H47019')) # exposure keratopathy
     # )
     # ),
     'ischemic_heart_disease' = list(
        'icd9' = expand_range( '410', '414'),
        'icd10' = expand_range( 'I20', 'I25') 
     ),
     # 'MI' = list(
     #    'icd9' = c('410','412') ,
     #    'icd10' =  c('I21','I22','I252')),
     'CHF' = list(
        # 'icd9' = c('39891','40201','40211','40291','40401','40403','40411','40413','40491','40493','4254','4255','4257','4258','4259','428') ,
         'icd9' = expand.each.code(c('39891','40201','40211','40291','40401','40411','40491','4254','4255','4257','4258','4259','428')) ,
        # 'icd10' =  c('I43','I50','I099','I110','I130','I132', 'I255','I420','I425','I426','I427','I428','I429','P290')),
         'icd10' =  expand.each.code(c('I43','I50','I099','I110','I130','I132', 'I420','I425','I426','I427','I428','I429','P290'))),
     'PVD' = list(
        'icd9' = expand.each.code(c('0930','440','441','4431','4432','4438','4439','4471','5571','5579','V434')),
        'icd10' =  expand.each.code(c('I70','I71', 'I731','I738','I739','I771', 'I790','K551','K558','K559','Z958','Z959'))), # removed 'I792',
     'CVD' = list(
        # 'icd9' = c('36234','430','431','432','433','434','435','436','437','438'),
         'icd9' = expand.each.code(c('430','431','432','433','434','435','436','437','438')),
        # 'icd10' =  c('G45','G46','I60','I61','I62','I63','I64','I65','I66','I67','I68','I69','H340')),
        'icd10' =  c( 'I64', expand.each.code(c('G45','G46','I60','I61','I62','I63','I65','I66','I67','I68','I69')))),
     'dementia' = list(
        'icd9' = expand.each.code(c('290','2941','3312')),
        'icd10' =  expand.each.code(c('F01','F02','F03','G30','G311'))),# removed F00*
     'COPD' = list(
        'icd9' = expand.each.code(c('4168','4169','490','491','492','493','494','495','496','500','501','502','503','504','505','5064','5081','5088')),
        'icd10' =  expand.each.code(c('J40','J41','J42','J43','J44','J45','J46','J47', 'J60','J61','J62','J63','J64','J65','J66','J67', 'I278','I279','J684','J701','J703'))),
     'PUD' = list(
        'icd9' = expand.each.code(c('531','532','533','534')),
        'icd10' =  expand.each.code(c('K25','K26','K27','K28'))),
     'MILDLD' = list(
        'icd9' = expand.each.code(c('07022','07023','07032','07033','07044','07054','0706','0709','570','571','5733','5734','5738','5739','V427')),
        'icd10' =  expand.each.code(c('B18','K73','K74','K700','K701','K702','K703','K709', 'K717','K713','K714','K715','K760','K762','K763','K764','K768','K769','Z944'))),
     'DIAB_UC' = list(
        'icd9' = expand.each.code(c('2500','2501','2502','2503','2508','2509')),
        'icd10' =  expand.each.code(c('E100','E101','E106','E108','E109','E110','E111','E116','E118','E119', 'E120','E121','E126','E128','E129', 'E130','E131','E136','E138','E139', 'E140','E141','E146','E148','E149'), icd9or10 = 'icd10')),
     'DIAB_C' = list(
        'icd9' = expand.each.code(c('2504','2505','2506','2507')),
        'icd10' =  expand.each.code(c('E102','E103','E104','E105','E107', 'E112','E113','E114','E115','E117', 'E122','E123','E124','E125','E127', 'E132','E133','E134','E135','E137', 'E142','E143','E144','E145','E147'), 'icd10')),
     'PARA' = list(
        'icd9' = expand.each.code(c('3341','342','343','3440','3441','3442','3443','3444','3445','3446','3449')),
        'icd10' =  expand.each.code(c('G81','G82','G041','G114','G801','G802', 'G830','G831','G832','G833','G834','G839'))),
     'RD' = list(
        'icd9' = expand.each.code(c('40301','40311','40391','40402','40403','40412','40413','40492','40493','582','5830','5831','5832','5834','5836','5837','585','586','5880','V420','V451','V56')),
        'icd10' =  expand.each.code(c('N18','N19','N052','N053','N054','N055','N056','N057', 'N250','I120','I131','N032','N033','N034','N035','N036','N037', 'Z490','Z491','Z492','Z940','Z992'))),
     'cancer_nonlung' = list(
        'icd9' = setdiff(expand.each.code(c('140','141','142','143','144','145','146','147','148','149','150','151','152','153','154','155','156','157','158','159','160','161','162','163','164','165','170','171','172','174','175','176','179','180','181','182','183','184','185','186','187','188','189','190','191','192','193','194','195','200','201','202','203','204','205','206','207','208','2386')),valid.dxs),
        'icd10' =  setdiff(expand.each.code(c('C00','C01','C02','C03','C04','C05','C06','C07','C08','C09', 'C10','C11','C12','C13','C14','C15','C16','C17','C18','C19', 'C20','C21','C22','C23','C24','C25','C26', 'C30','C31','C32','C33','C34','C37','C38','C39', 'C40','C41','C43','C45','C46','C47','C48','C49', 'C50','C51','C52','C53','C54','C55','C56','C57','C58', 'C60','C61','C62','C63','C64','C65','C66','C67','C68','C69', 'C70','C71','C72','C73','C74','C75','C76', 'C81','C82','C83','C84','C85','C88', 'C90','C91','C92','C93','C94','C95','C96','C97')),valid.dxs)),
     'MSLD' = list(
        'icd9' = expand.each.code(c('4560','4561','4562','5722','5723','5724','5728')),
        'icd10' =  expand.each.code(c('K704','K711','K721','K729','K765','K766','K767','I850','I859','I864','I982'))),
     'METS' = list(
        'icd9' = expand.each.code(c('196','197','198','199')),
        'icd10' =  expand.each.code(c('C77','C78','C79','C80'))),
     'HIV' = list(
        'icd9' =expand.each.code( c('042','043','044')),
        'icd10' = expand.each.code( c('B20','B21','B22','B24')))
)
dx.icd  <- c( dx.icd, 
     list ( 'mental_disorders' = list(
        'icd9' = setdiff( expand_range('290', '294'), dx.icd[['dementia']]$icd9) ,
        'icd10' = setdiff(expand_range('F01','F99'), dx.icd[['dementia']]$icd10)
     ),
     'nervous_system' = list(
        'icd9' = setdiff( expand_range( '320', '359'), c(dx.icd[['PARA']]$icd9,  dx.icd[['CVD']]$icd9,  dx.icd[['dementia']]$icd9  )),
        'icd10' =setdiff( expand_range('G00', 'G99'),c(  dx.icd[['PARA']]$icd10, dx.icd[['CVD']]$icd10, dx.icd[['dementia']]$icd10 )) ),
     'other_heart_disease' = list(
        'icd9' = setdiff ( expand_range( '390', '429'), 
                          c( dx.icd[['RD']]$icd9,dx.icd[['CHF']]$icd9,dx.icd[['ischemic_heart_disease']]$icd9 , dx.icd[['COPD']]$icd9)  ),
        'icd10' =setdiff ( expand_range('I00','I52'),  
                          c( dx.icd[['RD']]$icd10,dx.icd[['CHF']]$icd10 , dx.icd[['ischemic_heart_disease']]$icd10, dx.icd[['COPD']]$icd10  ))
     ),
     'veins_lymphatics_other_circulatory' = list(
        'icd9' = setdiff(  c(expand_range( '451', '459')), c (  dx.icd[['MSLD']]$icd9, dx.icd[['hemorrhoids']]$icd9 )),
        'icd10' =setdiff(  c( expand.each.code(c( 'I80',  'I86', 'I89', 'I95', 'I99', 'K64')), expand_range('I81', 'I83') ),  c (  dx.icd[['MSLD']]$icd10,dx.icd[['hemorrhoids']]$icd10) )
     ),
     'rheum' = list(
        'icd9' = setdiff( expand.each.code(c('4465','7100','7101','7102','7103','7104','7140','7141','7142','7148','725')),
                         dx.icd[['arthropathy']]$icd9),
        'icd10' =  setdiff(
                           expand.each.code(c('M05','M32','M33','M34','M06','M315','M351','M353','M360')),
                          dx.icd[['arthropathy']]$icd10 ))
     ))


# check for duplicates
dx.icd.codes  <-  unlist(dx.icd)
dx.icd.codes[dx.icd.codes %in%  dx.icd.codes[ duplicated(dx.icd.codes)]]


pet.scan.cpts <-c('78811', '78812', '78813', '78814', '78815', '78816', 'G0235')


procs  <- list (
            'hospital_beds_and_supplies' = expand_range_procs ( 'E0250', 'E0373', CPT_Codes),
            'wheelchairs_accessories' = c ( expand_range_procs ( 'K0001', 'K0462', CPT_Codes), 'K0669'),
            'walking_aids' = c ( expand_range_procs ( 'E0100', 'E0159', CPT_Codes)),
            'O2accessories' = c ( expand_range_procs ( 'E1353', 'E1406', CPT_Codes)),
            'other_supplies' =  expand_range_procs ( 'A4244', 'A4290', CPT_Codes),
            'diabetic_footwear' =  expand_range_procs ( 'A5500', 'A5513', CPT_Codes),
            'transportation_services' =  expand_range_procs ( 'A0021', 'A0999', CPT_Codes)
                )




sink('tbls/dx.icd.txt'); 
for (namei in (names(dx.icd))) { 
    cat(sprintf('Variable Name: %s', namei))
    print('ICD9')
    icd9.codes  <-  dx.icd[[namei]][['icd9']]
    cat(sprintf('%s\n', paste(icd9.codes, explain_code(condense=F, as.icd9(icd9.codes)))))
    print('ICD10')
    icd10.codes  <-  dx.icd[[namei]][['icd10']]
    cat(sprintf('%s\n', paste(icd10.codes, explain_code(condense=F, as.icd10cm(icd10.codes)))))
    cat('------------------------------\n\n\n\n')
} 
sink()
#TODO:Explain_code doesn't work for the ICD10 codes starting with E. The expand_range does, so the code all works, but the output above gives NA for the E* codes. need to use explain_icd9_10( as.icd10('E033'))

sink('tbls/proc.codes.txt'); 
for (namei in names(procs)) { 
    cat(sprintf('Variable Name: %s', namei))
    CPT_Codes %>% filter (HCPC %in% procs[[namei]] ) %>% select(HCPC, long_desc) %>% print(n=Inf)
    cat('------------------------------\n\n\n\n')
} 
sink()
