valid.dxs  <- c( expand_range('1622','1629'), expand_range(as.icd10('C34'), as.icd10('C349')))
sbrt.icds  <-  c('9230', '9231', '9232', '9233', '9239',
                    'DB22DZ', 'DB22HZZ', 'DB22JZZ')
sublobar.icds  <-  c(  '3230', '3239', '3220', '3229', '0BBC4ZX', '0BBC4ZZ', '0BBC0ZX', '0BBC0ZZ', '0BBD4ZX', '0BBD4ZZ', '0BBD0ZX', '0BBD0ZZ', '0BBF4ZX', '0BBF4ZZ', '0BBF0ZX', '0BBF0ZZ', '0BBG4ZX', '0BBG4ZZ', '0BBG0ZX', '0BBG0ZZ', '0BBH4ZX', '0BBH4ZZ', '0BBH0ZX', '0BBH0ZZ', '0BBJ4ZX', '0BBJ4ZZ', '0BBJ0ZX', '0BBJ0ZZ', '0BBK4ZX', '0BBK4ZZ', '0BBK0ZX', '0BBK0ZZ', '0BBL4ZX', '0BBL4ZZ', '0BBL0ZX', '0BBL0ZZ', '0BBM4ZX', '0BBM4ZZ', '0BBM0ZX', '0BBM0ZZ')
other.resection.icds  <-  c( '3240', '3241', '3249', '3260', '3250', '3259', '0BTC4ZX', '0BTC4ZZ', '0BTC0ZX', '0BTC0ZZ', '0BTD4ZX', '0BTD4ZZ', '0BTD0ZX', '0BTD0ZZ', '0BTF4ZX', '0BTF4ZZ', '0BTF0ZX', '0BTF0ZZ', '0BTG4ZX', '0BTG4ZZ', '0BTG0ZX', '0BTG0ZZ', '0BTH4ZX', '0BTH4ZZ', '0BTH0ZX', '0BTH0ZZ', '0BTJ4ZX', '0BTJ4ZZ', '0BTJ0ZX', '0BTJ0ZZ', '0BTK4ZX', '0BTK4ZZ', '0BTK0ZX', '0BTK0ZZ', '0BTL4ZX', '0BTL4ZZ', '0BTL0ZX', '0BTL0ZZ', '0BTM4ZX', '0BTM4ZZ', '0BTM0ZX', '0BTM0ZZ')

#https://cdn.jamanetwork.com/ama/content_public/journal/surg/931810/soi140036supp1_prod.pdf?Expires=1684345221&Signature=OSPftfHA9rDQPKLJrRnmTuVgeN80y9G1jEx1SFmWPu8xdjKv-GfJLKrub~yKR12CYlUhXH6JMBap3PqTII4naX9ZndV10ahKSFYdWMciw4diKk9daxaK0V~4m7pTmaMw~a7H6H-Xqcke5-g~XeFlSz18EX08a-d5MQsC5Os0Um9v8Ng~y7-xB8t2tJygyztL2UqJqhdejzcXdUWCERoQ~fWXJJh640Qer8NojGa8Va2YT8mEXhu3q3GEsVUJ~dEooFq5X17r5w4~qzSfIylot3LJVVGxo9eslbLjfnumtjUIlbdNlLx8zoguvAjP8VoHUEpSQI88n~wOwPpKht68Bg <- &Key-Pair-Id=APKAIE5G5CRDK6RD3PGA


sbrt.cpts  <-  c('77373', 'G0173', 'G0251', 'G0339', 'G0340', '61793',  '0082T' )
manual.comorbidities   <- list (
                                'smoking'       = list (
                                                        'icd9' = c( 'V1582', '3051'),
                                                        'icd10' = c('Z87891')
                                                        ),
                                'o2'       = list (
                                                        'icd9' = c('V462'),
                                                        'icd10' = c('Z9981')
                                                    )
                                )

pet.scan.cpts <-c('78811', '78812', '78813', '78814', '78815', '78816', 'G0235')


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
#    'cholelithiasis' = list(
#                            'icd9' = expand_range('5740', '57491'),
#                            'icd10' =  expand_range('K80', 'K8081')),
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
                            'icd10' =  expand_range('H00', 'H59') ),
    'optho2' = list(
                            'icd9' = c( 
                                       expand_range('362', '36218'),  # Diabetic, hypertensive, and other retinopathy
                                       expand_range('363', '36335') , # Uveitis
                                       expand_range('364', '3643')  ,
                                       expand_range('3623', '36237'),  # Retinal vascular occlusion                                     
                                       expand_range('37034', '37034'),  # exposure keratitis                                      
                                       expand_range('37741', '37741')  # ischemic optic neuropathy
                                       ),
                            'icd10' =  
                               c( 
                                 expand_range(as.icd10cm('E083'),as.icd10cm('E0839')), # Diabetic retinop
                                 expand_range(as.icd10cm('E093'),as.icd10cm('E0939')),
                                 expand_range(as.icd10cm('E103'),as.icd10cm('E1039')),
                                 expand_range(as.icd10cm('E113'),as.icd10cm('E1139')),
                                 expand_range(as.icd10cm('H35'),as.icd10cm('H3509')), # Other retinal disorders, including hypertensive retinopathy
                                 expand_range(as.icd10cm('H20'),as.icd10cm('H209')), # Uveitis
                                 expand_range(as.icd10cm('H30'),as.icd10cm('H309')), 
                                 expand_range(as.icd10cm('H4411'),as.icd10cm('H44119')),
                                 expand_range(as.icd10cm('H34'),as.icd10cm('H349')), # retinal vascular occlusiosn
                                 expand_range(as.icd10cm('H4701'),as.icd10cm('H47019')) # exposure keratopathy
                                 )
    )
    )


if (F) { 
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
}
