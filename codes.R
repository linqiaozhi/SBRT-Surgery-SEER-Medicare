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
  ),
  'organic psychotic conditions' = list(
    'icd9' = expand_range(
      '290', '295'
    ),
    'icd10' = c(
      'F01', 'F02', 'F03', 'F05', 'F06', expand_range('F1015', 'F1019'), expand_range('F1023', 'F1029'), expand_range('F1192', 'F1198'), expand_range('F1292', 'F1298'),  expand_range('F1392', 'F1398'),  expand_range('F1492', 'F1498'), expand_range('F1592', 'F1598'),  expand_range('F1692', 'F1698'), 'F17',  expand_range('F1892', 'F1898'),  expand_range('F1992', 'F1998')
    ) ), 
  'hereditary and degenerative diseases of the central nervous system' = list(
    'icd9' = expand_range(
      '330', '338'
    ),
    'icd10' = c(
      'G11', 'G12', 'G20', 'G23', 'G30', 'G31', 'G32', 'G89', 'G90', 'G95'
    )
  ),
  'other psychoses' = list(
    'icd9' = expand_range(
      '295', '299'
    ),
    'icd10' = c(
      'F20', 'F22', 'F28',  expand_range('F30', 'F39'), 'F84'
    ) 
  ),
  'other forms of heart disease' = list(
    'icd9' = expand_range(
      '410', '429'
    ),
    'icd10' = c(
      'I30', 'I31', 'I33', 'I34', 'I40', 'I42', 'I45', 'I49', 'I50', 'I51'
    ) 
  ),
  'open wound of lower limb' = list(
    'icd9' = expand_range(
      '890', '897'
    ),
    'icd10' = c(
      'S71', 'S81', 'S88', expand_range('S910', 'S913'), 'S98'
    ) 
  ),
  'ischemic heart disease' = list(
    'icd9' = expand_range(
      '410', '414'
    ),
    'icd10' = c(
      'I20', 'I21', 'I24', 'I252', 'I258'
    ) 
  ),
  'hypertensive disease' = list(
    'icd9' = expand_range(
      '401', '405'
    ),
    'icd10' = expand_range(
      'I10', 'I15'
    ) 
  ),
  'cerebrovascular disease' = list(
    'icd9' = expand_range(
      '430', '438'
    ),
    'icd10' = c(
      'I60', 'I61', 'I62', 'I65', 'I66', 'I67', 'I69', 'G45'
    ) 
  ),
  'neurotic disorders, personality disorders, and other nonpsychotic mental disorders' = list(
    'icd9' = expand_range(
      '300', '316'
    ),
    'icd10' = c(
      'F04', 'F07', 'F09', 'F1020', 'F1021', 'F111', 'F112', 'F121', 'F122', 'F131', 'F132', 'F141', 'F142', 'F151', 'F152', 'F161', 'F162', 'F181', 'F182', 'F191', 'F192', 'F32', expand_range('F40', 'F48'), 'F50', 'F51', 'F54', expand_range('F60', 'F69'), 'F80', 'F81', 'F91', 'F93', 'F95', 'F98', 'F99'
    )

  ),
  'arthropathies and related disorders' = list(
    'icd9' = expand_range(
      '710', '719'
    ),
    'icd10' = c(
      'M01', 'M06', expand_range('M1180', 'M1189'), 'M12', 'M148', 'M15', expand_range('M23', 'M25'), 'M35'
    ) 
  ),
  'chronic obstructive pulmonary disease and allied conditions' = list(
    'icd9' = expand_range(
      '490', '496'
    ),
    'icd10' = c(
      expand_range('J40', 'J45'), 'J47', 'J67'
    ) 
  ),
  'other bacterial diseases' = list(
    'icd9' = expand_range(
      '030', '041'
    ),
    'icd10' = c(
      'A30', 'A31', expand_range('A35', 'A42'), 'A46', 'A48', 'A49'
    ) 
  ),
  'diseases of veins and lymphatics, and other diseases of circulatory system' = list(
    'icd9' = expand_range(
      '451', '459'
    ),
    'icd10' = c(
      'I80', expand_range('I81', 'I83'), 'I86', 'I89', 'I95', 'I99', 'K64'
    ) 
  ),
  'pneumonia and influenza' = list(
    'icd9' = expand_range(
      '480', '487'
    ),
    'icd10' = c(
      expand_range('J09', 'J13'), expand_range('J15', 'J18')
    ) 
  ),
  'diseases of other endocrine glands' = list(
    'icd9' = expand_range(
      '250', '259'
    ),
    'icd10' = c(
      'E08', 'E09', 'E10', 'E11', 'E12', 'E13', 'E16', 'E21', 'E22', 'E23', 'E27', 'E28', 'E29', 'E31', 'E32', 'E34'
    )
  ),
  'other diseases of urinary system' = list(
    'icd9' = expand_range(
      '590', '599'
    ),
    'icd10' = c(
      expand_range('N10', 'N16'), 'N133', 'N20', 'N21', expand_range('N25', 'N29'), 'N30', 'N32', 'N34', 'N35', 'N39'
    ) 
  ),
  'ill-defined and unknown causes of morbidity and mortality' = list(
    'icd9' = expand_range(
      '797', '799'
    ),
    'icd10' = c(
      'R4181', 'R99'
    ) 
  ),
  'contusion with intact skin surface' = list(
    'icd9' = expand_range(
      '920', '924'
    ),
    'icd10' = expand_range(
      'S00', 'S99' 
    ) 
  ),
  'nephritis, nephrotic syndrome, and nephrosis' = list(
    'icd9' = expand_range(
      '580', '589'
    ),
    'icd10' = c(
      'N00', 'N03', 'N04', 'N05', 'N07', 'N17', 'N18', 'N19', 'N25', 'N269', 'N27'
    ) 
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
sink('tbls/manual.comorbidites.txt'); 
for (namei in (names(manual.comorbidities))) { 
    cat(sprintf('Variable Name: %s', namei))
    print('ICD9')
    icd9.codes  <-  manual.comorbidities[[namei]][['icd9']]
    cat(sprintf('%s\n', paste(icd9.codes, explain_code(condense=F, as.icd9(icd9.codes)))))
    print('ICD10')
    icd10.codes  <-  manual.comorbidities[[namei]][['icd10']]
    cat(sprintf('%s\n', paste(icd10.codes, explain_code(condense=F, as.icd10(icd10.codes)))))
    cat('------------------------------\n\n\n\n')
} 
sink()
}
