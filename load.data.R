load.data   <-  function (data.num, subset.name,  nc_time_days = 90, W1s.global = F ) {
    if (W1s.global ){
        W1s  <-  list(
                      'death' = 'death.other.cause', 
                      'death.lc.specific' = 'death.other.cause', 
                      'death.other.cause' = 'death.other.cause',
                      'death.copd' = 'death.noncopd', 
                      'death.heart' = 'death.nonheart',
                      'death.stroke' = 'death.nonstroke',
                      'death.noncopd.nonheart.nonstroke' = 'death.copd.heart.stroke' ,
                      'death.other.1' = 'death.other.non.1', 
                      'death.other.2' = 'death.other.non.2',
                      'death.other.3' = 'death.other.non.3',
                      'death.other.4' = 'death.other.non.4'
        )
    }else {
        W1s  <-  list(
                      'death' = 'death.other.cause', 
                      'death.lc.specific' = 'death.other.cause', 
                      'death.other.cause' = 'death.other.cause', 
                      'death.copd' = 'death.noncopd', 
                      'death.heart' = 'death.nonheart',
                      'death.stroke' = 'death.nonstroke',
                      'death.noncopd.nonheart.nonstroke' = 'death.copd.heart.stroke' 
        )
    }
    outcome.names  <- names(W1s)

    ################################
    # Load data 
    ################################
    filename.in  <-  sprintf('data/A.final%d.%s.RDS',data.num, subset.name)
    A.final  <-  readRDS(filename.in)  %>% 
        mutate(treatment.year = year(tx.date),
               death.90.day = if_else ( ninety.day.mortality, death, as.Date(NA_character_)),
               death.lc.specific = if_else ( lc.specific.mortality == 'Death', death, as.Date(NA_character_)),
               death.other.cause = if_else ( other.cause.mortality == 'Death' , death, as.Date(NA_character_)),
               death.copd = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE == '50130',  death, as.Date(NA_character_)),
               death.noncopd = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50130' , death, as.Date(NA_character_)),
               death.heart = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE == '50060', death, as.Date(NA_character_)),
               death.nonheart = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50060' , death, as.Date(NA_character_)),
               death.other = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE == '50300', death, as.Date(NA_character_)),
               death.nonother = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50300', death, as.Date(NA_character_)),
               death.stroke = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE == '50080', death, as.Date(NA_character_)),
               death.nonstroke = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50080', death, as.Date(NA_character_)),
               death.noncopd.nonheart.nonstroke = if_else ( other.cause.mortality == 'Death' & COD_TO_SITE_RECODE != '50130' & COD_TO_SITE_RECODE != '50060' & COD_TO_SITE_RECODE != '50080' , death, as.Date(NA_character_)),
               death.copd.heart.stroke = if_else ( other.cause.mortality == 'Death' & ( COD_TO_SITE_RECODE == '50130' | COD_TO_SITE_RECODE == '50060' | COD_TO_SITE_RECODE == '50080' ) , death, as.Date(NA_character_)),
        )
    A.final  <- A.final %>% 
        mutate (
                #randomly split into 4 groups
                placebo.partition = sample(1:4, nrow(A.final), replace = T),
                death.other.1  = if_else (other.cause.mortality == 'Death' & placebo.partition == 1, death, as.Date(NA_character_)),
                death.other.non.1  = if_else (other.cause.mortality == 'Death' & placebo.partition != 1, death, as.Date(NA_character_)),
                death.other.2  = if_else (other.cause.mortality == 'Death' & placebo.partition == 2, death, as.Date(NA_character_)),
                death.other.non.2  = if_else (other.cause.mortality == 'Death' & placebo.partition != 2, death, as.Date(NA_character_)),
                death.other.3  = if_else (other.cause.mortality == 'Death' & placebo.partition == 3, death, as.Date(NA_character_)),
                death.other.non.3  = if_else (other.cause.mortality == 'Death' & placebo.partition != 3, death, as.Date(NA_character_)),
                death.other.4  = if_else (other.cause.mortality == 'Death' & placebo.partition == 4, death, as.Date(NA_character_)),
                death.other.non.4  = if_else (other.cause.mortality == 'Death' & placebo.partition != 4, death, as.Date(NA_character_)),
                )


    A.final %>%   
        mutate(COD = ifelse( COD_TO_SITE_RECODE %in% cod.df$COD_TO_SITE_RECODE, COD_TO_SITE_RECODE, '50300')) %>%
        group_by(COD) %>% summarise(  n = n(), p = n()/nrow(.) ) %>% arrange(desc(n)) %>%   left_join(cod.df, by = c('COD'='COD_TO_SITE_RECODE'))  %>% mutate (np = sprintf('%d (%.1f%%)', n ,p*100 ))%>% select( COD, Name,  np) %>% print(n =10)
    A.final %>% summarise(sum(other.cause.mortality == 'Death'),  mean(other.cause.mortality == 'Death') *100) 

    table( A.final$tx, useNA="ifany")
    # For the sensitivity analysis, some node positive patients are included
    A.final$tnm.n[is.na(A.final$tnm.n)]  <- 'X'
    A.final  <- A.final  %>% filter (tnm.n %in% c('0', '1', '2')) 
    A.final  <- A.final %>% filter (tx == 'sbrt' |
                                    ( tx == 'sublobar' & tnm.n == '0' ) | 
                                    (tx == 'sublobar' & tnm.n != '0' & REGIONAL_NODES_EXAMINED_1988 != '00' ) )
    print(sprintf('%.3f%% of the sublobar patients are N+', 100*sum(A.final$tnm.n != 0 & A.final$tx == 'sublobar')/ sum(A.final$tx == 'sublobar')))
    # print(table( A.final$REGIONAL_NODES_EXAMINED_1988, A.final$tnm.n, useNA="ifany"))
    A.final$tx  <-  droplevels(A.final$tx)

    # preprocessing
    A.final <- A.final %>% mutate( race2 = ifelse (race == 'White' , 'White', 'Other'),
                                  treatment.year2 = as.character(treatment.year),
                                  treatment.year2 = (ifelse (treatment.year2 %in% c('2019', '2020'), '2019_2020', treatment.year2)),
                                  treatment.year2 = (ifelse (treatment.year2 %in% c('2010', '2011'), '2010_2011', treatment.year2)),
                                  histology2 = ifelse (grepl('Adeno', histology), 'Adenocarcinoma', 'Squamous cell'))

    # Define X
    X1.factor  <-  c('sex', 'race2', 'treatment.year2' ) 
    X1.numeric  <-  c('age') 
    X1s  <-  c(sprintf('%s_z', X1.numeric),  X1.factor)

    X2.factor  <-  c('sex', 'race2',  'histology2', 'treatment.year2' ) 
    X2.numeric  <-  c('age', 'size') 
    X2s  <-  c(sprintf('%s_z', X2.numeric),  X2.factor)

    X.numeric  <- unique(c(X1.numeric,X2.numeric))
    X.factor  <- unique(c(X1.factor,X2.factor))

    Z.count  <- c('O2accessories', 'mobility_aids' , 'transportation_services', 'other_supplies',   'pressure_ulcer', 'ischemic_heart_disease', 'CHF', 'PVD', 'CVD',    'LD', 'DIAB_UC', 'DIAB_C',  'RD', 'mental_disorders', 'nervous_system',    'echo','stress_test','pft',   'Anticoags',  'smoking', 'o2',  'pneumonia_and_influenza','asthma','interstitial_lung', 'COPD')
    Z.count.unscaled = sprintf( '%s_pre_12months_unique_count', Z.count )
    Zs  <-   sprintf( '%s_pre_12months_unique_count', Z.count )

    A.final  <- A.final %>% mutate( 
                                   time.offset = pre.tx.months, 
                                   across( all_of(c(Z.count.unscaled)), function(x) (x >0), .names = "{.col}_bool" ),
                                   # across( all_of(c(Z.count.unscaled)), function(x) quartile(x), .names = "{.col}_s" ),
                                   across( all_of(c(X.numeric)), scale_, .names = "{.col}_z" ))

    if (F) {
        print('X1')
        for (i in 1:length(X1s)) cat(i, Xs[i], '\n')
        print('Z')
        for (i in 1:length(Zs)) cat(sprintf('%s\n ', Zs[i]))
    }


    out  <- list(
                 A.final = A.final,
                 W1s = W1s,
                 outcome.names = outcome.names,
                 X1s = X1s,
                 X2s = X2s,
                 X.numeric = X.numeric,
                 X.factor = X.factor,
                 Zs = Zs
    )
    return(out)

}
