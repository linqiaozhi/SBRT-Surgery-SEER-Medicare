two.step  <-  function(data.mat, outcome.name,Xs,Zs, Y.count.bool = F, verbose = F, W1  =  'death.other.cause.gt90day', offst =T){
    A.temp  <-  data.mat %>% mutate( 
                                    W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt),
                                    W1.time  = if_else (W1.time == 0, 0.5, W1.time)/365,
                                    W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
    )
    A_ = (A.temp$tx == 'sbrt')*1.0
    X_  <-  model.matrix(as.formula(sprintf('~ %s', paste(Xs, collapse = '+'))),  A.temp)[,-1]
    W_ <- as.matrix( A.temp$W1.time)
    D2_  <- A.temp$W1.bool*1.0
    Z_  <- A.temp[,Zs]
    N_  <- dim(A.temp)[1]
    if (offst)  {
        offset_  <- log(A.temp$time.offset)
    }else {
        offset_  <-  rep(0, nrow(A.temp))
    }
    Xw =  list(  as.matrix(cbind(  A_, X_, Z_ )))
    if (Y.count.bool) {
        A.temp  <- A.temp %>% mutate( Y.count  = !!rlang::sym(outcome.name),)
        Y_ = A.temp$Y.count
        out.model  <- p2sls.negbin(Y = Y_, offset= offset_,  
                                 A = A_, 
                                 X = X_, 
                                 W = W_,
				Z = Z_,
                                 Xw = Xw,        
                                 nco_type = c("ah"),
                                 nco_args =  list(list(offset = rep(0, N_), event = D2_)),
                                 verbose =F)
    #TODO: Fix the verbose statement
    }
    if(!Y.count.bool){
        A.temp  <- A.temp %>% mutate(
                                     Y.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                                     Y.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        )
        Y_ = A.temp$Y.time
        D_ = A.temp$Y.bool*1.0
	out.model  <- p2sls.ah(Y = Y_, D = D_,  A = A_, X = X_,
			      W = W_, Z = Z_,
			      Xw = Xw,        
			      nco_type = c("ah"),
			      nco_args = list(list(offset = rep(0, N_), event = D2_)) )
    
    }
    return(out.model)
}

two.step.loglinear  <-  function(data.mat, outcome.name,Xs,Zs, Y.count.bool = F, verbose = F, W1  =  'death.other.cause.gt90day'){
    A.temp  <-  data.mat %>% mutate( 
                                    W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt),
                                    W1.time  = if_else (W1.time == 0, 0.5, W1.time)/365,
                                    W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
    )
    A_ = (A.temp$tx == 'sbrt')*1.0
    X_  <-  model.matrix(as.formula(sprintf('~ %s', paste(Xs, collapse = '+'))),  A.temp)[,-1]
    # browser()
    W_ <- as.matrix( A.temp$W1.time)
    D2_  <- A.temp$W1.bool*1.0
    Z_  <- A.temp[,Zs]
    N_  <- dim(A.temp)[1]
    offset_  <- log(A.temp$time.offset)
    Xw =  list(  as.matrix(cbind(  A_, X_, Z_ )))
    if (Y.count.bool) {
        A.temp  <- A.temp %>% mutate( Y.count  = !!rlang::sym(outcome.name),)
        Y_ = A.temp$Y.count
        out.model  <- p2sls.loglin(Y = Y_,   
                                 A = A_, 
                                 X = X_, 
                                 W = W_,
				Z = Z_,
                                 Xw = Xw,        
                                 nco_type = c("ah"),
                                 nco_args =  list(list(offset = rep(0, N_), event = D2_)),
                                 verbose =F)
    #TODO: Fix the verbose statement
    }
    if(!Y.count.bool){
        A.temp  <- A.temp %>% mutate(
                                     Y.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                                     Y.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        )
        Y_ = A.temp$Y.time
        D_ = A.temp$Y.bool*1.0
	out.model  <- p2sls.ah(Y = Y_, D = D_,  A = A_, X = X_,
			      W = W_, Z = Z_,
			      Xw = Xw,        
			      nco_type = c("ah"),
			      nco_args = list(list(offset = rep(0, N_), event = D2_)) )
    
    }
    return(out.model)
}

two.step.survfuncs  <-  function(data.mat, outcome.name,Xs,Zs, W1.selected.variables, Y.count.bool = F, verbose = F, W1  =  'death.other.cause.gt90day'){
    A.temp  <-  data.mat %>% mutate( 
                                    W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt),
                                    W1.time  = if_else (W1.time == 0, 0.5, W1.time)/365,
                                    W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
    )
    A_ = (A.temp$tx == 'sbrt')*1.0
    X_  <-  model.matrix(as.formula(sprintf('~ %s', paste(Xs, collapse = '+'))),  A.temp)[,-1]
    # browser()
    W_ <- as.matrix( A.temp$W1.time)
    D2_  <- A.temp$W1.bool*1.0
    Z_  <- A.temp[,Zs]
    N_  <- dim(A.temp)[1]
    offset_  <- log(A.temp$time.offset)
    Xw =  list(  as.matrix(cbind(  A_, X_[,colnames(X_) %in% W1.selected.variables], Z_[,colnames(Z_) %in% W1.selected.variables] )))
    A.temp  <- A.temp %>% mutate(
                                 Y.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                                 Y.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
    )
    Y_ = A.temp$Y.time
    D_ = A.temp$Y.bool*1.0
    survfunc_a1  <- p2sls.ah.survfunc( Y = Y_, D = D_,  A = A_, a = 1, X = X_,
                                      W = W_, Z = Z_,
                                      Xw = Xw,        
                                      nco_type = c("ah"),
                                      nco_args = list(list(offset = rep(0, N_), event = D2_)) )
    survfunc_a0  <- p2sls.ah.survfunc( Y = Y_, D = D_,  A = A_, a = 0, X = X_,
                                      W = W_, Z = Z_,
                                      Xw = Xw,        
                                      nco_type = c("ah"),
                                      nco_args = list(list(offset = rep(0, N_), event = D2_)) )
    survfunc_a0$strata <- 'Surgery'
    survfunc_a1$strata <- 'SBRT'
    out.a  <-  rbind (a1 = survfunc_a1, a0 = survfunc_a0)
    return(out.a)
}

two.step.loglinear.variable.selection  <-  function(data.mat, outcome.name,Xs,Zs, W1.selected.variables, Y.count.bool = F, verbose = F, W1  =  'death.other.cause.gt90day'){
    A.temp  <-  data.mat %>% mutate( 
                                    W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt),
                                    W1.time  = if_else (W1.time == 0, 0.5, W1.time)/365,
                                    W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
    )
    A_ = (A.temp$tx == 'sbrt')*1.0
    X_  <-  model.matrix(as.formula(sprintf('~ %s', paste(Xs, collapse = '+'))),  A.temp)[,-1]
    # browser()
    W_ <- as.matrix( A.temp$W1.time)
    D2_  <- A.temp$W1.bool*1.0
    Z_  <- A.temp[,Zs]
    N_  <- dim(A.temp)[1]
    offset_  <- log(A.temp$time.offset)
    Xw =  list(  as.matrix(cbind(  A_, X_[,colnames(X_) %in% W1.selected.variables], Z_[,colnames(Z_) %in% W1.selected.variables] )))
    if (Y.count.bool) {
        A.temp  <- A.temp %>% mutate( Y.count  = !!rlang::sym(outcome.name),)
        Y_ = A.temp$Y.count
        out.model  <- p2sls.loglin(Y = Y_,   
                                 A = A_, 
                                 X = X_, 
                                 W = W_,
				Z = Z_,
                                 Xw = Xw,        
                                 nco_type = c("ah"),
                                 nco_args =  list(list(offset = rep(0, N_), event = D2_)),
                                 verbose =F)
    #TODO: Fix the verbose statement
    }
    if(!Y.count.bool){
        A.temp  <- A.temp %>% mutate(
                                     Y.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt)/365,
                                     Y.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F),
        )
        Y_ = A.temp$Y.time
        D_ = A.temp$Y.bool*1.0
	out.model  <- p2sls.ah(Y = Y_, D = D_,  A = A_, X = X_,
			      W = W_, Z = Z_,
			      Xw = Xw,        
			      nco_type = c("ah"),
			      nco_args = list(list(offset = rep(0, N_), event = D2_)) )
    
    }
    return(out.model)
}

two.step.variable.selection  <-  function(data.mat, 
                                          outcome.name, 
                                          Xw, # Stage 1 matrix of selected A,X,Z for each W
                                          W.selected.Y,
                                          X.selected.Y,
                                          Y.count.bool = F, verbose = F,  skip.W1 = F){
    # print( sprintf('Outcome: %s, with the following Ws:', outcome.name))
    # for (W_i in W.selected.Y) {
        # print( sprintf('-->W: %s', W_i))
    # print(sprintf('-->-->X, Z: %s', paste(colnames(Xw[[W_i]]), collapse = ',')))
    # }
    # Construct the NOC args
    W1  <-  'death.other.cause.gt90day'
    A.temp  <-  data.mat %>% mutate( 
                                    W1.time  = if_else (nna(!!rlang::sym(W1)), as.numeric( !!rlang::sym(W1) - tx.date, units = "days" ), tt),
                                    W1.time  = if_else (W1.time == 0, 0.5, W1.time)/365,
                                    W1.bool = ifelse( nna(!!rlang::sym(W1)), T, F),
    )
    if (length(X.selected.Y) == 0 ) {
        X_  <- NULL
        warning('No X are selected')
    }else{
  X_  <-  model.matrix(as.formula(sprintf('~ %s', paste(X.selected.Y, collapse = '+'))),  A.temp)[,-1]
    }
    A_ = (A.temp$tx == 'sbrt')*1.0
    D2_  <- A.temp$W1.bool*1.0
    N_  <- dim(A.temp)[1]
    # Only select the Xw with names in selected.Ws
    offset_  <- log(A.temp$time.offset)
    if (Y.count.bool) {
        W_ <- cbind( A.temp$W1.time, A.temp[,W.selected.Y])
        Xw_  <- Xw[names(Xw) %in% c(W1,W.selected.Y)]
        A.temp  <- A.temp %>% mutate( Y.count  = !!rlang::sym(outcome.name),)
        Y_ = A.temp$Y.count
        out.model  <- p2sls.negbin(Y = Y_, offset= offset_,  
                                 A = A_, 
                                 X = X_, 
                                 W = W_,
                                 Xw = Xw_,        
                                 nco_type = c("ah", rep("negbin", length(W.selected.Y))),
                                 nco_args = append( list(list(offset = rep(0, N_), event = D2_)),
                                                   replicate (length(W.selected.Y),  list(offset = offset_, init = NA), simplify=F)),
                                 verbose =verbose
        )
    }else{
         A.temp  <- A.temp %>% mutate(
                                      Y.time  = if_else (nna(!!rlang::sym(outcome.name)), as.numeric( !!rlang::sym(outcome.name) - tx.date, units = "days" ), tt),
                                      Y.time  = if_else (Y.time == 0, 0.5, Y.time)/365,
                                      Y.bool = ifelse( nna(!!rlang::sym(outcome.name)), T, F)
         )
         Y_ = A.temp$Y.time
         D_ = A.temp$Y.bool*1.0
         if (skip.W1) {
             print('Skipping W1')
             nco_type =  rep("negbin", length(W.selected.Y))
             nco_args =  replicate (length(W.selected.Y),  list(offset = offset_, init = NA), simplify=F)
             W_ <- cbind( A.temp[,W.selected.Y])
             Xw_  <- Xw[names(Xw) %in% c(W.selected.Y)]
         }else {
             nco_type = c("ah", rep("negbin", length(W.selected.Y)))
             nco_args = append( list(list(offset = rep(0, N_), event = D2_)),  replicate (length(W.selected.Y),  list(offset = offset_, init = NA), simplify=F))
             W_ <- cbind( A.temp$W1.time, A.temp[,W.selected.Y])
             Xw_  <- Xw[names(Xw) %in% c(W1,W.selected.Y)]
         }
         # if(verbose) cat(sprintf('X: %s, W: %s, Xw: %s\n', paste(colnames(X_), collapse = ','), paste(colnames(W_), collapse = ','), paste(names(Xw_), collapse = ',')))
         out.model  <- p2sls.ah(Y = Y_, 
                              D = D_, 
                              A = A_, 
                              X = X_,
                              W = W_,
                              Xw = Xw_,        
                              nco_type = nco_type,
                              nco_args = nco_args,
                              verbose =verbose)
    }
    return(out.model)
}

