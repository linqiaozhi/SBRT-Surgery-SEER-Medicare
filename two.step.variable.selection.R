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
    X_  <-  model.matrix(as.formula(sprintf('~ %s', paste(X.selected.Y, collapse = '+'))),  A.temp)[,-1]
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
