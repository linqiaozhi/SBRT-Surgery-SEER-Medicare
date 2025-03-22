two.step <-  function(data.mat, outcome.name,Xs,Zs, W1.selected.variables, Y.count.bool = F, verbose = F, W1  =  'death.other.cause.gt90day'){
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


