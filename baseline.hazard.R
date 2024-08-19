################################################################################
#=========== semiparametric additive hazards model =============================
################################################################################
# Semiparametric additive hazards model in Lin & Ying 1997 with time-independent
# covariates
# time: vector of time-to-event
# event: vector of event indicator; 0 = censored
# covariates: design matrix of covariates, intercept should not be included.


#' Additive hazards regression model
#'
#' Fit semiparametric additive hazards model in Lin & Ying 1997 with time-independent
#' covariates.
#' This is a utility function.
#'
#' @param time vector of time-to-event
#' @param event vector of event indicator; 0 = censored
#' @param offset n-vector of offset
#' @param covariates design matrix of covariates; intercept should not be included
#' @param variance whether to return the variance and components for the calculation
#' @examples
#' x <- rep(1:3, each = 10)
#' t_uncensored <- rexp(30, rate = 1 / exp(0.5 -0.1 * x))
#' t_censored <- pmin(5, t_uncensored)
#' event <- as.numeric(t_uncensored <= 5)
#'
#' ah_result <- lin_ah(time = t_censored, event = event, covariates = x)
#' ah_result$summary
#' @export
baseline.cum.hazard <- function(time, event, covariates, ESTIMATE, weights = NULL,
                               offset = rep(0, length(time)), variance = TRUE) {
  t <- as.numeric(time)
  d <- as.numeric(event)
  x <- as.matrix(covariates)
  eta <- as.numeric(offset)

  nn <- length(t) # sample size

  # assign weights
  if (is.null(weights)) weights <- rep(1, nn)

  w <- as.numeric(weights)

  # sort the observations by time in increasing order
  dat1 <- cbind(t, d, w, eta, x)[order(t),]

  t1 <- dat1[, 1]
  d1 <- dat1[, 2]
  w1 <- dat1[, 3]
  eta1 <- dat1[, 4]
  x1 <- dat1[, -(1:4), drop = F]

  t2 <- unique(t1)
  ntime <- length(t2) # length of unique time points

  # obtain the order of time
  o1 <- dplyr::left_join(data.frame(t1 = t1),
                  data.frame(t1 = t2, o1 = 1:length(unique(t1))),
                  by = "t1")$o1

  # the minimum index of those who have t1 equals the ith ordered t1
  tmin <- sapply(1:ntime, function(tt) min(which(o1 == tt)))
  dtime <- c(t2[1], diff(t2))

    # baseline cumulative hazard function
    haz <- 1:ntime * NA

    for (k in 1:ntime) {
      haz[k] <- sum(d1 * (o1 == k) - (eta1 + x1 %*% ESTIMATE) * (o1 >= k) * dtime[k]) /
        sum(o1 >= k)
    }
    cumhaz <- cumsum(haz)


    unadj.haz <- 1:ntime * NA
    for (k in 1:ntime) {
      unadj.haz[k] <- sum(d1 * (o1 == k) ) /
        sum(o1 >= k)
    }
    cum.unadj.haz <- cumsum(unadj.haz)
    return(list(cumhaz = cumhaz,
                cum.unadj.haz = cum.unadj.haz,
                t_=t2))
    
}







