library(tidyverse)
library(survival)
cox_nc_est <- function(data, A, X, Z, W, Y, event) {
  A <- data[, A]
  Y <- data[, Y]
  X <- as.matrix(data[, X])
  Z <- data[, Z]
  W <- data[, W]
  event <- data[, event]
  qfunc <- function(tau) {
    c(exp(tau[1] + tau[2] * Z + X %*% tau[3:(2 + ncol(X))]))
  }
  U_q <- function(tau) {
    q_est <- qfunc(tau)
    cbind(1, X, W) * c((1 - A) * q_est - A)
  }
  G <- function(tau) {
    sum(colMeans(U_q(tau)) ^ 2)
  }
  tau_init <- coef(glm(A ~ Z + X, family = binomial))
  nlm_out <- nlm(f = G, p = tau_init, iterlim = 500)
  G(rnorm(19))
  tau_est <- nlm_out$estimate
  # print(sprintf('Number of iterations: %d, cost %f', nlm_out$iter, nlm_out$minimum))
  # print(rbind(tau_init, tau_est))
  q_est <- qfunc(tau_est)
  nc_weights <- q_est ^ (1 - A)
  nc_ph <- coxph(Surv(Y, event) ~ A, weights = nc_weights)
   # glm( A  ~ W[,1] , family = poisson) %>% summary
   # glm( A ~ W[,1], family = poisson, weights = nc_weights) %>% summary
  return(coef(nc_ph)[1])
}

cox_nc <- function(data, A, X, Z, W, Y, event, B = 1000, ncores = 1) {
  est <- cox_nc_est(data, A, X, Z, W, Y, event)
  est_boot <- parallel::mclapply(1:B, function(bb){
    boot_data <- data[sample(nrow(data), replace = T),]
    nc_temp <- try(cox_nc_est(boot_data, A, X, Z, W, Y, event), silent = T)
    if (class(nc_temp) == "try-error") {
      return(NA)
    } else {
      return(nc_temp)
    }
  }, mc.cores = ncores)
  se <- sd(unlist(est_boot), na.rm = T)
  return(list(estimate = est, se = se))
}
make_data <- function(N = 5000) {
  U <- runif(N, 0, 1)
  X1 <- runif(N, 0, 1)
  X2 <- runif(N, 0, 1)
  A <- rbinom(N, 1, 0.2 + 0.6 * U - 0.1 * X1 + 0.2 * X2)
  Z <- rbinom(N, 1, 0.2 + 0.6 * U - 0.2 * X1 + 0.1 * X1)
  TT <- rexp(N, rate = exp(-0.5 + 1 * U - 1.2 * A + 1.2 * X1 - 0.05 * X2))
  CC <- 10
  Y <- round(pmin(TT, CC) * 100) / 100
  event <- as.numeric(TT <= CC)
  W1 <- rpois(N, Y * exp(0.2 + 5 * U))
  W2 <- rpois(N, Y * exp(0.5 + 3 * U))
  return(data.frame(U = U, A = A, Z = Z, Y = Y, W1 = W1, W2 = W2,
                    X1 = X1, X2 = X2, event = event))
}

# df2  <-  make_data(N = 10000)
# # Make a survival plot  using plot_surv
# library(survival)
# library(ggfortify)
# library(ggplot2)
# fit  <-  survfit(Surv(Y, event) ~ A, data = df2)
# autoplot(fit)


# Nits  <-  200
# ests  <-  rep(NA, Nits )



# set.seed(3)
# ests <- parallel::mclapply(1:1000, function(bb){
#                    data2 <- make_data(N = 5000)
#                    nc_temp <- try( cox_nc_est(data2, 'A', c("X1", "X2"), "Z", c("W1", "W2"), "Y", "event"), silent=T)
#                    if (class(nc_temp) == "try-error") {
#                        print(nc_temp)
#                        ests[i] = NA
#                    } else {
#                        ests[i] = nc_temp['A']
#                    }
#   }, mc.cores = 10) %>% unlist

# cor(data2$W1, data2$U)
# cor(data2$Z, data2$U)
# summary(ests)
# quantile(ests, c(0.025, 0.975))

# # Use the function coxph to estimate the causal effect of A on Y
# cox.out  <- coxph(Surv(Y, event) ~ A + X1 + X2 + U, data = df2)
# summary(cox.out)
# cox.out  <- coxph(Surv(Y, event) ~ A + X1 + X2 + W1 + W2, data = df2)
# summary(cox.out)

################################
# Diagnostics 
################################
