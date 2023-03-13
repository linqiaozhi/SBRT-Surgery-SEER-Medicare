# Pair of negative controls
N  <-  5E5
e1  <-  rnorm(N); sig1  <-  2;
e2  <-  rnorm(N); sig2  <-  1
e3  <-  rnorm(N); sig3  <-  1
e4  <-  rnorm(N); sig4  <-  2
alpha1  <- 1; alpha2  <-  1; alpha3  <- 1; alpha4  <-  1
beta_  <- 1.3 # Causal parameter

U  <-  rnorm (N )
W  <-  alpha1 * U + sig1 * e1 # Negative control exposure
Z  <-  alpha4 * U + sig4 * e4 # Negative control outcome

A  <-  alpha2*U + sig2*e2     # Treatment
Y  <-  beta_*A + alpha3 * U+ sig3*e3  # Outcome

# Without controlling
lmout <- lm ( Y~ A)
abs(beta_ - (coef(lmout)[['A']]))/ beta_  #0.46

# Naively adding both controls as additional covariates
lmout <- lm ( Y~ A  + W + Z)
abs(beta_ - (coef(lmout)[['A']]))/ beta_  # 0.31, slight improvemenrt in bias

# Proximal identification using two stage least squares
step1  <-  lm(W ~ A + Z)
W_hat  <-  predict(step1)
step2  <-  lm( Y ~ W_hat + A )
abs(beta_ - (coef(step2)[['A']])) / beta_ # 0.02, essentially complete removal of bias


cor(W_hat, U)
cor(W, U)
1/sqrt(1+4)
var(W)

