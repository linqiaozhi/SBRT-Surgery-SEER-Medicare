library(dplyr)
# Pair of negative controls
N  <-  1E4
e1  <-  rnorm(N); sig1  <-  2;
e2  <-  rnorm(N); sig2  <-  1
e3  <-  rnorm(N); sig3  <-  1
e4  <-  rnorm(N); sig4  <-  2
e5  <-  rnorm(N); sig5  <-  2
alpha1  <- 1; alpha2  <-  1; alpha3  <- 1; alpha4  <-  1
beta_  <- 1.3 # Causal parameter

U  <-  rnorm (N )
W  <-  alpha1 * U + sig1 * e1 # Negative control exposure
W2  <-  alpha1 * U + sig5 * e5 # Negative control exposure
Z  <-  alpha4 * U + sig4 * e4 # Negative control outcome

A  <-  alpha2*U + sig2*e2     # Treatment
Y  <-  beta_*A + alpha3 * U+ sig3*e3  # Outcome

# Without controlling
lmout <- lm ( Y~ A)
abs(beta_ - (coef(lmout)[['A']]))/ beta_  #0.46

# Naively adding both controls as additional covariates
lmout <- lm ( Y~ A  + W+ Z)
abs(beta_ - (coef(lmout)[['A']]))/ beta_  # 0.31, slight improvemenrt in bias

# Proximal identification using two stage least squares
# step1  <-  lm(W ~ A + Z)
step1  <-  lm(W ~ A + Z)
W_hat  <-  predict(step1)
step2  <-  lm( Y ~ W_hat + A )
abs(beta_ - (coef(step2)[['A']])) / beta_ # 0.02, essentially complete removal of bias


calc.err  <-  function (cox.out, beta_A, var.name = 'A')  100*(coef(cox.out)[var.name] - beta_A)  / beta_A



set.seed(3)
Nits = 200
outs  <-  data.frame()
for (i in 1:Nits) { 
    U  <-  rnorm (N )
    U2  <-  rnorm (N )
    W  <-   U + 4*rnorm(N) 
    W2  <-  U2 + 4*rnorm(N)
    Z  <-   U + 4*rnorm(N)
    A  <-   2*U -3*U2+ 4*rnorm(N)     # Treatment
    beta_  <- 0.3
    Y  <-  beta_*A +U + U2+ rnorm(N) # Outcome
    step1  <-  lm(W ~ A + Z)
    W_hat  <-  predict(step1)
    step1  <-  lm(W2 ~ A + Z)
    W2_hat  <-  predict(step1)
    step2_single  <-  lm( Y ~ A +  W_hat)
    step2  <-  lm( Y ~ A +  W_hat + W2_hat  )
    print(summary(step2_single))
    print(summary(step2))
    outs  <-  rbind (outs, as.data.frame( calc.err(step2, beta_ ) ))
}
out.means  <-  apply( outs , 2, mean)
out.se  <-  apply( outs , 2, FUN = function(x) 1.96*sd(x)/sqrt(length(x)))
out.mae  <-  apply( outs , 2, FUN = function( x) mean(abs(x)))
print(sprintf('Mean %.2f%% +/- %.2f%%, mean absolute error: %.2f ', out.means, out.se, out.mae ))

Z1  <-  cbind( A, Z, W)[1:7,]
Z1.svd  <-  svd(Z1)
Z2  <-  cbind( A, Z, W, 2*A + 3* Z - 1*W)[1:7,]
Z2.svd  <-  svd(Z2)
str(Z1.svd)
str(Z2.svd)
sum(Z1.svd$u[,1:3] != Z2.svd$u[,1:3] )

cor(Z1.svd$u[,1], Z2.svd$u[,1])

hatmat1  <-  Z1.svd$v %*% diag(Z1.svd$d) %*% t(Z1.svd$u)
hatmat2  <-  Z2.svd$v %*% diag(Z2.svd$d) %*% t(Z2.svd$u)

beta_1  <-  hatmat1 %*% Y[1:7]
beta_2  <-  hatmat2 %*% Y[1:7]

fofo  <-  t(Z1.svd$u) %*% 




# Low noise, using both
"Aalen's with two stage: 0.01% +/- 0.03% error"


# High noise, using both
[1] "Aalen's with two stage: -3.24% +/- 2.38% error"

# High noise, using just one
[1] "Aalen's with two stage: -3.62% +/- 2.18% error"


# Tons of noise, just one
[1] "Aalen's with two stage: 37.81% +/- 2.40% error"




head(outs)


cor(W_hat, U)
cor(W, U)
1/sqrt(1+4)
var(W)



library(survival)
lung$age2  <-  lung$age + lung$sex
lfit <- aareg(Surv(time, status) ~ age + sex + ph.ecog + age2, data=lung, nmin=1)
