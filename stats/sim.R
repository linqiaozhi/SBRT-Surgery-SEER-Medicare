set.seed(2)
library(dplyr)
# Pair of negative controls
N  <-  1E4
e1  <-  rnorm(N); sig1  <-  1;
e2  <-  rnorm(N); sig2  <-  1
e3  <-  rnorm(N); sig3  <-  1
e4  <-  rnorm(N); sig4  <-  1
e5  <-  rnorm(N); sig5  <-  1
alpha1  <- 1; alpha2  <-  1; alpha3  <- 1; alpha4  <-  1
beta_  <- 1.3 # Causal parameter

U  <-  rnorm (N )
Usvd  <- svd(cbind(U))
M_U = diag(N) - Usvd$u %*% t(Usvd$u)
# orthogonal_noise = M_U%*% rnorm(N)
W  <-  alpha1 * U + sig1 * e1 # Negative control exposure
 # W2  <-  alpha1 * U + sig5 *orthogonal_noise/norm(orthogonal_noise, type='2')*norm(e1, type='2') # Negative control exposure
 W2  <-  alpha1 * U + rnorm(N)
Z  <-  alpha4 * U + sig4 * e4 # Negative control outcome

A  <-  alpha2*U + sig2*e2     # Treatment
# Y  <-  beta_*A + alpha3 * U+ sig3*e3  # Outcome
 Y  <-  beta_*A + alpha3 * U

# Without controlling
lmout <- lm ( Y~ A)
abs(beta_ - (coef(lmout)[['A']]))/ beta_  #0.38

lmout <- lm ( Y~ e2)
abs(beta_ - (coef(lmout)[['A']]))/ beta_  #0.38

# Oracle
lmout <- lm ( Y~ A + U)
abs(beta_ - (coef(lmout)[['A']]))/ beta_  #0.004


# e2svd  <- svd(cbind(e2))
# M2 = diag(N) - e2svd$u %*% t(e2svd$u)
# W3  <-  alpha1 * U + M2%*%rnorm(N)
# cor(W2,e2)
# # Just with W
# lmout <- lm ( Y~ A + W2)
# abs(beta_ - (coef(lmout)[['A']]))/ beta_  #0.25

################################
# Begin actual experiments 
################################
make.M  <- function (mat) {
    svdout  <- svd(mat)
    return(diag(N) - svdout$u %*% t(svdout$u))
}

# Proximal identification using W 
# step1  <-  lm(W ~ A + Z)
set.seed(3)
O  <-  rnorm(N) # Negative control exposure
step1_O  <-  lm(U ~ A + W-1)
W_hat_O  <-  predict(step1_O)
step2_O  <-  lm( Y ~ W_hat_O + A )
abs(beta_ - (coef(step2_O)[['A']])) / beta_ 
cor(make.M(cbind(W_hat_O))%*%A, U)
# cor(make.M(cbind(A+O))%*%A, U)

# Proximal identification using W 
# step1  <-  lm(W ~ A + Z)
set.seed(3)
O  <-  rnorm(N)*0.01 # Negative control exposure
step1_O  <-  lm(U ~ W + Z-1)
W_hat_O  <-  predict(step1_O)
step2_O  <-  lm( Y ~ W_hat_O + A )
abs(beta_ - (coef(step2_O)[['A']])) / beta_ 
cor( make.M(cbind(W_hat_O))%*%A, U)

#U Is this 10,000 dimensional thing that has been corrupted by noise. I am now projecting it in a certain direction. 
# U1 is if I project it toward A. 
# U2 is if I project it in a random direction.
# U1 is better because it is in the general direction of S. But why does that matter?
# We are trying to remove U from S+U. Why do I want it to be more toward S?
# IT is not because there is U in the direction of S. 
# A is a corrupted version of U. Just like W. Now, what's the best way to remove U from A? 
# U1 is almost colinear with A.k
# Yeah this is not a good model. Ther eason why <MU1, A> = 0 is just because U1 is almost colinear with A.


# WAnt to project A to be orthogonal U. But we don't have U. We have two options, to project into the complement of span (W,V), or to project to the complement of W = \alpha A + \beta Z
# Why is the latter so much better? Well, because it is in the direction of A.
# Think of the cone of error around U. A is in there. So is W, and Z. 
# Now, we want to set the component of A toward U to be zero. That is,we want it to be orthogonal to U.
# Which direction to we project away from? 
# \hat{W} is the closest point to W on the plane span(A,Z)
# \hat{W2} is the closest point to W on the plane (W2,Z)
# The goal is to remove U from A. If I had just U, that would be fine. But I don't. I have A and noisy version of U.
# Instead of removing a random noisy version of U from A, I can do better. I can project a noisy version of U onto the plane spanned by A and Z. 
# This is much closer to the plane spanned by A and U.
# That is the direction I care about. In other words, if I can choose from the error vectors in that cone, I want one that is in the diection of A, because that's all I care about.
# Consider the plane made by span (U, A). 
# \hat{W} is between A and U. 


set.seed(3)
O  <-  rnorm(N) # Negative control exposure
step1_O  <-  lm(W ~ A + Z-1)
W_hat_O  <-  predict(step1_O)
step2_O  <-  lm( Y ~ W_hat_O + A )
abs(beta_ - (coef(step2_O)[['A']])) / beta_ 
fofo  <- make.M(cbind(W_hat_O))%*%A
cor(make.M(cbind(W_hat_O))%*%A, U)
cor(make.M(cbind(W_hat_O))%*%W2, U)
cor( make.M(cbind(W_hat_O))%*%Z, U)

# Proximal identification using W 
# step1  <-  lm(W ~ A + Z)
set.seed(3)
# O  <-  rnorm(N)*0.01 # Negative control exposure
step1_O  <-  lm(W ~ W2 + Z-1)
W_hat_O  <-  predict(step1_O)
step2_O  <-  lm( Y ~ W_hat_O + A )
abs(beta_ - (coef(step2_O)[['A']])) / beta_ 
fofo  <- make.M(cbind(W_hat_O))%*%A
cor(fofo, U)

cor( make.M(cbind(W_hat_O))%*%Z, U)
cor(fofo, U)



################################
# End 
################################

# Of all the vectors that are within this cone of error around U, you want the one that is closest to A. But why? 


# Naively adding both controls as additional covariates
lmout <- lm ( Y~ A  + W+ Z)
abs(beta_ - (coef(lmout)[['A']]))/ beta_  # 0.31, slight improvemenrt in bias

# Proximal identification using two stage least squares
# step1  <-  lm(W ~ A + Z)
step1  <-  lm(W ~ A + Z -1)
W_hat  <-  predict(step1)
step2  <-  lm( Y ~ W_hat + A )
abs(beta_ - (coef(step2)[['A']])) / beta_ # 0.02, essentially complete removal of bias



# Proximal identification using two stage least squares, where we use U instead of Z
# step1  <-  lm(W ~ A + Z)
step1  <-  lm(W ~ Z + A -1)
W_hat  <-  predict(step1)
step2  <-  lm( Y ~ W_hat + A )
abs(beta_ - (coef(step2)[['A']])) / beta_ # 0.02, essentially complete removal of bias



# Proximal identification using two stage least squares, where we use U instead of W
# step1  <-  lm(W ~ A + Z)
step1  <-  lm(U ~ A + Z -1)
W_hat  <-  predict(step1)
step2  <-  lm( Y ~ W_hat + A )
abs(beta_ - (coef(step2)[['A']])) / beta_ # 0.02, essentially complete removal of bias



# Proximal identification using two stage least squares, where we use U instead of W
step1  <-  lm(U ~ A + O -1)
W_hat  <-  predict(step1)
step2  <-  lm( Y ~ W_hat + A )
summary(step2)
sum(resid(step2)^2)
step1  <-  lm(U ~ W + O -1)
W_hat  <-  predict(step1)
step2  <-  lm( Y ~ W_hat + A )
summary(step2)
sum(resid(step2)^2)
step2  <-  lm( Y ~ U + A )
summary(step2)
sum(resid(step2)^2)
step2  <-  lm( Y ~ e2 )
summary(step2)
sum(resid(step2)^2)
# It's a worse fit when you use the W_hat. But it's the truth. 

abs(beta_ - (coef(step2)[['A']])) / beta_ # 0.02, essentially complete removal of bias



################################
# Denoising experiments 
################################
# Denoising experiments


# Project span of W,Z
WZsvd  <- svd(cbind(W, Z))
M_WZ = diag(N) - WZsvd$u %*% t(WZsvd$u)


# Two stage
What = predict(lm(W ~ A + Z))
Whatsvd  <- svd(cbind(What))
M_What = diag(N) - Whatsvd$u %*% t(Whatsvd$u)
Ynew = M_What%*%Y
Anew = M_What%*%A
summary(lm(Anew ~ U))
summary(lm(Ynew ~ U))
summary(lm ( Ynew ~Anew))
beta_


cor(e2,U)

# Waht's the differnce between projecting to space orthogonal to What and space othogonal to span(W,Z)?







# Two stage with O
O  <-  rnorm(N)*0.01 # Negative control exposure
WhatO = predict(lm(W ~ A + O-1))
# WhatO = predict(lm(W ~ A))
WhatOsvd  <- svd(cbind(WhatO))
MO_What = diag(N) - WhatOsvd$u %*% t(WhatOsvd$u)
cor(MO_What%*%A, U)



WhatA = predict(lm(W ~ A -1))
WhatAsvd  <- svd(cbind(WhatA))
MA_What = diag(N) - WhatAsvd$u %*% t(WhatAsvd$u)

WhatAU = predict(lm(U ~ A +O-1))
cor(WhatAU, e2)
WhatAU = predict(lm(U ~ W +O-1))
cor(WhatAU, e2)

WhatAUsvd  <- svd(cbind(WhatAU))
MAU_What = diag(N) - WhatAUsvd$u %*% t(WhatAUsvd$u)
cor(e2, U)
cor(MAU_What%*%A, U)
cor(MAU_What%*%A, e2)
cor(MAU_What%*%Y, e2)

cor(MAU_What%*%A, A)
summary(lm(MAU_What%*%A ~ U))

# My hypothesis is that it somehow removes more U that is shared with e2.
# The answer is in e2. It must be. 




# evaluate
cor(A, U)
cor(M_What%*%A, U)


# evaluate
cor(M_WZ%*%A, U)
cor(A, U)
cor(M_What%*%A, U)


cor(M_U%*%Y, U)
cor(M_WZ%*%Y, U)
cor(M_What%*%Y, U)


cor(M_U%*%Y, A)
cor(M_WZ%*%Y, A)
cor(M_What%*%Y, A)

cor(M_U%*%Y, M_U%*%A)
cor(M_WZ%*%Y, M_WZ%*%A)
cor(M_What%*%Y, M_What%*%A)


# % If you project U onto A vs if you project U onto W, how different is the relationship with A? Obviously big. in the former, it has 1. In the latter it will be much smaller.

# These qre equally good denoising procedures.
cor (What, U)
cor (WZsvd$u %*% t(WZsvd$u) %*% U, U)
# Why is the first so much better at debiasing?
norm(WZsvd$u %*% t(WZsvd$u) %*% U- U, type = '2')
norm (What- U, type = '2')

cor(U,A)
A.noU  <- M_What %*% A
summary(lm(A.noU ~ U ))
cor(M_What%*%U, A)

# evaluate
ur  <- rnorm(N)
ur.svd  <- svd(ur)
M_ur = diag(N) - ur.svd$u %*% t(ur.svd$u)

cor(M_ur%*%A, U)
cor(M_ur%*%U, U)


cor(M_What%*%Y, U)
cor(MO_What%*%Y, U)

cor(M_What%*%Z, U)
cor(M_What%*%W2, U)

summary(lm(U ~  Z + W + A))
summary(lm(A ~  -1+W))
A %*% U
angle_between <- function(a,b) acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )

angle_between(A,U)

angle_between(rnorm(N),rnorm(N))

summary(lm(W ~ A + Z))


# We don't need all of $U$. We nly need the part that is correlated with A. 
# Write 





# ENd denoising experiments



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
