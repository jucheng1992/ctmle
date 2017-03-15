
<!--
[![](http://cranlogs.r-pkg.org/badges/simcausal)](https://CRAN.R-project.org/package=ctmle) [![](http://cranlogs.r-pkg.org/badges/grand-total/simcausal)](https://CRAN.R-project.org/package=ctmle)
-->
<!-- README.md is generated from README.Rmd. Please edit that file -->
Installation
============

To install the CRAN release version of `ctmle`:

``` r
install.packages('ctmle')
```

To install the development version (requires the devtools package):

``` r
devtools::install_github('jucheng1992/ctmle')
```

Collaborative Targeted Maximum Likelihood Estimation
====================================================

In this package, we implemented the general template of C-TMLE, for estimation of average additive treatment effect (ATE). The package also offers the functions for discrete C-TMLE, which could be used for variable selection, and C-TMLE for model selection of LASSO.

C-TMLE for variable selection
-----------------------------

In this section, we start with examples of discrete C-TMLE for variable selection, using greedy forward searhcing, and scalable discrete C-TMLE with pre-ordering option.

``` r
library(ctmle)
#> Loading required package: SuperLearner
#> Loading required package: nnls
#> Super Learner
#> Version: 2.0-21
#> Package created on 2016-11-11
#> Loading required package: tmle
#> Welcome to the tmle package, version 1.2.0-5
#> 
#> Use tmleNews() to see details on changes and bug fixes
#> Loading required package: glmnet
#> Loading required package: Matrix
#> Loading required package: foreach
#> Loaded glmnet 2.0-5
set.seed(123)

N <- 1000
p = 10
Wmat <- matrix(rnorm(N * p), ncol = p)
beta1 <- 4+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]+2*Wmat[,6]+2*Wmat[,8]
beta0 <- 2+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]+2*Wmat[,6]+2*Wmat[,8]
tau <- 2
gcoef <- matrix(c(-1,-1,rep(-(3/((p)-2)),(p)-2)),ncol=1)
W <- as.matrix(Wmat)

g <- 1/(1+exp(W%*%gcoef /3))
A <- rbinom(N, 1, prob = g)

epsilon <-rnorm(N, 0, 1)
Y  <- beta0 + tau * A + epsilon

# With initial estimate of Q
Q <- cbind(rep(mean(Y[A == 0]), N), rep(mean(Y[A == 1]), N))

time_greedy <- system.time(
      ctmle_discrete_fit1 <- ctmleDiscrete(Y = Y, A = A, W = data.frame(Wmat), Q = Q,
                                           preOrder = FALSE, detailed = TRUE)
)
ctmle_discrete_fit2 <- ctmleDiscrete(Y = Y, A = A, W = data.frame(Wmat),
                                     preOrder = FALSE, detailed = TRUE)


time_preorder <- system.time(
      ctmle_discrete_fit3 <- ctmleDiscrete(Y = Y, A = A, W = data.frame(Wmat), Q = Q,
                                           preOrder = TRUE,
                                           order = rev(1:p), detailed = TRUE)
)


time_greedy
#>    user  system elapsed 
#>   4.733   0.098   4.892
time_preorder
#>    user  system elapsed 
#>   1.855   0.027   1.891

# Show the detailed results from greedy CTMLE
ctmle_discrete_fit1
#> C-TMLE result:
#>  parameter estimate:  1.91732 
#>  estimated variance:  0.00549 
#>             p-value:  <2e-16 
#>   95% conf interval: (1.77215, 2.06249)
summary(ctmle_discrete_fit1)
#> 
#> Number of candidate TMLE estimators created:  11 
#> A candidate TMLE estimator was created at each move, as each new term
#> was incorporated into the model for g.
#> ---------------------------------------------------------------------- 
#>          term added cleverCovar estimate cv-RSS cv-varIC cv-penRSS
#> cand 1  (intercept)           1     3.92   21.4   0.0859     20418
#> cand 2           X2           1     3.07   21.2   0.0911     20249
#> cand 3           X1           1     2.56   21.0   0.0908     19990
#> cand 4           X6           1     2.32   20.8   0.0917     19827
#> cand 5           X8           1     2.13   20.6   0.0909     19689
#> cand 6           X5           1     1.92   20.4   0.0919     19489
#> cand 7           X7           1     1.88   20.5   0.0942     19514
#> cand 8          X10           1     1.86   20.5   0.0941     19536
#> cand 9           X3           1     1.84   20.5   0.0925     19553
#> cand 10          X4           2     1.85   20.5   0.0961     19573
#> cand 11          X9           2     1.86   20.5   0.1017     19569
#> ---------------------------------------------------------------------- 
#> Selected TMLE estimator is candidate 6 
#> 
#> Each TMLE candidate was created by fluctuating the initial fit, Q0(A,W)=E[Y|A,W], obtained in stage 1.
#> 
#>  cand 1: Q1(A,W) = Q0(A,W) + epsilon1a * h1a 
#>              h1a is based on an intercept-only model for treatment mechanism g(A,W)
#> 
#>      cand 2: Q2(A,W) = Q0(A,W) + epsilon1b * h1b 
#>              h1b is based on a treatment mechanism model containing covariates X2
#> 
#>      cand 3: Q3(A,W) = Q0(A,W) + epsilon1c * h1c 
#>              h1c is based on a treatment mechanism model containing covariates X2, X1
#> 
#>      cand 4: Q4(A,W) = Q0(A,W) + epsilon1d * h1d 
#>              h1d is based on a treatment mechanism model containing covariates X2, X1, X6
#> 
#>      cand 5: Q5(A,W) = Q0(A,W) + epsilon1e * h1e 
#>              h1e is based on a treatment mechanism model containing covariates X2, X1, X6, X8
#> 
#>      cand 6: Q6(A,W) = Q0(A,W) + epsilon1f * h1f 
#>              h1f is based on a treatment mechanism model containing covariates X2, X1, X6, X8, X5
#> 
#>      cand 7: Q7(A,W) = Q0(A,W) + epsilon1g * h1g 
#>              h1g is based on a treatment mechanism model containing covariates X2, X1, X6, X8, X5, X7
#> 
#>      cand 8: Q8(A,W) = Q0(A,W) + epsilon1h * h1h 
#>              h1h is based on a treatment mechanism model containing covariates X2, X1, X6, X8, X5, X7, X10
#> 
#>      cand 9: Q9(A,W) = Q0(A,W) + epsilon1i * h1i 
#>              h1i is based on a treatment mechanism model containing covariates X2, X1, X6, X8, X5, X7, X10, X3
#> 
#>      cand 10: Q10(A,W) = Q0(A,W) + epsilon1i * h1i + epsilon2a * h2a                     = Q9(A,W) + epsilon2a * h2a,
#>              h2a is based on a treatment mechanism model containing covariates X2, X1, X6, X8, X5, X7, X10, X3, X4
#> 
#>      cand 11: Q11(A,W) = Q0(A,W) + epsilon1i * h1i + epsilon2b * h2b                     = Q9(A,W) + epsilon2b * h2b,
#>              h2b is based on a treatment mechanism model containing covariates X2, X1, X6, X8, X5, X7, X10, X3, X4, X9
#> 
#> ---------- 
#> C-TMLE result:
#>  parameter estimate:  1.91732 
#>  estimated variance:  0.00549 
#>             p-value:  <2e-16 
#>   95% conf interval: (1.77215, 2.06249)
```

C-TMLE LASSO for model selection of LASSO
-----------------------------------------

In this section, we introduce the C-TMLE algorithms for model selection of LASSO in the estimation of Propensity Score, and for simplicity we call them C-TMLE LASSO algorithm. We have three variacions of C-TMLE LASSO algorithms, see technical details in !!!.

``` r
library(ctmle)
# Generate high-dimensional data
set.seed(123)

N <- 1000
p = 100
Wmat <- matrix(rnorm(N * p), ncol = p)
beta1 <- 4+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]+2*Wmat[,6]+2*Wmat[,8]
beta0 <- 2+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]+2*Wmat[,6]+2*Wmat[,8]
tau <- 2
gcoef <- matrix(c(-1,-1,rep(-(3/((p)-2)),(p)-2)),ncol=1)
W <- as.matrix(Wmat)

g <- 1/(1+exp(W%*%gcoef /3))
A <- rbinom(N, 1, prob = g)

epsilon <-rnorm(N, 0, 1)
Y  <- beta0 + tau * A + epsilon

# With initial estimate of Q
Q <- cbind(rep(mean(Y[A == 0]), N), rep(mean(Y[A == 1]), N))


glmnet_fit <- cv.glmnet(y = A, x = W, family = 'binomial', nlambda = 10)

# We suggest start build a sequence of lambdas from the lambda selected by cross-validation
lambdas <-glmnet_fit$lambda[(which(glmnet_fit$lambda==glmnet_fit$lambda.min)):length(glmnet_fit$lambda)]

# We fit C-TMLE1 algorithm
time_ctmlelasso1 <- system.time(
      ctmle_fit1 <- ctmleGlmnet(Y = Y, A = A,
                                W = data.frame(W = W),
                                Q = Q, lambdas = lambdas, ctmletype=1, 
                                family="gaussian",gbound=0.025, V=5)
)

# We fit C-TMLE2 algorithm
time_ctmlelasso2 <- system.time(
      ctmle_fit2 <- ctmleGlmnet(Y = Y, A = A,
                                W = data.frame(W = W),
                                Q = Q, lambdas = lambdas, ctmletype=2, 
                                family="gaussian",gbound=0.025, V=5)
)
# For C-TMLE3, we need two gn estimators, one with lambda selected by cross-validation, 
# and the other with lambda slightly different from the selected lambda
gcv <- predict.cv.glmnet(glmnet_fit, newx=W, s="lambda.min",type="response")
gcv <- bound(gcv,c(0.025,0.975))

s_prev <- glmnet_fit$lambda[(which(glmnet_fit$lambda == glmnet_fit$lambda.min))] + 1e-4
gcvPrev <- predict.cv.glmnet(glmnet_fit,newx = W,s = s_prev,type="response")
gcvPrev <- bound(gcvPrev,c(0.025,0.975))

time_ctmlelasso3 <- system.time(
      ctmle_fit3 <- ctmleGlmnet(Y = Y, A = A, W = W, Q = Q,
                                ctmletype=3, g1W = gcv, g1WPrev = gcvPrev,
                                family="gaussian",
                                gbound=0.025, V = 5)
)


time_ctmlelasso1
#>    user  system elapsed 
#>   8.927   0.033   9.011
time_ctmlelasso2
#>    user  system elapsed 
#>   9.967   0.074  10.187
time_ctmlelasso3
#>    user  system elapsed 
#>   0.005   0.000   0.005

ctmle_fit1
#> C-TMLE result:
#>  parameter estimate:  2.19459 
#>  estimated variance:  0.10154 
#>             p-value:  5.6933e-12 
#>   95% conf interval: (1.57003, 2.81915)
ctmle_fit2
#> C-TMLE result:
#>  parameter estimate:  2.16923 
#>  estimated variance:  0.05691 
#>             p-value:  <2e-16 
#>   95% conf interval: (1.70164, 2.63681)
ctmle_fit3
#> C-TMLE result:
#>  parameter estimate:  2.0329 
#>  estimated variance:  0.05025 
#>             p-value:  <2e-16 
#>   95% conf interval: (1.59354, 2.47225)

# Show which regularization parameter (lambda) is selected by C-TMLE1:
lambdas[ctmle_fit1$best_k]
#> [1] 0.003751394
```

Advanced topic: the general template of C-TMLE
----------------------------------------------

``` r
library(ctmle)
# Generate high-dimensional data for the general template of C-TMLE
set.seed(123)

N <- 1000
p = 100
Wmat <- matrix(rnorm(N * p), ncol = p)
beta1 <- 4+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]+2*Wmat[,6]+2*Wmat[,8]
beta0 <- 2+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]+2*Wmat[,6]+2*Wmat[,8]
tau <- 2
gcoef <- matrix(c(-1,-1,rep(-(3/((p)-2)),(p)-2)),ncol=1)
W <- as.matrix(Wmat)

g <- 1/(1+exp(W%*%gcoef /3))
A <- rbinom(N, 1, prob = g)

epsilon <-rnorm(N, 0, 1)
Y  <- beta0 + tau * A + epsilon

# With initial estimate of Q
Q <- cbind(rep(mean(Y[A == 0]), N), rep(mean(Y[A == 1]), N))
```

Citation
--------

If you used `ctmle` package in your research, please cite:

> Ju, Cheng; Susan, Gruber; van der Laan, Mark J.; ctmle: Collaborative Targeted Maximum Likelihood Estimation for Variable and Model Selection in Causal Inference

References
----------

### C-TMLE LASSO and C-TMLE for Model Selection

TBD

#### Scalable Discrete C-TMLE with Pre-ordering

> Ju, Cheng; Gruber, Susan; Lendle, Samuel D.; Chambaz, Antoine; Franklin, Jessica M.; Wyss, Richard; Schneeweiss, Sebastian; and van der Laan, Mark J., "Scalable Collaborative Targeted Learning for High-dimensional Data" (June 2016). U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 352. <http://biostats.bepress.com/ucbbiostat/paper352>

#### Discrete C-TMLE with Greedy Search

> Susan, Gruber, and van der Laan, Mark J.. "An Application of Collaborative Targeted Maximum Likelihood Estimation in Causal Inference and Genomics." The International Journal of Biostatistics 6.1 (2010): 1-31.

#### General Template of C-TMLE

> van der Laan, Mark J., and Susan Gruber. "Collaborative double robust targeted maximum likelihood estimation." The international journal of biostatistics 6.1 (2010): 1-71.
