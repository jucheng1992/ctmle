
<!-- README.md is generated from README.Rmd. Please edit that file -->
Collaborative Targeted Maximum Likelihood Estimation (C-TMLE)
=============================================================

C-TMLE for variable selection
-----------------------------

We first generate some code for later showcases.

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
```

Here is some examples on discrete C-TMLE for variable selection, with greedy selection and pre-ordering option.

``` r

time_greedy <- system.time(
ctmle_discrete_fit1 <- ctmleDiscrete(Y = Y, A = A, W = data.frame(Wmat), Q = Q,
                                     preOrder = FALSE)
)
ctmle_discrete_fit2 <- ctmleDiscrete(Y = Y, A = A, W = data.frame(Wmat),
                                     preOrder = FALSE, detailed = TRUE)

ctmle_discrete_fit2$candidates
#> $candidates
#> $candidates$terms
#>  [1] "1"   "X7"  "X3"  "X8"  "X6"  "X5"  "X10" "X9"  "X4"  "X2"  "X1" 
#> 
#> $candidates$covar
#>  [1] 1 1 2 2 3 4 5 6 7 8 9
#> 
#> $candidates$epsilon
#> Tx$h[, "hAW"]                                                         
#>  6.987164e-04  7.752665e-04 -4.612903e-05 -7.287867e-05 -1.236389e-04 
#>                                                                       
#> -2.236605e-05 -5.377358e-05  2.395094e-06  1.947543e-05  1.866394e-04 
#>               
#>  6.002121e-05 
#> 
#> $candidates$earlyStop
#> [1] FALSE
#> 
#> 
#> $results.all
#> $results.all$est
#>  [1] 1.988054 1.990549 1.989122 1.988294 1.984429 1.983729 1.982037
#>  [8] 1.982113 1.982734 1.988911 1.990944
#> 
#> $results.all$likelihood
#>  [1] 0.9632404 0.9632261 0.9632256 0.9632249 0.9632201 0.9632197 0.9632183
#>  [8] 0.9632183 0.9632185 0.9632132 0.9632138
#> 
#> $results.all$varDstar
#>  [1] 0.003871799 0.003824677 0.003816245 0.003809873 0.003826331
#>  [6] 0.003838100 0.003860491 0.003893353 0.003963457 0.004063018
#> [11] 0.004187879
#> 
#> $results.all$varIC
#>  [1] 0.003871799 0.003823223 0.003814907 0.003806124 0.003822576
#>  [6] 0.003834103 0.003853560 0.003887060 0.003957072 0.004051785
#> [11] 0.004168582
#> 
#> $results.all$varIC_origscale
#>  [1] 4.250123 4.196800 4.187672 4.178031 4.196090 4.208744 4.230101
#>  [8] 4.266875 4.343728 4.447696 4.575906
#> 
#> $results.all$varDstar_origscale
#>  [1] 4.250123 4.198396 4.189141 4.182146 4.200212 4.213131 4.237710
#>  [8] 4.273783 4.350737 4.460027 4.597088
#> 
#> 
#> $terms
#>  [1] "1"   "X7"  "X3"  "X8"  "X6"  "X5"  "X10" "X9"  "X4"  "X2"  "X1"

time_preorder <- system.time(
ctmle_discrete_fit3 <- ctmleDiscrete(Y = Y, A = A, W = data.frame(Wmat), Q = Q,
                                     preOrder = TRUE,
                                     order = rev(1:p), detailed = TRUE)
)
```

C-TMLE for model selection of LASSO
-----------------------------------

Advanced topic: the general template of C-TMLE
----------------------------------------------
