
<!-- README.md is generated from README.Rmd. Please edit that file -->
Collaborative Targeted Maximum Likelihood Estimation
====================================================

C-TMLE for variable selection
-----------------------------

We first generate some code for later showcases.

``` r
library(ctmle)
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
