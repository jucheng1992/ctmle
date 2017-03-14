
<!-- README.md is generated from README.Rmd. Please edit that file -->
``` r
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
```

Here is examples on discrete C-TMLE for variable selection, with greedy selection and pre-ordering option.
