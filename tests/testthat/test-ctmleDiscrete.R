N <- 100
p = 5

Wmat <- matrix(rnorm(N * p), ncol = p)
beta1 <- 4+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]
beta0 <- 2+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]

tauW <- 2
tau <- 2

gcoef <- matrix(c(-1,-1,rep(-(3/((p)-2)),(p)-2)),ncol=1)

Wm <- as.matrix(Wmat)

g <- 1/(1+exp(Wm%*%gcoef))
A <- rbinom(N, 1, prob = g)

sigma <- 1
epsilon <-rnorm(N,0,sigma)
Y  <- beta0 + tauW*A + epsilon
#----------------------------------------------------------------
#----------------------Test for discrete CTMLE-------------------
#----------------------------------------------------------------
# With initial estimate of Q
Q <- cbind(rep(mean(Y[A == 0]), N), rep(mean(Y[A == 1]), N))

time_greedy <- system.time(
ctmle_discrete_fit1 <- ctmleDiscrete(Y = Y, A = A, W = data.frame(Wmat), Q = Q,
                                     preOrder = FALSE)
)
ctmle_discrete_fit2 <- ctmleDiscrete(Y = Y, A = A, W = data.frame(Wmat),
                                     preOrder = FALSE, detailed = TRUE)

ctmle_discrete_fit2$candidates

time_preorder <- system.time(
ctmle_discrete_fit3 <- ctmleDiscrete(Y = Y, A = A, W = data.frame(Wmat), Q = Q,
                                     preOrder = TRUE,
                                     order = rev(1:p), detailed = TRUE)
)


V = 5
folds <- by(sample(1:N,N), rep(1:V, length=N), list)


ctmle_discrete_fit4 <- ctmleDiscrete(Y = Y, A = A, W = data.frame(Wmat), Q = Q,
                                     preOrder = TRUE, V = 5, folds = folds,
                                     order = 1:p, detailed = TRUE)

# No user-specified order, should be a warning
expect_warning(
ctmle_discrete_fit5 <- ctmleDiscrete(Y = Y, A = A, W = data.frame(Wmat),Q = Q,
                                     preOrder = TRUE, V = 5, folds = folds,
                                     detailed = TRUE)
)

# Actually used same order and cv-folds, should have same estimate
expect_equal(ctmle_discrete_fit4$est, ctmle_discrete_fit5$est)
expect_that(ctmle_discrete_fit4, is_a("ctmle"))
# Greedy takes longer time
expect_true(time_greedy[3] > time_preorder[3])




