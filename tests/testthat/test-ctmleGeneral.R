N <- 100
p = 10
Wmat <- matrix(rnorm(N * p), ncol = p)
beta1 <- 4+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]+2*Wmat[,6]+2*Wmat[,8]
beta0 <- 2+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]+2*Wmat[,6]+2*Wmat[,8]

tau <- 2

gcoef <- matrix(c(-1,-1,rep(-(3/((p)-2)),(p)-2)),ncol=1)

W <- as.matrix(Wmat)

g <- 1/(1+exp(W%*%gcoef /3))
A <- rbinom(N, 1, prob = g)

sigma <- 1
epsilon <-rnorm(N,0,sigma)
Y  <- beta0 + tau * A + epsilon
#----------------------------------------------------------------
#----------------------Test for glmnet CTMLE---------------------
#----------------------------------------------------------------

# cgmleGlmnet must provide user-specified Q
W_tmp <- data.frame(W[,1:3])
treated<- W_tmp[which(A==1),]
untreated<-W_tmp[which(A==0),]
Y1<-Y[which(A==1)]
Y0<-Y[which(A==0)]

#Initial Q-estimate
beta1hat <- predict(lm(Y1~.,
                       data=treated),newdata=W_tmp)
beta0hat <- predict(lm(Y0~.,
                       data=untreated),newdata=W_tmp)

Q <- matrix(c(beta0hat,beta1hat),ncol=2)

W = as.matrix(W)
glmnet_fit <- cv.glmnet(y = A, x = W, family = 'binomial', nlambda = 40)

lambdas <-glmnet_fit$lambda[(which(glmnet_fit$lambda==glmnet_fit$lambda.min)):length(glmnet_fit$lambda)]
gn_candidates <- predict(glmnet_fit, s = lambdas, newx = W)

# No gn_candidate_cv, expect a warning
expect_warning(
      ctmle_general_fit1 <- ctmleGeneral(Y=Y, A=A,
                                         W=data.frame(W=W), Q = Q,
                                         gn_candidates = gn_candidates,
                                         ctmletype=1, alpha=.995,
                                         family="gaussian",
                                         gbound=0.025,like_type="loglik" ,
                                         fluctuation="logistic", verbose=FALSE,
                                         detailed=FALSE, PEN=FALSE,
                                         V=5, stopFactor=10^6)
)

# No corresponding folds, expect a warning

expect_warning(
      ctmle_general_fit2 <- ctmleGeneral(Y=Y, A=A,
                                         W=data.frame(W=W), Q = Q,
                                         gn_candidates = gn_candidates,
                                         gn_candidates_cv = gn_candidates,
                                         ctmletype=1, alpha=.995,
                                         family="gaussian",
                                         gbound=0.025,like_type="loglik" ,
                                         fluctuation="logistic", verbose=FALSE,
                                         detailed=FALSE, PEN=FALSE,
                                         V=5, stopFactor=10^6)
)


gcv <- stats::predict(glmnet_fit, newx=W, s="lambda.min",type="response")
gcv <- bound(gcv,c(0.025,0.975))

s_prev <- glmnet_fit$lambda[(1:which(glmnet_fit$lambda == glmnet_fit$lambda.min))]
gcvPrev <- stats::predict(glmnet_fit,newx=W,s=s_prev[length(s_prev)],type="response")
gcvPrev <- bound(gcvPrev,c(0.025,0.975))

tlme_fit <- tmle::tmle(Y = Y, A = A, W = W, Q = Q, g1W = gcv)

ctmle_general_fit3 <- ctmleGeneral(Y=Y, A=A,
                                   W=W, Q = Q,
                                   ctmletype=3, g1W = gcv, g1WPrev = gcvPrev,
                                   alpha=.995, family="gaussian",
                                   gbound=0.025,like_type="loglik" ,
                                   fluctuation="logistic", verbose=FALSE,
                                   detailed=FALSE, PEN=FALSE,
                                   V=5, stopFactor=10^6)


ctmle_general_fit1
tlme_fit
ctmle_general_fit3


expect_that(ctmle_general_fit1, is_a("ctmle"))
expect_that(ctmle_general_fit2, is_a("ctmle"))
expect_that(ctmle_general_fit3, is_a("ctmle"))
