#' General Template for Collaborative Targeted Maximum Likelihood Estimation
#'
#' This function computes the Collaborative Targeted Maximum Likelihood Estimator.
#'
#'@import glmnet
#'@import SuperLearner
#'@import tmle
#'@import stats
#' @export
#' @param Y continuous or binary outcome variable
#' @param A binary treatment indicator, 1 for treatment, 0 for control
#' @param W vector, matrix, or dataframe containing baseline covariates for Q bar
#' @param Wg vector, matrix, or dataframe containing baseline covariates for propensity score model
#' (defaults to W if not supplied by user)
#' @param Q n by 2 matrix of initial values for Q0W, Q1W in columns 1 and 2, respectively.
#' Current version does not support SL for automatic initial estimation of Q bar
#' @param ctmletype 1 or 3. Type of general C-TMLE. Type 1 uses cross-validation to select best gn,
#'  while Type 3 directly solves extra clever covariates.
#' @param gn_candidates matrix or dataframe, each column stand for a estimate of propensity score.
#'  Estimate in the column with larger index should have smaller empirical loss
#' @param gn_candidates_cv matrix or dataframe, each column stand for a the cross-validated estimate.
#' For example, the (i,j)-th element is the predicted propensity score by j-th estimator,
#'  for i-th observation, when it is in the validation set
#' @param folds The list of indices for cross-validation step. We recommend the cv-splits in C-TMLE matchs that in gn_candidate_cv
#' @param V Number of folds. Only used if folds is not specified
#' @param alpha used to keep predicted initial values bounded away from (0,1) for logistic fluctuation,
#'  0.995 (default)
#' @param family family specification for working regression models, generally 'gaussian' for continuous
#' outcomes (default), 'binomial' for binary outcomes
#' @param like_type 'RSS' or 'loglike'. The metric to use for forward selection and cross-validation
#' @param gbound bound on P(A=1|W), defaults to 0.025
#' @param fluctuation 'logistic' (default) or 'linear', for targeting step
#' @param verbose print status messages if TRUE
#' @param PEN boolean. If true, penalized loss is used in cross-validation step
#' @param g1W Only used when type is 3. a user-supplied propensity score estimate.
#' @param g1WPrev Only used when type is 3. a user-supplied propensity score estimate,
#' with small fluctuation compared to g1W.
#' @param stopFactor Numerical value with default 1e6.
#' If the current empirical likelihood is stopFactor times larger than the best previous one,
#' the construction would stop
#' @param detailed boolean number. If it is TRUE, return more detailed results
#' @return best_k  the index of estimate that selected by cross-validation
#' @return est estimate of psi_0
#' @return CI  IC-based 95% confidence interval for parameter estimate
#' @return pvalue pvalue for the null hypothesis that Psi = 0
#' @return likelihood sum of squared residuals, based on selected estimator evaluated on all obs or, logistic loglikelihood if like_type != "RSS"
#' @return varIC empirical variance of the influence curve adjusted for estimation of g
#' @return varDstar empirical variance of the influence curve
#' @return var.psi variance of the estimate
#' @return varIC.cv cross-validated variance of the influence curve
#' @return penlikelihood.cv penalized cross-validated likelihood
#' @return cv.res all cross-validation results for each fold
#' @examples
#'N <- 1000
#'p = 100
#'V = 5
#'Wmat <- matrix(rnorm(N * p), ncol = p)
#'gcoef <- matrix(c(-1,-1,rep(-(3/((p)-2)),(p)-2)),ncol=1)
#'
#'W <- as.data.frame(Wmat)
#'g <- 1/(1+exp(Wmat%*%gcoef / 3))
#'A <- rbinom(N, 1, prob = g)
#'
#'# Build potential outcome pairs, and the observed outcome Y
#'beta1 <- 4+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]+2*Wmat[,6]+2*Wmat[,8]
#'beta0 <- 2+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]+2*Wmat[,6]+2*Wmat[,8]
#'
#'tau = 2
#'sigma <- 1
#'epsilon <-rnorm(N,0,sigma)
#'Y  <- beta0 + tau * A + epsilon
#'# Initial estimate of Q
#'Q <- cbind(rep(mean(Y[A == 1]), N), rep(mean(Y[A == 0]), N))
#'
#'folds <-by(sample(1:N,N), rep(1:V, length=N), list)
#'
#'lasso_fit <- cv.glmnet(x = as.matrix(W), y = A, alpha = 1, nlambda = 100, nfolds = 10)
#'lasso_lambdas <- lasso_fit$lambda[lasso_fit$lambda <= lasso_fit$lambda.min][1:5]

#'# Build template for glmnet
#'SL.glmnet_new <- function (Y, X, newX, family, obsWeights, id, alpha = 1,
#'                           nlambda = 100, lambda = 0,...)
#'{
#'     # browser()
#'     if (!is.matrix(X)) {
#'           X <- model.matrix(~-1 + ., X)
#'          newX <- model.matrix(~-1 + ., newX)
#'    }
#'    fit <- glmnet::glmnet(x = X, y = Y,
#'                          lambda = lambda,
#'                          family = family$family, alpha = alpha)
#'    pred <- predict(fit, newx = newX, type = "response")
#'      fit <- list(object = fit)
#'    class(fit) <- "SL.glmnet"
#'    out <- list(pred = pred, fit = fit)
#'    return(out)
#'}
#'
#'# Use a sequence of LASSO estimators to build gn sequence:
#'SL.cv1lasso <- function (... , alpha = 1, lambda = lasso_lambdas[1]){
#'    SL.glmnet_new(... , alpha = alpha, lambda = lambda)
#'}
#'
#'SL.cv2lasso <- function (... , alpha = 1, lambda = lasso_lambdas[2]){
#'     SL.glmnet_new(... , alpha = alpha, lambda = lambda)
#'}
#'
#'SL.cv3lasso <- function (... , alpha = 1, lambda = lasso_lambdas[3]){
#'     SL.glmnet_new(... , alpha = alpha, lambda = lambda)
#'}
#'
#'SL.cv4lasso <- function (... , alpha = 1, lambda = lasso_lambdas[4]){
#'      SL.glmnet_new(... , alpha = alpha, lambda = lambda)
#'}
#'
#'SL.library = c('SL.cv1lasso', 'SL.cv2lasso', 'SL.cv3lasso', 'SL.cv4lasso', 'SL.glm')
#'
#'# Build the sequence. See more details in the function build_gn_seq, and package SuperLearner
#'gn_seq <- build_gn_seq(A = A, W = W, SL.library = SL.library, folds = folds)
#'
#'
#'# Use the output of build_gn_seq for ctmle general templates
#'ctmle_fit <- ctmleGeneral(Y = Y, A = A, W = W, Q = Q, ctmletype = 1,
#'                      gn_candidates = gn_seq$gn_candidates,
#'                      gn_candidates_cv = gn_seq$gn_candidates_cv,
#'                      folds = gn_seq$folds, V = length(folds))

ctmleGeneral <- function(Y, A, W, Wg = W, Q,
                         ctmletype,
                         gn_candidates,
                         gn_candidates_cv = NULL,
                         alpha=.995, family= "gaussian",
                         gbound=.025, like_type = "RSS",
                         fluctuation="logistic", verbose=FALSE,
                         detailed=FALSE, PEN=FALSE,
                         g1W = NULL, g1WPrev = NULL,
                         V=5, folds = NULL,
                         stopFactor=10^6) {

      if(ctmletype != 3){
            if(is.null(gn_candidates_cv)){
                  gn_candidates_cv <- gn_candidates
                  warning('All gn candidates are only estimated on whole data.
                    Cross-validation might not be objective enough!')
            }else if(is.null(folds)){
                  warning('It is more objective to provide folds, corresponding to gn_candidate_cv.')
            }
      }

      n <- length(Y)
      maptoYstar <- fluctuation=="logistic"
      Qbounds <- range(Y)
      QAW <- (1 - A) * Q[, 1] + A * Q[, 2]
      Q <- cbind(QAW, Q0W = Q[, 1], Q1W = Q[, 2])
      Q <- bound(Q, Qbounds)
      Ystar <- bound(Y, Qbounds)
      ab <- range(Ystar, na.rm = TRUE)
      Ystar[is.na(Ystar)] <- 0
      Ystar <- (Ystar - ab[1])/diff(ab)
      Q <- (Q - ab[1])/diff(ab)
      Qbounds <- c(alpha, 1 - alpha)
      Q <- bound(Q, Qbounds)
      Q <- qlogis(Q)
      family.stage2 <-  "binomial"
      Qinit<-list(Q=Q, Y=Ystar, Qbounds=Qbounds, ab=ab, family=family.stage2)


      if(ctmletype == 3){
            #------------------------------------------------------------------------
            #-----------------------ctmleGeneral type 3-------------------------------
            #------------------------------------------------------------------------
            if(is.null(g1W) | is.null(g1WPrev)){
                  stop("For cglmeGlmnet type 3 needs inute g1W and g1WPrev")
            }

            g <- g1W
            Y=Qinit$Y
            X=data.frame(A,Wg)
            Q=Qinit$Q

            family=Qinit$family
            Qbounds=Qinit$Qbounds
            ab=Qinit$ab
            gbound=gbound

            g1W.total1 <- bound(g1WPrev, c(gbound, 1-gbound))
            g0W.total1 <- 1 - bound(g1WPrev, c(gbound, 1-gbound))
            g1W.total <- bound(g, c(gbound, 1-gbound))
            g0W.total <- 1 - bound(g, c(gbound, 1-gbound))

            H1W <- X[,1]/g1W.total
            H0W <- (1 - X[,1])/g0W.total


            ddg1<-(g1W.total-g1W.total1)
            ddg1[which(ddg1==0)]<-1e-10
            ddg0<-(g0W.total-g0W.total1)
            ddg0[which(ddg0==0)]<-1e-10

            ddH1W <- (X[,1]/(g1W.total^2)) * ddg1
            ddH0W <- ((1 - X[,1])/(g0W.total^2)) * ddg0

            suppressWarnings(epsilon<- coef(glm(Y ~ -1 +
                                                      offset(Q[, "QAW"])+
                                                      H0W+H1W+
                                                      ddH0W + ddH1W, family = family)))

            Qstar <- Q + cbind((epsilon[1] * H0W+ epsilon[2] * H1W),
                               epsilon[1]/g0W.total, epsilon[2]/g1W.total)+
                  cbind((epsilon[3] * ddH0W+ epsilon[4] * ddH1W),
                        (epsilon[3]/(g0W.total^2))*ddg0, (epsilon[4]/(g1W.total^2))*ddg1)

            if(family[1]=="binomial"){
                  Qstar <- plogis(Qstar)
            }

            mu1 <- mean(Qstar[,"Q1W"])
            mu0 <- mean(Qstar[,"Q0W"])
            est  <- (mean(Qstar[,"Q1W"]) - mean(Qstar[,"Q0W"]))*diff(ab)


            Dstar <-  (A/g1W.total - (1 - A)/(1-g1W.total))  * (Y - Qstar[, "QAW"]) +
                  Qstar[, "Q1W"] - Qstar[, "Q0W"] -  (mean(Qstar[,"Q1W"]) - mean(Qstar[,"Q0W"]))


            var.psi <- var(Dstar)/length(Y)  * (diff(ab)^2)
            CI <- est + c(-1.96,1.96) * sqrt(var.psi)

            # !!! Next step: add more details

            emptyInfo <- 'Not Applicable'
            returnVal <- (list(candidates=emptyInfo,
                               best_k=emptyInfo,
                               est= est,
                               CI = CI,
                               Qinit=Qinit,
                               ab = Qinit$ab,
                               # Next step: provide likelihood.
                               likelihood= emptyInfo,
                               varIC = emptyInfo,
                               varDstar = emptyInfo,
                               varIC.cv = emptyInfo,
                               var.psi = var.psi,
                               penlikelihood.cv = emptyInfo,
                               cv.res = emptyInfo)
            )

            returnVal$pvalue <- 2 * pnorm(-abs(returnVal$est)/sqrt(returnVal$var.psi))

            # !!!
            class(returnVal) <- "ctmle"
            return(returnVal)

      }else{
            # Construct by gn estimated with whole data
            # !!Notice here gns are gn_candidates
            candidates <- stage2_general(Y=Qinit$Y, X=data.frame(A,Wg), Q=Qinit$Q,
                                         gn_candidates = gn_candidates,
                                         ctmletype=ctmletype,
                                         family=Qinit$family,
                                         Qbounds=Qinit$Qbounds,ab=Qinit$ab,
                                         # Use whole data as training
                                         training_set=rep(TRUE,length(Y)),
                                         like_type=like_type,gbound=gbound, verbose=verbose,
                                         stopFactor=stopFactor)

            if(verbose) {cat("\t\t-----Cross-validation to select estimator-----\n")}

            # Select best k by CV
            # !!Notice here gns are gn_candidates_cv
            cv.res <- cv_general(Y=Qinit$Y,X=data.frame(A,Wg),
                                 est.all=candidates$results.all$est,
                                 Q=Qinit$Q,
                                 gn_candidates = gn_candidates_cv,
                                 ctmletype=ctmletype,
                                 family=Qinit$family,
                                 Qbounds=Qinit$Qbounds,
                                 ab=Qinit$ab,
                                 like_type=like_type,gbound=gbound,
                                 verbose=verbose, PEN=PEN , V=V,
                                 folds = folds)

            Q.init <- Qinit$Q

            if(family == "binomial" | maptoYstar){Q.init <- plogis(Q.init)}

            psi.init <- mean(Q.init[,"Q1W"] - Q.init[,"Q0W"]) * diff(Qinit$ab)

            if(verbose) {
                  cat("\tEstimate based on initial, untargeted Q:", round(psi.init,4),"\n")
                  cat(paste("\tCross validation complete.\n"))
            }

            candidates$results.all$varIC_origscale <-
                  candidates$results.all$varIC * ifelse(maptoYstar, (diff(Qinit$ab))^2, 1)
            candidates$results.all$varDstar_origscale <-
                  candidates$results.all$varDstar * ifelse(maptoYstar, (diff(Qinit$ab))^2, 1)

            # browser()
            returnVal <- (list(candidates=candidates,
                               best_k=cv.res$best_k,
                               est= candidates$results.all$est[cv.res$best_k],
                               CI = candidates$results.all$est[cv.res$best_k] +
                                     c(-1.96,1.96) * sqrt(candidates$results.all$varIC_origscale[cv.res$best_k]/n),
                               Qinit=Qinit,
                               ab = Qinit$ab,
                               likelihood= candidates$results.all$likelihood[cv.res$best_k],
                               varIC = candidates$results.all$varIC_origscale[cv.res$best_k],
                               varDstar = candidates$results.all$varDstar_origscale[cv.res$best_k],
                               varIC.cv=cv.res$varIC[cv.res$best_k],
                               var.psi=candidates$results.all$varIC_origscale[cv.res$best_k]/n,
                               penlikelihood.cv = cv.res$penlikelihood[cv.res$best_k],
                               cv.res=cv.res))
            returnVal$pvalue <- 2*pnorm(-abs(returnVal$est)/sqrt(returnVal$var.psi))
            just_enough <- c(2:4,6,8,9,11,14)
            if (!detailed) {
                  returnVal <- returnVal[just_enough]
            }
      }
      class(returnVal) <- "ctmle"
      return(returnVal)

}











