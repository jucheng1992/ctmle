#' Discrete Collaborative Targeted Minimum-loss based Estimation
#'
#' This function computes the discrete Collaborative Targeted Minimum-loss based Estimator for variable selection.
#'It includes the greedy C-TMLE algorithm (Gruber and van der Laan 2010), and scalable C-TMLE algorithm
#'(Ju, Gruber, and Lendle et al. 2016) with a user-specified order.
#'
#' @export
#' @param Y continuous or binary outcome variable
#' @param A binary treatment indicator, 1 for treatment, 0 for control
#' @param W vector, matrix, or dataframe containing baseline covariates for Q bar
#' @param Wg vector, matrix, or dataframe containing baseline covariates for propensity score model (defaults to W if not supplied by user)
#' @param Q  n by 2 matrix of initial values for Q0W, Q1W in columns 1 and 2, respectively. Current version does not support SL for automatic initial estimation of Q bar
#' @param Qform optional regression formula for estimating initial Q
#' @param SL.library optional vector of prediction algorithms for data adaptive estimation of Q, defaults to glm, and glmnet
#' @param Qbounds bound on initial Y and predicted values for Q.
#' @param cvQinit if TRUE, cross-validate initial values for Q to avoid overfits
#' @param preOrder boolean indicator for using scalable C-TMLE algorithm or not
#' @param order the use-specified order of covariables. Only used when (preOrder = TRUE). If not supplied by user,
#' it would automatically order covariates from W_1 to W_p
#' @param patience a number to stop early when the score in the CV function does not improve after so many covariates. Used only when (preOrder = TRUE)
#' @param folds The list of indices for cross-validation step. We recommend the cv-splits in C-TMLE matchs that in gn_candidate_cv
#' @param V Number of folds. Only used if folds is not specified
#' @param alpha used to keep predicted initial values bounded away from (0,1) for logistic fluctuation, 0.995 (default)
#' @param family family specification for working regression models, generally 'gaussian' for continuous outcomes (default), 'binomial' for binary outcomes
#' @param like_type  'RSS' or 'loglike'. The metric to use for forward selection and cross-validation
#' @param gbound bound on P(A=1|W), defaults to 0.025
#' @param fluctuation 'logistic' (default) or 'linear', for targeting step
#' @param verbose print status messages if TRUE
#' @param PEN boolean. If true, penalized loss is used in cross-validation step
#' @param stopFactor Numerical value with default 1e6. If the current empirical likelihood is stopFactor times larger than the best previous one, the construction would stop
#' @param detailed boolean number. If it is TRUE, return more detailed results
#' @return best_k  the index of estimate that selected by cross-validation
#' @return est estimate of psi_0
#' @return CI  IC-based 95% confidence interval for parameter estimate
#' @return pvalue pvalue for the null hypothesis that Psi = 0
#' @return likelihood sum of squared residuals, based on selected estimator evaluated on all obs or, logistic loglikelihood if like_type != 'RSS'
#' @return varIC empirical variance of the influence curve adjusted for estimation of g
#' @return varDstar empirical variance of the influence curve
#' @return var.psi variance of the estimate
#' @return varIC.cv cross-validated variance of the influence curve
#' @return penlikelihood.cv penalized cross-validated likelihood
#' @return cv.res all cross-validation results for each fold
#' @examples
#' \dontrun{
#'N <- 1000
#'p = 10
#'Wmat <- matrix(rnorm(N * p), ncol = p)
#'beta1 <- 4+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]+2*Wmat[,6]+2*Wmat[,8]
#'beta0 <- 2+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]+2*Wmat[,6]+2*Wmat[,8]
#'tauW <- 2
#'tau <- 2
#'gcoef <- matrix(c(-1,-1,rep(-(3/((p)-2)),(p)-2)),ncol=1)
#'Wm <- as.matrix(Wmat)
#'g <- 1/(1+exp(Wm%*%gcoef))
#'A <- rbinom(N, 1, prob = g)
#'sigma <- 1
#'epsilon <-rnorm(N,0,sigma)
#'Y  <- beta0 + tauW*A + epsilon
#'
#'# Initial estimate of Q
#'Q <- cbind(rep(mean(Y[A == 0]), N), rep(mean(Y[A == 1]), N))
#'
#'# User-suplied initial estimate
#'time_greedy <- system.time(
#'ctmle_discrete_fit1 <- ctmleDiscrete(Y = Y, A = A, W = data.frame(Wmat), Q = Q,
#'                                     preOrder = FALSE)
#')
#'
#'# If there is no input Q, then intial Q would be estimated by SL with Sl.library
#'ctmle_discrete_fit2 <- ctmleDiscrete(Y = Y, A = A, W = data.frame(Wmat),
#'                                     preOrder = FALSE, detailed = TRUE)
#'
#'# scalable C-TMLE with pre-order option; order is user-specified,
#'# If 'order' is  not specified takes order from W1 to Wp.
#'time_preorder <- system.time(
#'ctmle_discrete_fit3 <- ctmleDiscrete(Y = Y, A = A, W = data.frame(Wmat), Q = Q,
#'                                     preOrder = TRUE,
#'                                     order = rev(1:p), detailed = TRUE)
#')
#'
#'# Compare the running time
#'time_greedy
#'time_preorder
#'}
ctmleDiscrete <- function(Y, A, W, Wg=W, Q=NULL,
                  preOrder=FALSE, order = NULL, patience=FALSE,
                  Qbounds=NULL, cvQinit=FALSE, Qform=NULL, SL.library=NULL,
                  alpha=.995, family= "gaussian",
                  gbound=.025,like_type = "RSS",
                  fluctuation="logistic",
                  verbose=FALSE,  detailed=FALSE,
                  PEN=FALSE, V=5, folds = NULL,
                  stopFactor=10^6
                  ) {


      n <- length(Y)
      maptoYstar <- fluctuation=="logistic"

      if(is.null(colnames(W))){
            colnames(W) <-  paste("W",1:NCOL(W), sep="")
      }
      invalid.name <- nchar(colnames(W)) == 0
      if(any(invalid.name)){
            colnames(W)[invalid.name] <- paste(".internalW", which(invalid.name), sep="")
      }

      if(is.null(colnames(Wg))){
            colnames(Wg) <-  paste("W",1:NCOL(Wg), sep="")
      }
      invalid.name <- nchar(colnames(Wg)) == 0
      if(any(invalid.name)){
            colnames(Wg)[invalid.name] <- paste(".internalW", which(invalid.name), sep="")
      }

      if(!is.null(Q)){
            if(NCOL(Q)==2){
                  Q <- cbind(QAW=(1-A)*Q[,1] + A*Q[,2], Q1W=Q[,2], Q0W=Q[,1])
            }}
      Qinit <- stage1(Y,X=cbind(A,W), Qform=Qform, family, Qbounds=Qbounds, Qinit=Q, alpha=alpha,
                      SL.library, cvQinit=cvQinit,maptoYstar=maptoYstar, verbose)

      ##If preOrder is TRUE, then we run cv first to determine the optimal number of covariates to run CTMLE on. Saves Time and allows for stopping early by setting patience parameter
      if(preOrder==TRUE){
            if(is.null(order)){
                  warning("The order is not specified. Naive order starting from W1 to Wp would be used")
            }else{
                  Wg <- Wg[, order]
            }

            if(verbose) {cat("\t\t-----Cross-validation to select the number of variables to run CTMLE on-----\n")}
            cv.res <- cv(Y=Qinit$Y,X=data.frame(A,Wg),
                         Q=Qinit$Q, family=Qinit$family,
                         Qbounds=Qinit$Qbounds, ab=Qinit$ab,
                         like_type=like_type,gbound=gbound,
                         verbose=verbose, PEN=PEN ,
                         V=V, patience=patience, preOrder=preOrder,
                         folds = folds)

            best_k<- cv.res$best_k
            if(cv.res$best_k==1) best_k<- cv.res$best_k + 1 ##we add one because if best_k=1 then stage2 function gets error (requires at least one covariate). So we will construct one more estimator than needed but we'll end up taking the best one based on cv

            #note: we add one to best_k for X input to make sure we take at least one variable. This will calculate one more estimate than needed
            candidates <- stage2(Y=Qinit$Y, X=data.frame(A,Wg)[,1:best_k], Qinit$Q, family=Qinit$family, Qinit$Qbounds,ab=Qinit$ab,
                                 ncandidates=best_k-1, training_set=rep(TRUE,length(Y)),
                                 like_type=like_type,gbound=gbound, verbose=verbose,
                                 stopFactor=stopFactor, preOrder=preOrder)
      }
      if(preOrder==FALSE){
            candidates <- stage2(Y=Qinit$Y, X=data.frame(A,Wg), Qinit$Q,
                                 family=Qinit$family, Qinit$Qbounds,ab=Qinit$ab,
                                 ncandidates=NCOL(Wg), training_set=rep(TRUE,length(Y)),
                                 like_type=like_type,gbound=gbound, verbose=verbose,
                                 stopFactor=stopFactor, preOrder=preOrder)

            if(verbose) {cat("\t\t-----Cross-validation to select the number of variables to run CTMLE on-----\n")}
            cv.res <- cv(Y=Qinit$Y,X=data.frame(A,Wg), Q=Qinit$Q,
                         family=Qinit$family, Qbounds=Qinit$Qbounds,
                         ab=Qinit$ab,  like_type=like_type,gbound=gbound,
                         verbose=verbose, PEN=PEN , V=V,
                         patience=patience, preOrder=preOrder)
      }
      if(verbose){
            cat("\n\n")
            cat("\t terms included in final estimate: ", candidates$terms, "\n")
            cat("\t all estimates: ", candidates$results.all$est, "\n")
            cat("\t best_k: ", cv.res$best_k, "\n")
            cat("\t best estimate: ", candidates$results.all$est[cv.res$best_k], "\n\n")
      }


      Q.init <- Qinit$Q
      if(family == "binomial" | maptoYstar){Q.init <- plogis(Q.init)}

      psi.init <- mean(Q.init[,"Q1W"] - Q.init[,"Q0W"]) * diff(Qinit$ab)
      if(verbose) {
            cat("\t Estimate based on initial, untargeted Q:", round(psi.init,4),"\n")
            cat(paste("\t Cross validation complete.\n"))
      }

      # add varDstar and varIC on the original scale
      candidates$results.all$varIC_origscale <- candidates$results.all$varIC * ifelse(maptoYstar, (diff(Qinit$ab))^2, 1)
      candidates$results.all$varDstar_origscale <- candidates$results.all$varDstar * ifelse(maptoYstar, (diff(Qinit$ab))^2, 1)
      returnVal <- (list(candidates=candidates, best_k=cv.res$best_k,
                         est= candidates$results.all$est[cv.res$best_k],
                         CI = candidates$results.all$est[cv.res$best_k]+c(-1.96,1.96) * sqrt(candidates$results.all$varIC_origscale[cv.res$best_k]/n),
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
      class(returnVal) <- "ctmle"
      return(returnVal)
}







