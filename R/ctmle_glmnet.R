#' Collaborative Targeted Maximum Likelihood Estimation for hyper-parameter tuning of LASSO
#'
#' This function computes the Collaborative Maximum Likelihood Estimation for hyper-parameter tuning of LASSO.
#'
#' @export
#' @param Y continuous or binary outcome variable
#' @param A binary treatment indicator, 1 for treatment, 0 for control
#' @param W vector, matrix, or dataframe containing baseline covariates for Q bar
#' @param Wg vector, matrix, or dataframe containing baseline covariates for propensity score model (defaults to W if not supplied by user)
#' @param Q  n by 2 matrix of initial values for Q0W, Q1W in columns 1 and 2, respectively. Current version does not support SL for automatic initial estimation of Q bar
#' @param ctmletype 1, 2 or 3. Type of general C-TMLE. Type 1 uses cross-validation to select best gn, Type 3 directly solves extra clever covariates,
#' and Type 2 uses both cross-validation and extra covariate. See more details in !!!
#' @param lambdas numeric vector of lambdas (regularization parameter) for glmnet estimation of propensity score, with decreasing order. We recommend the
#' first lambda is selected by external cross-validation.
#' @param folds The list of indices for cross-validation step. We recommend the cv-splits in C-TMLE matchs that in gn_candidate_cv
#' @param V Number of folds. Only used if folds is not specified
#' @param alpha used to keep predicted initial values bounded away from (0,1) for logistic fluctuation, 0.995 (default)
#' @param family family specification for working regression models,
#' generally 'gaussian' for continuous outcomes (default), 'binomial' for binary outcomes
#' @param like_type 'RSS' or 'loglike'. The metric to use for forward selection and cross-validation
#' @param gbound bound on P(A=1|W), defaults to 0.025
#' @param fluctuation 'logistic' (default) or 'linear', for targeting step
#' @param verbose print status messages if TRUE
#' @param PEN boolean. If true, penalized loss is used in cross-validation step
#' @param g1W Only used when type is 3. a user-supplied propensity score estimate.
#' @param g1WPrev Only used when type is 3. a user-supplied propensity score estimate, with small fluctuation compared to g1W.
#' @param stopFactor Numerical value with default 1e6. If the current empirical likelihood is stopFactor times larger than the best previous one, the construction would stop
#' @param detailed boolean number. If it is TRUE, return more detailed results
#' @return best_k  the index of estimate that selected by cross-validation
#' @return est estimate of psi_0
#' @return CI  IC-based 95% confidence interval for parameter estimate
#' @return pvalue pvalue for the null hypothesis that Psi = 0
#' @return likelihood sum of squared residuals, based on selected estimator evaluated on all obs or,
#'  logistic loglikelihood if like_type != 'RSS'
#' @return varIC empirical variance of the influence curve adjusted for estimation of g
#' @return varDstar empirical variance of the influence curve
#' @return var.psi variance of the estimate
#' @return varIC.cv cross-validated variance of the influence curve
#' @return penlikelihood.cv penalized cross-validatedlikelihood
#' @return cv.res all cross-validation results for each fold
#' @examples
#' \dontrun{
#' set.seed(123)
#' N <- 1000
#' p = 10
#' Wmat <- matrix(rnorm(N * p), ncol = p)
#' beta1 <- 4+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]+2*Wmat[,6]+2*Wmat[,8]
#' beta0 <- 2+2*Wmat[,1]+2*Wmat[,2]+2*Wmat[,5]+2*Wmat[,6]+2*Wmat[,8]
#' tau <- 2
#' gcoef <- matrix(c(-1,-1,rep(0,(p)-2)),ncol=1)
#' Wm <- as.matrix(Wmat)
#' g <- 1/(1+exp(Wm%*%gcoef / 3))
#' A <- rbinom(N, 1, prob = g)
#' sigma <- 1
#' epsilon <-rnorm(N,0,sigma)
#' Y  <- beta0 + tau * A + epsilon
#' # ctmleGlmnet must provide user-specified Q
#' W_tmp <- data.frame(Wm[,1:3])
#' treated<- W_tmp[which(A==1),]
#' untreated<-W_tmp[which(A==0),]
#' Y1<-Y[which(A==1)]
#' Y0<-Y[which(A==0)]
#' # Initial Q-estimate
#' beta1hat <- predict(lm(Y1~.,data=treated),newdata=W_tmp)
#' beta0hat <- predict(lm(Y0~., data=untreated),newdata=W_tmp)
#' Q <- matrix(c(beta0hat,beta1hat),ncol=2)
#' W = Wm
#' glmnet_fit <- cv.glmnet(y = A, x = Wm,
#'                        family = 'binomial', nlambda = 40)
#' start = which(glmnet_fit$lambda==glmnet_fit$lambda.min))
#' end = length(glmnet_fit$lambda)
#' lambdas <-glmnet_fit$lambda[start:end]
#' ctmle_fit1 <- ctmleGlmnet(Y=Y, A=A,
#'                          W=data.frame(W=W),
#'                          Q = Q, lambdas = lambdas,
#'                          ctmletype=1, alpha=.995,
#'                          family="gaussian",
#'                          gbound=0.025,like_type="loglik" ,
#'                          fluctuation="logistic",
#'                          verbose=FALSE,
#'                          detailed=FALSE, PEN=FALSE,
#'                          V=5, stopFactor=10^6)
#'}


ctmleGlmnet <- function(Y, A, W, Wg = W,
                        Q, lambdas = NULL,
                        ctmletype,
                        V=5, folds = NULL,
                        alpha=.995, family= "gaussian",
                        gbound=.025, like_type = "RSS",
                        fluctuation="logistic", verbose=FALSE,
                        detailed=FALSE, PEN=FALSE,
                        g1W = NULL, g1WPrev = NULL,
                        stopFactor=10^6) {

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
            #-----------------------ctmleGlmnet type 3-------------------------------
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


            var.psi <- var(Dstar)/length(Y) * (diff(ab)^2)
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
                               cv.res = emptyInfo
            )
            )

            returnVal$pvalue <- 2*pnorm(-abs(returnVal$est)/sqrt(returnVal$var.psi))

            # !!!
            class(returnVal) <- "ctmle"
            return(returnVal)

      }else{

            #------------------------------------------------------------------------
            #-----------------------ctmleGlmnet type 1 or 2--------------------------
            #------------------------------------------------------------------------
            if(is.null(lambdas)){
                  stop('For ctmleGlmnet type 1 and 2, lambdas must be provided')
            }


            candidates <- stage2_glmnet(Y=Qinit$Y, X=data.frame(A,Wg), Q=Qinit$Q,
                                        lambdas = lambdas,ctmletype=ctmletype,
                                        family=Qinit$family,
                                        Qbounds=Qinit$Qbounds,ab=Qinit$ab,
                                        training_set=rep(TRUE,length(Y)),
                                        like_type=like_type,gbound=gbound, verbose=verbose,
                                        stopFactor=stopFactor)

            if(verbose) {cat("\t\t-----Cross-validation to select estimator-----\n")}


            cv.res <- cv_glmnet(Y=Qinit$Y,X=data.frame(A,Wg),
                                est.all=candidates$results.all$est,
                                Q=Qinit$Q,lambdas=lambdas,ctmletype=ctmletype,
                                family=Qinit$family,
                                Qbounds=Qinit$Qbounds, ab=Qinit$ab,
                                like_type=like_type,gbound=gbound,
                                verbose=verbose, PEN=PEN , V=V, folds = folds)
            Q.init <- Qinit$Q


            if(family == "binomial" | maptoYstar){Q.init <- plogis(Q.init)}

            psi.init <- mean(Q.init[,"Q1W"] - Q.init[,"Q0W"]) * diff(Qinit$ab)

            if(verbose) {
                  cat("\tEstimate based on initial, untargeted Q:", round(psi.init,4),"\n")
                  cat(paste("\tCross validation complete.\n"))
            }

            #-------------------------------------------------------------

            # browser()
            if(ctmletype == 1){

                  Y=Qinit$Y
                  X=data.frame(A,Wg)
                  family=Qinit$family
                  Qbounds=Qinit$Qbounds
                  ab=Qinit$ab
                  gbound=gbound
                  # Cheng: This extra step guarantee the critival equation is solved!
                  # Modify the stage2_glmnet function to return the corresponding Qstar
                  candidates_detail <- stage2_glmnet(Y=Qinit$Y, X=data.frame(A,Wg), Q=Qinit$Q,
                                                     lambdas=lambdas,ctmletype=ctmletype,
                                                     family=Qinit$family,
                                                     Qbounds=Qinit$Qbounds,ab=Qinit$ab,
                                                     training_set=rep(TRUE,length(Y)),
                                                     like_type=like_type,gbound=gbound, verbose=verbose,
                                                     stopFactor=stopFactor,
                                                     best_k = cv.res$best_k)

                  Qstar_best <- candidates_detail$results.all$Qstar_best
                  # Print the best result before solving critical equation!
                  # mean(Qstar_best[,3] - Qstar_best[,2]) * diff(ab)
                  #------------------------------------------------

                  lambdas_remain <- lambdas[lambdas < lambdas[cv.res$best_k]]

                  f1 <- eval(paste("A ~ ", paste(paste(names(X[,-1]), sep=''), collapse=" + ")))
                  f2 <- as.formula(f1)
                  lassox <- model.matrix(f2, X)[,-1]
                  lassoy <-as.factor(X[,1])

                  glmnet_remain <- glmnet(x=lassox,y=lassoy,family="binomial", lambda= lambdas_remain)
                  gn_remain <- predict(glmnet_remain, newx = lassox,
                                       type="response")

                  gn_remain <-bound(gn_remain, c(gbound, 1- gbound))

                  finalTMLE <- function(g_tmp, detail = FALSE){

                        g1W.total <- bound(g_tmp, c(gbound, 1-gbound))
                        g0W.total <- 1 - bound(g_tmp, c(gbound, 1-gbound))

                        H1W <- X[,1]/g1W.total
                        H0W <- (1 - X[,1])/g0W.total

                        suppressWarnings(glm_tmp <- glm(Y ~ -1 + offset(qlogis(Qstar_best[, "QAW"]))
                                                        + H0W+H1W, family = 'binomial'))
                        epsilon<- coef(glm_tmp)

                        Qstar_tmp <- qlogis(Qstar_best) + cbind((epsilon[1] * H0W+ epsilon[2] * H1W),
                                                                epsilon[1]/g0W.total, epsilon[2]/g1W.total)


                        Qstar_tmp <- plogis(Qstar_tmp)


                        mu1 <- mean(Qstar_tmp[,"Q1W"])
                        mu0 <- mean(Qstar_tmp[,"Q0W"])
                        est  <- (mean(Qstar_tmp[,"Q1W"]) - mean(Qstar_tmp[,"Q0W"]))*diff(ab)

                        Dstar <-  (A/g1W.total - (1 - A)/(1-g1W.total))  * (Y - Qstar_tmp[, "QAW"]) +
                              Qstar_tmp[, "Q1W"] -Qstar_tmp[, "Q0W"] -
                              (mean(Qstar_tmp[,"Q1W"]) - mean(Qstar_tmp[,"Q0W"]))
                        # mean(Dstar)
                        # browser()
                        if(!detail){
                              return(rbind(est
                                           , -logLik(glm_tmp)[1]))
                        }else{

                              var.psi_tmp <- var(Dstar)/length(Y) *(diff(ab)^2)
                              CI <- est + c(-1.96,1.96) * sqrt(var.psi_tmp)
                              CI
                              return(list(est, CI, var.psi = var.psi_tmp))
                        }
                  }

                  tmle_loss <- t(apply(gn_remain, 2, finalTMLE))[,2]
                  # !! Find the minimizer of the empirical loss
                  ind_final <- which.min(tmle_loss)

                  res <- finalTMLE(gn_remain[, ind_final], detail = TRUE)

                  #                   # Estimate that without last extra fluctuation
                  #                   #-------------------------------------------------------
                  #                   est_old <- mean(Qstar_best[,3] - Qstar_best[,2]) * diff(ab)
                  #                   candidates$results.all$varIC_origscale <- candidates$results.all$varIC * ifelse(maptoYstar, (diff(Qinit$ab))^2, 1)
                  #                   candidates$results.all$varDstar_origscale <- candidates$results.all$varDstar * ifelse(maptoYstar, (diff(Qinit$ab))^2, 1)
                  #                   CI_old <- candidates$results.all$est[cv.res$best_k]+c(-1.96,1.96) * sqrt(candidates$results.all$varIC_origscale[cv.res$best_k]/n)
                  #                   #-------------------------------------------------------
                  #                   returnVal <- list(est = res[[1]], CI = res[[2]],
                  #                                     est_old = est_old, CI_old = CI_old)
                  #                   class(returnVal) <- "ctmle"
                  #
                  #                   return(returnVal)
            }
            #-------------------------------------------------------------

            candidates$results.all$varIC_origscale <- candidates$results.all$varIC *
                  ifelse(maptoYstar, (diff(Qinit$ab))^2, 1)
            candidates$results.all$varDstar_origscale <- candidates$results.all$varDstar *
                  ifelse(maptoYstar, (diff(Qinit$ab))^2, 1)

            returnVal <- (list(candidates=candidates,
                               best_k=cv.res$best_k,
                               est= candidates$results.all$est[cv.res$best_k],
                               CI = candidates$results.all$est[cv.res$best_k]+c(-1.96,1.96) *
                                     sqrt(candidates$results.all$varIC_origscale[cv.res$best_k]/n),
                               Qinit=Qinit,
                               ab = Qinit$ab,
                               likelihood= candidates$results.all$likelihood[cv.res$best_k],
                               varIC = candidates$results.all$varIC_origscale[cv.res$best_k],
                               varDstar = candidates$results.all$varDstar_origscale[cv.res$best_k],
                               varIC.cv=cv.res$varIC[cv.res$best_k],
                               var.psi = candidates$results.all$varIC_origscale[cv.res$best_k]/n,
                               penlikelihood.cv = cv.res$penlikelihood[cv.res$best_k],
                               cv.res=cv.res))
            if(ctmletype == 1){
                  returnVal$est = res[[1]]
                  returnVal$CI = res[[2]]
                  returnVal$var.psi = res[[3]]
            }

            returnVal$pvalue <- 2*pnorm(-abs(returnVal$est)/sqrt(returnVal$var.psi))

            just_enough <- c(2:4,6,8,9,11,14)
            if (!detailed) {
                  returnVal <- returnVal[just_enough]
            }
      }

      class(returnVal) <- "ctmle"
      return(returnVal)
}


