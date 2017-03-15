
#---------- function Tx_mech -------
# purpose: calculate h and g1W given a model for g and
#          treatment indicator, A.
# arguments:
#  g -  glm (or other object) with a prediction method
#       for P(A=1|W)
#  X -  design matrix, A in column 1, plus covariates W
#  gbound - bound for g
# returns:
#  h - an nx3 matrix.
#      h[,"hAW"] is used to fit epsilon, and update initial Q
#      h[,"h1W"] is used to update counterfactual Q(1,W)
#      h[,"h0W"] is used to update counterfactual Q(0,W)
#  g1W - P(A=1, W)
#----------------------------------------
Tx_mech <- function (g,X, gbound) {
      g1W <- bound(predict(g, newdata=X, type="response"), c(gbound, 1-gbound))
      h <- matrix(data=NA, nrow=length(g1W), ncol=3)
      colnames(h) <- c("hAW", "h1W", "h0W")
      h[,"hAW"]   <- h[,"h1W"] <- 1/g1W
      h[,"h0W"]   <- -1/(1-g1W)
      h[X$A==0, "hAW"] <- h[X$A==0, "h0W"]
      return(list(h=h, g1W=g1W))
}

#----------------function calc_varIC--------------
# purpose: calculate the variance of the IC,
# arguments:
# 	Y - outcome
#   Q - response, not linear predictors
#   h - clever covariate
#   A - target var
#   W - only the terms in g
#   g1W - p(A=1|W)
#   (Inference should take ICg term into account)
# returns: variance of the influence curve
#-------------------------------------------------
calc_varIC <- function(Y, Q, h = NULL, g1W=NULL, A=NULL, W=NULL, ICg=FALSE) {
      if(is.null(h)){
            Dstar <-  (A/g1W - (1 - A)/(1-g1W))  * (Y - Q[, "QAW"]) +
                  Q[, "Q1W"] - Q[, "Q0W"] -  (mean(Q[,"Q1W"]) - mean(Q[,"Q0W"]))
      }else{
            Dstar <- (Y-Q[,"QAW"]) * h + Q[,"Q1W"] - Q[,"Q0W"] - (mean(Q[,"Q1W"]) - mean(Q[,"Q0W"]))
      }
      varDstar <- var(Dstar)
      varIC <- varDstar
      if(ICg & !is.null(W) & !is.null(A)){
            if(NCOL(W) > 0){
                  W <- cbind(rep(1,nrow(as.matrix(W))), W)
                  term1 <-colMeans( -(Y-Q[,"QAW"]) *  W  *  (A * (1-g1W)/g1W + (1-A) * g1W/(1-g1W)))
                  term2 <- try(solve(matrix(rowMeans(apply(cbind(g1W, W), 1,
                                                           function(v){v[-1] %o% v[-1]  *  v[1]  *  (1-v[1])})),byrow=TRUE, nrow=ncol(W))),
                               silent=TRUE)
                  if(is.matrix(term2)){
                        IC <- as.vector(Dstar + term1 %*% term2 %*% t((A-g1W) * W))
                        varIC <- var(IC)
                  }
            } }
      return(c(varDstar, varIC))
}


#--------------function stage1----------------
# purpose:  Use Superlearner to get the stage 1 estimate of the density
# Note: SL library should include many more prediction algorithms in practice
# arguments:
#  Y - outcome var, only include obs where Y and A both observe
#  X - design matrix used to fit the models, A in first column
#  Qform - optional regression formula to predict Q
#  family
#  Qbounds - bound on initial Y and predicted values for Q.
#  Qinit - optional user-supplied values for $Q$
#  alpha: if we're doing the Y* transformation, keep Ystar values away from (0,1)
#  SL.library - only used if estimating Q with SL
#  maptoYstar - boolean - whether or not to map to Ystar and fluctuate on logistic scale
# returns:
#  Q - an nx3 matrix of predicted values:
#	 QAW - fitted Y corresponding to observed (A,W)
#	 Q1W - fitted Y as a result of intervening on node A by setting A=1
#    Q0W - fitted Y, A set to 0
#  Note: values in Q are on the logit scale for binary outcomes, Y
#--------------------------------------------------
stage1 <- function(Y,X, Qform=NULL, family, Qbounds,Qinit, alpha=.995, SL.library, cvQinit, maptoYstar=FALSE,verbose=FALSE){
      ab <- c(0,1)
      if(family == "binomial"){
            Qbounds <- sort(c(1-alpha, alpha))
      } else {
            if (is.null(Qbounds)){
                  Qbounds <- range(Y)
            }
            if(verbose){
                  cat("\tTruncating Y at (", paste(round(Qbounds,3), collapse=", "), ")\n\n")
            }
            Y <- bound(Y, Qbounds)
            if(!is.null(Qinit)){
                  Qinit <- bound(Qinit, Qbounds)
                  colnames(Qinit) <- c("QAW", "Q1W", "Q0W")
            }
      }
      if(maptoYstar) {
            ab <- range(Y)
            Y <- (Y-ab[1])/diff(ab)
            if(!is.null(Qinit)){
                  Qinit <- (Qinit-ab[1])/diff(ab)
            }
            Qbounds <- sort(c(1-alpha, alpha))
            if(verbose){
                  cat("\tMapping Y to Y* (a=", round(ab[1],3), ", b=", round(ab[2],3),") \n\n")
            }
      }
      if(is.null(Qinit)){
            if(verbose){
                  cat("\t\t-----Stage 1: Estimating the conditional mean of Y given A and W-----\n\n")
            }
            n <- length(Y)
            Q_SL <- NULL
            Q <- matrix(data=NA, nrow=n, ncol=3)
            colnames(Q) <- c("QAW", "Q1W", "Q0W")
            X1 <- data.frame(A=1, X[,-1])
            X0 <- data.frame(A=0, X[,-1])
            colnames(X0) <- colnames(X1) <- colnames(X)

            if(is.null(Qform)){

                  if(is.null(SL.library)){
                        SL.library <-  c("SL.glm", 'SL.glmnet')
                  }
                  if(cvQinit){
                        Q_SL <- try(estQ_cvSL.nomissing(Y, X, newX=rbind(X, X1, X0),
                                                        SL.library, family=family,
                                                        Qbounds=Qbounds, verbose=verbose))
                        if(class(Q_SL) != "try-error"){
                              Q <- Q_SL
                        }
                  } else {
                        if(verbose){ cat("\tStage 1: Estimating initial Q with Super Learner\n") }
                        Q_SL <- try(SuperLearner(Y=Y, X = X, newX=rbind(X, X1, X0),
                                                 SL.library=SL.library,
                                                 family=family))
                        if(identical(class(Q_SL), "SuperLearner")){
                              Q[1:(3*n)] <- Q_SL$SL.predict
                        }
                  }

                  if(identical(class(Q_SL), "try-error") | identical(class(Q_SL), "NULL")) {
                        if(verbose){
                              cat("\tSuper Learner prediction not available. Using main terms linear regression to estimate initial Q\n")
                        }
                        Qform <- paste ("Y~1+", paste(colnames(X), collapse="+"))
                        m <- glm(as.formula(Qform), data=data.frame(Y,X), family=family)
                        Q[,"QAW"] <- predict(m, newdata=data.frame(X), type="response")
                        Q[,"Q1W"] <- predict(m, newdata=X1, type="response")
                        Q[,"Q0W"] <- predict(m, newdata=X0, type="response")
                  }
            } else {
                  if(verbose){cat("\tEstimating initial Q using user-supplied regression formula\n")}
                  m <- glm(as.formula(Qform), data=data.frame(Y,X), family=family)
                  Q[,"QAW"] <- predict(m, newdata=data.frame(X),  type="response")
                  Q[,"Q1W"] <- predict(m, newdata=X1, type="response")
                  Q[,"Q0W"] <- predict(m, newdata=X0, type="response")
            }

      } else {
            Q <- bound(Qinit,Qbounds)
      }

      if (family=="binomial" | maptoYstar){
            Q <- qlogis(bound(Q, Qbounds))
            if(verbose){
                  cat("\n\tBounding predicted values for Q away from (0,1) at (", paste(round(Qbounds,3), collapse=", "),")\n\n")
            }
      } else {
            Q <- bound(Q, Qbounds)
      }

      family.stage2 <- c(ifelse(maptoYstar, "binomial", family), family)

      # project stage 1 estimate onto betaA + offset(QAW)
      #coefA <- suppressWarnings(coef(glm(Y~-1+X[,"A"]+offset(Q[,"QAW"]), family=family.stage2[1])))
      #coefA.mat <- cbind(coefA*X[,"A"], coefA, 0)

      #if(identical(family.stage2[1], binomial) | identical(family.stage2[1],"binomial")){
      #	Q <- qlogis(bound(plogis(Q+coefA.mat), c(alpha, 1-alpha)))
      #} else {
      #	Q <- bound(Q+coefA.mat, Qbounds)
      #}

      return(list(Q=Q, Y=Y, Qbounds=Qbounds, ab=ab, family=family.stage2))
}


