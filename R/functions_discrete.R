#--------------function select_terms----------
# purpose:  add terms to the current model for A
# 			as long as they continue to increase the
#			(penalized) likelihood
# arguments:
#  Y - outcome
#  X - design matrix, A in column 1
#  g.dataset - dataset containing A and W for all observations
#  Q: nx3 matrix for Q(A|W), Q(1|W), Q(0|W)
#  family - binomial for binary outcomes, gaussian for continuous outcomes
#  Qbounds - respect constraints on Y
#  Aform - starting formula for regression of A on W
#		   (e.g. terms forced in from previous iterations)
#  terms_remaining - terms to consider adding to the model for g
#  nterms - bound on number of terms to add
#  training_set - identifies which values of Tx$h correspond to
#                 the observations in Y and X (since g is fit on all obs).
#  gbound - bounds for predicted values of g
#  preOrder - added by Richie Wyss to allow for pre-ordering of variables
#--------------------------------------------------
select_terms <- function(Y,X,g.dataset=X, Q, family, Aform,
                         terms_remaining, nterms, training_set, covar, like_type,
                         gbound, verbose=TRUE, preOrder) {
      bestScore <- Inf
      improving <- TRUE
      DONE <- length(terms_remaining) == 0
      allInf <- FALSE
      counter <- 0
      epsilon <- NULL
      # first covar- special check for intercept-only model
      if(covar==1){
            g  <- glm(Aform,data=g.dataset, family=binomial)
            Tx <- Tx_mech(g, g.dataset[training_set,], gbound)
            epsilon <- suppressWarnings(coef(glm(Y~-1+offset(Q[,"QAW"])+Tx$h[,"hAW"], start = 0, family=family)))
            Qstar <- Q + epsilon  *  Tx$h
            if(family=="binomial"){
                  Qstar <- plogis(Qstar)
            }

            varIC <- calc_varIC(Y, Qstar, Tx$h[,"hAW"])[1] * sum(training_set)/length(training_set)
            if(like_type == "RSS"){
                  bestScore <-  sum((Y - Qstar[,"QAW"])^2) + varIC[1]
            } else {
                  bestScore <- -sum(Y*log(Qstar[,"QAW"]) + (1-Y)*log(1-Qstar[,"QAW"])) + varIC[1]
            }

            if(verbose){cat("\n\t\t",paste("penalized", like_type, "(intercept only):"), round(bestScore,5), "\n")}
      }

      while (!DONE) {
            counter <- counter +1
            score <- rep(NA, length(terms_remaining))
            eps_try <- rep(NA, length(terms_remaining))
            if(verbose){ cat("\n\t  Selecting best term to add to model...\n") }

            ###########################################
            ##	Addition to allow for pre-ordering
            if(preOrder==FALSE){
                  for (i in 1:length(terms_remaining)){
                        Aform_try <- paste(Aform, "+",terms_remaining[i])
                        if(verbose) { cat("\t\tTrying", Aform_try,"\n") }
                        g  <- glm(Aform_try,data=g.dataset, family=binomial)

                        Tx <- Tx_mech(g, g.dataset[training_set,], gbound)
                        eps_try[i] <- suppressWarnings(coef(glm(Y~-1+offset(Q[,"QAW"])+Tx$h[,"hAW"],start = 0, family=family)))
                        Qstar <- Q + eps_try[i]  *  Tx$h
                        if(family=="binomial"){
                              Qstar <- plogis(Qstar)
                        }
                        varIC <- calc_varIC(Y, Qstar, Tx$h[,"hAW"])[1] * sum(training_set)/length(training_set)
                        if(like_type == "RSS"){
                              score[i] <-  sum((Y - Qstar[,"QAW"])^2) + varIC[1]
                        } else {
                              score[i] <- -sum(Y*log(Qstar[,"QAW"]) + (1-Y)*log(1-Qstar[,"QAW"])) + varIC[1]
                        }
                  }
            }
            if(preOrder==TRUE){
                  Aform_try <- paste(Aform, "+",terms_remaining[1])
                  if(verbose) { cat("\t\tTrying", Aform_try,"\n") }
                  g  <- glm(Aform_try,data=g.dataset, family=binomial)
                  Tx <- Tx_mech(g, g.dataset[training_set,], gbound)
                  eps_try <- suppressWarnings(coef(glm(Y~-1+offset(Q[,"QAW"])+Tx$h[,"hAW"],start = 0, family=family)))
                  Qstar <- Q + eps_try  *  Tx$h
                  if(family=="binomial"){
                        Qstar <- plogis(Qstar)
                  }
                  varIC <- calc_varIC(Y, Qstar, Tx$h[,"hAW"])[1] * sum(training_set)/length(training_set)
                  if(like_type == "RSS"){
                        score <-  sum((Y - Qstar[,"QAW"])^2) + varIC[1]
                  } else {
                        score <- -sum(Y*log(Qstar[,"QAW"]) + (1-Y)*log(1-Qstar[,"QAW"])) + varIC[1]
                  }
            }
            if(verbose) {cat("\t\t",paste("penalized", like_type,":"), round(score,5), "\n")}
            score[is.nan(score)] <- Inf
            best <- which.min(abs(score))
            if(verbose) {cat("\t\tbest choice:", terms_remaining[best], "\n")}
            if (length(best)==0) {
                  allInf    <- TRUE
                  improving <- FALSE
                  # next lines assure an assignment later doesn't crash
                  best <- 1
                  score[best]   <- Inf
                  eps_try[best] <- 0
            }
            else if (is.infinite(score[best])) {
                  allInf    <- TRUE
                  improving <- FALSE
            } else {
                  improving <- (bestScore[length(bestScore)]>score[best] )
            }
            prevAform <- Aform
            Aform     <- paste(Aform,"+", terms_remaining[best])
            bestScore <- c(bestScore, score[best])
            epsilon   <- c(epsilon, eps_try[best])
            terms_remaining <- terms_remaining[-best]
            DONE      <-  !improving | length(terms_remaining) == 0 | nterms == counter

      }
      if(counter==1 & allInf) {
            if(length(terms_remaining)>1){
                  Aform <- paste(Aform," + ", paste(terms_remaining[1:(nterms-counter)], collapse=" + "), sep="")
            }
            epsilon <- c(epsilon, rep(0,(nterms-counter)))
            bestScore <- Inf
      } else if (!improving & (counter > 1 | allInf) ){
            if(verbose){cat("\t\t(However, no additional term improves the penalized", like_type, ", so none of these\n\t\twill be included in the model for the current clever covariate.)\n\n")}
            Aform     <- prevAform
            epsilon   <- epsilon[-length(epsilon)]
            bestScore <- bestScore[length(bestScore)-1]
      } else {
            bestScore <- bestScore[length(bestScore)]
      }

      return(list(Aform=Aform, epsilon=epsilon, score=bestScore))
}
#-------------function evaluate_candidates-------
# purpose:  calculate RSS, var, estimate for each candidate estimator
#           adds Gcomp on initial as 1st estimator, so k=2 corresponds
#			to tmle estimator using first g-hat to calculcate h.
# arguments:
#   Y - outcome var
#   X - design matrix  (when called from cv, this should be test set observations)
#   Q - nx2 matrix of values corresponding to initial Q[QAW, Q1W]
#   family - binomial for binary outcomes, gaussian for continuous outcomes
#   Qbounds - respect constraints on Y
# 	g.dataset - typically all observations, used to fit g
#   candidates:  a list with one entry for each candidate estimator
# 	  term - term name, in order of incorporation into clever covariates
#	  covar - a number indicating which clever covariate the term is a member of
#     epsilon - epsilon corresponding to each clever covar
#  ncandidates - number of candidates to evaluate
#  like_type "RSS" or "loglike"
#  gbound - bound on predicted value for g
#  training_set - obs for model fits
#  returns: likelihood, var(Dstar+IC_g), estimate for each candidate
#  note: each term "within" the same clever covar has the same offset,
#		 but its epsilon will be different
#-------------------------------------------------
evaluate_candidates <- function(Y,X,Q, family, g.dataset, ab, candidates, ncandidates, like_type, gbound,training_set=1:length(Y)) {

      #---------------function which.cols----------------
      # given a vector of column names and a formula
      # return indices of the names of terms on rhs of formula
      #--------------------------------------------------
      which.cols <- function(cnames, form){
            fterms <- attr(terms(as.formula(form)),"term.labels")
            which(cnames %in% fterms)
      }
      Qinit <- Q  # don't need this?
      test_set <- !training_set | all(training_set)

      Aform = "A~1"
      nterms <- sum(!is.na(candidates$terms))
      cum_nterms <- cumsum(table(candidates$covar))				##Richie comment: cum_nterms contains as many elements as clever covariates. e.g., c(3,5,6) means first clever covariate has 3 variables (including intercept) next clever covariate has 5 variables and final clever covariate has 6 variables (including intercept).

      ncovar <- length(candidates$epsilon)
      est <- likelihood <- varIC <- varDstar <- rep(Inf, ncandidates+1)
      covar <- 1
      for (i in 1:nterms){
            Aform <- paste(Aform, "+",candidates$terms[i], sep="")

            g <- glm(Aform, data=g.dataset, family=binomial) #fit g on all obs
            Tx <- Tx_mech(g, g.dataset, gbound)
            Qstar <- Q + candidates$epsilon[i]  *  Tx$h
            if(family=="binomial"){
                  Qstar <- plogis(Qstar)
            }
            # track est on training_set to use for bias calculation in cv function

            est[i]   <- (mean(Qstar[training_set,"Q1W"]) - mean(Qstar[training_set,"Q0W"]))*diff(ab)

            if(like_type=="RSS") {
                  likelihood[i]   <- sum((Y[test_set]-Qstar[test_set,"QAW"])^2)
            } else {
                  likelihood[i] <- -sum(Y[test_set]*log(Qstar[test_set,"QAW"]) + (1-Y[test_set])*log(1-Qstar[test_set,"QAW"]))
            }
            temp 	 <- calc_varIC(Y[test_set], Q=Qstar[test_set,], h=Tx$h[test_set,"hAW"], A=X[test_set,1],
                                 W=X[test_set, which.cols(colnames(X), Aform)], g1W=Tx$g1W[test_set], ICg=TRUE)
            varDstar[i] <- temp[1]
            varIC[i]    <- temp[2]
            if(is.nan(est[i])|is.infinite(est[i])){est[i] <- NA}
            if(is.nan(likelihood[i])|is.infinite(likelihood[i])){likelihood[i] <- NA}
            if(is.nan(varIC[i])|is.infinite(varIC[i])){varIC[i] <- NA}
            if(is.nan(varDstar[i])|is.infinite(varDstar[i])){varDstar[i] <- NA}

            #jump to next clever covariate at pre-determined model sizes
            if (i == cum_nterms[covar]) {
                  covar <- covar+1
                  Q <- Q + candidates$epsilon[i]  *  Tx$h
            }

      }

      # this estimate is on the original scale, but varIC and varDstar are not
      return(list(est=est, likelihood=likelihood, varDstar=varDstar, varIC=varIC))
}

#--------------function construct_candidates------
# purpose - build increasingly non-p candidate g-hats
## construct clever covariates using forward selection
# for each clever covar add 1 or more variables to terms in the previous clever covar,
# until no improvement in RSS.  Iterate until all terms have been added
# For stability, calculate g based on data from full sample.
# arguments:
#  Y - outcome vector
#  X - design matrix (A in first column)
#  Q - stage 1 estimate of Q
#  family - binomial for binary outcomes, gaussian for continuous outcomes
#  ncandidates - max # candidates to construct - i think ignored for now.
#  terms - names of covariates to consider adding to model for g
#  training_set - indicates which obs in Tx$h correspond to obs in Y, X
#  gbound - bounds for predicted values for g1W
#  like_type - RSS or loglike
# return:
#  terms - a vector containing the covariates in the order that they were added
#		to the model for g,
#  covar - which clever covariate each term is in,
#  epsilon -  corresponding epsilon values.
#  preOrder - added by Richie Wyss to allow for pre-ordering
# this function modified for early stopping
#--------------------------------------------------
construct_candidates <- function(Y, X, Q, family, ncandidates=ncol(X)-1,terms,training_set, gbound, like_type, verbose=FALSE, preOrder, stopFactor=10^6){

      g.dataset <- X
      terms_remaining <- terms
      nterms <- length(terms_remaining)
      Aform <- "A~1"
      epsilon <- NULL
      covar <- NULL
      minScore <- Inf
      earlyStop <- FALSE
      i <- 0
      DONE <- nterms <= 0

      while(!DONE) {
            i <- i+1
            if(verbose) {cat("\tBeginning construction of clever covariate", i,"\n")}

            nextCandidates <- select_terms(Y[training_set],X[training_set,],
                                           g.dataset, Q[training_set,], family,  Aform, terms_remaining,
                                           nterms=length(terms_remaining), training_set, covar=i, like_type=like_type,
                                           gbound=gbound, verbose=verbose, preOrder=preOrder)
            earlyStop <- minScore*stopFactor < nextCandidates$score

            #browser()
            if(earlyStop & verbose){
                  cat("Stopping early because loss function of current candidate >", stopFactor, "times the best candidate seen so far\n")
                  cat("The ratio of best to current (",i,"), is ", round(nextCandidates$score/minScore, 2), "\n")
            }
            minScore <- ifelse(minScore<nextCandidates$score, minScore, nextCandidates$score)

            newAform <- nextCandidates$Aform							  ##Richie comment: newAform is the form of latest clever covariate

            epsilon  <- c(epsilon, nextCandidates$epsilon)					  ##epsilon contains multiple values. Each variable and intercept has an epsilon. nextCandidates$epsilon also has multiple epsilons for each variable selected in that round for that particular clever covariate.

            newterms <- setdiff(attr(terms(as.formula(newAform)),"term.labels"),
                                attr(terms(as.formula(Aform)),"term.labels"))			  ##newterms are the new variables selected for the next clever covariate
            if(i==1) {covar <- 1}
            terms_remaining <- terms_remaining[-which(terms_remaining %in% newterms)] ##terms_remaining is variables that haven't been selected yet
            covar <- c(covar, rep(i, length(newterms)))					  ##covar is a vector of numbers e.g., c(1,1,1,2,2,2,etc) where each number corresponds to the clever covariate the variable was first selected. Includes the intercept with a value of 1.
            g  <- glm(newAform, data=g.dataset, family="binomial")
            Tx <- Tx_mech(g, g.dataset[training_set,], gbound)
            Q[training_set,]  <- Q[training_set,] + epsilon[length(epsilon)] * Tx$h

            Aform <- newAform
            DONE  <- length(terms_remaining)==0 | ncandidates <= length(covar)-1 | earlyStop

            if(verbose){
                  cat("\tThe model for clever covariate", i, "is complete.\n")
                  cat("\tConstructed regression equation for estimating g(A,W)=p(A=1|W):\n\t\t\t", Aform, "\n\n")
                  cat("\t\t...Calculating h(A,W) based on candidate g-estimator", i,"\n")
                  cat("\t\t...Running a logistic regression to fit epsilon\n")
                  cat("\t\t...Updating estimate of Q(A,W) = Q(A,W) + epsilon * h(A,W)\n")
                  if(DONE){
                        cat("\tAll candidate TMLE estimators have been constructed\n\n")
                  } else {
                        cat("\n\tReady to use the updated estimate of Q(A,W) to construct the next clever covariate.\n")
                        cat("\t(The base model for clever covariate", i+1, "contains all terms included in clever covariate", i,")\n\n")
                  }
            }
      }
      terms=c(attr(terms(as.formula(Aform)), "term.labels")) 			##Richie comment: terms is all variables in order as selected for final clever covariate. Note: 1 is added for intercept in return object below.

      end <- min(length(terms), ncandidates)+1			 			##Richie comment: end is number of variables plus 1
      return(list(terms=c('1',terms[1:(end-1)]), covar=covar[1:end], epsilon=epsilon[1:end], earlyStop=earlyStop))
}

#----- function estQ_cvSL ----
# purpose: Obtain cross-validated estimates for initial Q using discrete SL.
# 	This function will return cross-validated predicted initial values
#   corresponding to the best algorithm in the library, or the convex combination
#   calculated by SL itself, as chosen by cvRSS.
#   The fitted value for each observation i, is predicted from a fit based on a training set
#	excluding the fold containing observation i.
# arguments:
# 	Y - outcome
# 	X - design matrix with A in first column
# 	SL.library - prediction algorithms
# 	V - number of outer cv-folds for disscrete SL
# 	V_SL - number of folds for internal SL cross-validation
# 	family - binomial or gaussian
# returns:
#  Q - nx3 matrix of predicted values on linear scale (e.g. logit if Qfamily=binomial)
#-----------------------------------------------
estQ_cvSL.nomissing <- function(Y,X,newX=X,SL.library, V=5, V_SL=5, family="gaussian", Qbounds, verbose){
      Q <- cvRSS <- NULL
      n <- length(Y)
      fold <- by(sample(1:n),rep(1:V, length.out=n),function(x){x})
      n_predictors <-length(SL.library)

      if(verbose){ cat("\tStage 1: Estimating cross-validated initial Q with Super Learner\n") }

      # We'll create a matrix of predictions - one column for each predictor in the library
      # plus one more for SL itself, with 3*n predicted values per column, corresponding to newX.
      predictions <- matrix(nrow=3*n, ncol=length(SL.library)+1)
      m_SL <- NULL
      for (v in 1:V) {
            fold_rows <- c(fold[[v]], fold[[v]]+n, fold[[v]]+2*n)
            if(class(m_SL) != "try-error"){
                  m_SL <- try(SuperLearner(Y=Y[-fold[[v]]], X=X[-fold[[v]],], newX=newX[fold_rows,],
                                           # !!! maybe come from the changing of SL package
                                           control = list(V=V_SL), family=family,SL.library=SL.library))
            }
            if(class(m_SL) != "try-error"){
                  predictions[fold_rows,1] <- m_SL$SL.predict
                  for (s in 1:n_predictors){
                        predictions[fold_rows,s+1] <- m_SL$library.predict[,s]
                  }
                  predictions <- bound(predictions, Qbounds)
            }
      }
      cvRSS <- colSums((Y-predictions[1:n,])^2)
      best <- which.min(cvRSS)
      best_alg <- c("SL", SL.library)[best]
      Q <- matrix(data=predictions[,best], nrow=n, ncol=3, byrow=FALSE)
      colnames(Q) <- c("QAW", "Q1W", "Q0W")
      if(verbose){cat("\tDiscrete SL: best algorithm = ", best_alg,"\n")}
      if (is.null(Q) | class(m_SL) == "try-error"){
            Q <- 0
            class(Q) <- "try-error"
      }
      return(Q)
}


#--------------function stage2-----------------
# purpose:  construct candidate cTMLE estimators
# based on Stage 1 fit and candidate g-hat nuisance parameter
# estimators constructed within this function
# arguments:
#  Y - outcome vector (training set passed in from cv)
#  X - design matrix (training set for cv)
#  Q - predicted values on training set
#  family - binomial for binary outcomes, gaussian for continuous outcomes
# 			a vector of length 2 here, though.  The first is the family to use for
#		 	fluctuating the initial estimate of Q.  The second is the original family
#			setting. These are only different when the inputs are mapped to Ystar
#			In this case fluctuations are carried out on the logistic scale, but
#			the initial projection onto A + offset(QAW) is carried out on the gaussian.
#			One could argue that this projection operation properly belongs in stage1,
#			but then one would be unable to cross-validate the projection, so one would
#		    be correct, but worse off.
#  Qbounds - respect constraints on Y after projecting on A + QAW, but do not enforce
# 			after targeting with TMLE, as this would not then solve the eff IC equation.
#  ncandidates - max number of candidates to construct
#  training_set
#  like_type = RSS or loglike
#  enrich - boolean, whether or not to enrich covariate set
#  gbound - bounds on predicted values for g1W
#  preOrder - addition by Richie Wyss to allow for pre-ordering of variables.
# returns:
#  candidates - candidates constructed by forward selection
#  results.all- values of RSS, varIC, estimate, for each candidate, including unadj.
# this function modified for early stopping
#----------------------------------------------
stage2 <- function(Y, X, Q, family, Qbounds, ab, ncandidates,
                   training_set=rep(T,length(Y)), like_type,  gbound,verbose=FALSE, preOrder,
                   stopFactor=10^6) {

      if(verbose) { cat ("\n\t\t-----Stage 2: Constructing TMLE estimators on full data-----\n\n") }

      candidates  <- construct_candidates(Y, X, Q, family[1],
                                          ncandidates=ncandidates, terms=colnames(X)[-1], training_set=training_set,
                                          gbound=gbound, like_type=like_type, verbose=verbose, preOrder=preOrder,
                                          stopFactor=stopFactor)
      results.all <- evaluate_candidates(Y, X, Q, family[1], g.dataset=X, candidates, ab=ab, ncandidates=ncandidates,
                                         like_type=like_type, gbound=gbound, training_set=training_set)
      return(list(candidates=candidates, results.all=results.all, terms=candidates$terms))
}

#--------------function cv-----------------
# purpose:  cross-validate stage 2 of the C-TMLE procedure
# split data into folds and track est, RSS, varIC, bias
# for each fold
# arguments:
#  Y - outcome
#  X - design matrix, A in column 1
#  est.all - vector of estimates corresponding to candidates fitted on full model
#  Q - nx3 matrix of initial fitted values, Q0[QAW, Q1W, Q0W]
#  family - binomial for binary outcomes, gaussian for continuous outcomes
#  Qbounds - respect constraints on Y
#  ncandidates - max # candidate nuisance parameter estimators to construct
#  like_type - "RSS" or "log_like", for negloglikelihood
#  gbound - truncation level for g
# returns:  est, cv-likelihood, cv-varIC, cv-bias for each estimator, and selected candidate
# added an argument for number of folds - V=5 by default
#---------------------------------------------
# changed from varIC to varDstar, 12/7/11

################################################################################################################################################
# modified by Richie Wyss to allow investigator to specify patience parameter which will stop cv early if the score does not improve after a
# certain number of covariates. Also modified to allow for pre-ordering of variables.

cv <- function(Y,X, Q, family, Qbounds, ab, ncandidates, like_type,
               gbound, verbose=FALSE, PEN, V, stopFactor=10^6,
               terms=colnames(X)[-1], patience, preOrder, folds = NULL) {

      n <- length(Y)
      nconstructed <- dim(X)[2]-1
      totalEst<- dim(X)[2]

      #set.seed(1001)
      if(is.null(folds)){
            folds <- by(sample(1:n,n), rep(1:V, length=n), list)
      }else{
            if(length(folds) != V){
                  stop("The number of user-specified folds information does not match V")
            }else if(mean(sort(unlist(folds)) - 1:n)!= 0){
                  stop("Error in the indices of the user-specified folds")
            }
      }

      likelihood <- varDstar <- bias <- list() #c(rep(0,nconstructed+1))
      est <- matrix(data=NA, nrow=nconstructed+1, ncol=V)
      est.temp <- likelihood.temp <- varIC.temp <- varDstar.temp <- ncandidates.temp<- list()
      Q<- list(Q,Q,Q,Q,Q)
      Aform.evaluate<- list()
      g.dataset <- X
      terms_remaining <- list(terms,terms,terms,terms,terms)
      nterms1 <- c(1,1,1,1,1)
      Aform <- list("A~1","A~1","A~1","A~1","A~1")
      epsilon <- list(NULL,NULL,NULL,NULL,NULL)
      covar <- list(NULL,NULL,NULL,NULL,NULL)
      minScore <- list(Inf,Inf,Inf,Inf,Inf)
      earlyStop <- FALSE
      i <- 0
      patience.count<- 0
      DONE <- FALSE

      while(!DONE) {

            i <- i+1
            likelihood[[i]] <- varDstar[[i]] <- bias[[i]] <- rep(0,totalEst)

            if(verbose) {cat("\tBeginning construction and cross validation of clever covariate", i,"for cross validation \n")}

            for (v in 1:V){
                  #if(verbose & nterms1[v]<=totalEst) {cat("\tfold: ", v, "\n")}
                  #if(verbose & nterms1[v]>totalEst) {cat("\tfold: ", v, "is done \n")}

                  test_set <- folds[[v]]
                  training_set <- rep(TRUE, n)
                  training_set[folds[[v]]] <- FALSE

                  ##comment: nterms1[v] keeps track of what variable we are on within each fold
                  if(nterms1[v]<=totalEst){

                        nextCandidates <- select_terms(Y[training_set],X[training_set,],
                                                       g.dataset, Q[[v]][training_set,], family[1],  Aform[[v]], terms_remaining[[v]],
                                                       nterms=length(terms_remaining[[v]]), training_set, covar=i, like_type=like_type,
                                                       gbound=gbound, verbose=FALSE, preOrder=preOrder)
                        earlyStop <- minScore[[v]]*stopFactor < nextCandidates$score
                        if(earlyStop & verbose){
                              cat("Stopping cross-validation early because loss function of current candidate >", stopFactor, "times the best candidate seen so far\n")
                              cat("The ratio of best to current (",i,"), is ", round(nextCandidates$score/minScore[[v]], 2), "\n")
                        }
                        minScore[[v]] <- ifelse(minScore[[v]]<nextCandidates$score, minScore[[v]], nextCandidates$score)

                        newAform <- nextCandidates$Aform							  ##Richie comment: newAform is the form for new, or latest clever covariate
                        epsilon[[v]]  <- c(epsilon[[v]], nextCandidates$epsilon)			  ##epsilon contains multiple values. Each variable and intercept has an epsilon. nextCandidates$epsilon also has multiple epsilons for each variable selected in that round for that particular clever covariate.
                        newterms <- setdiff(attr(terms(as.formula(newAform)),"term.labels"), attr(terms(as.formula(Aform[[v]])),"term.labels"))    ##newterms are the new variables selected for the next clever covariate
                        if(i==1) {covar[[v]] <- 1}
                        terms_remaining[[v]] <- terms_remaining[[v]][-which(terms_remaining[[v]] %in% newterms)] ##terms_remaining is variables that haven't been selected yet
                        covar[[v]] <- c(covar[[v]], rep(i, length(newterms)))					  ##covar is a vector of numbers e.g., c(1,1,1,2,2,2,etc) where each number corresponds to the clever covariate the variable was first selected. Includes the intercept with a value of 1.
                        g  <- glm(newAform, data=g.dataset, family="binomial")
                        Tx <- Tx_mech(g, g.dataset[training_set,], gbound)
                        Aform[[v]] <- newAform

                        ##note: when testing remove pounds to see what is going on with comments below
                        #if(verbose){
                        #	cat("\tThe model for clever covariate", i, "for fold", v, "is complete.\n")
                        #	cat("\tConstructed regression equation for cross validation estimating g(A,W)=p(A=1|W):\n\t\t\t", Aform[[v]], "\n\n")
                        #	cat("\t\t...Calculating h(A,W) for cross validation based on candidate g-estimator", i,"\n")
                        #	cat("\t\t...Running a logistic regression to fit epsilon for cross validation \n")
                        #	if(DONE){
                        #		cat("\tAll candidate TMLE estimators have been constructed for each fold in cross validation\n\n")
                        #	} else {
                        #		cat("\n\tReady to use the updated estimate of Q(A,W) to construct the next clever covariate for cross validation.\n")
                        #		cat("\t(The base model for clever covariate", i+1, " in cross validation contains all terms included in clever covariate", i,")\n\n")
                        #	}
                        #}

                        #---------------function which.cols----------------
                        # given a vector of column names and a formula
                        # return indices of the names of terms on rhs of formula
                        #--------------------------------------------------
                        which.cols <- function(cnames, form){
                              fterms <- attr(terms(as.formula(form)),"term.labels")
                              which(cnames %in% fterms)
                        }
                        test_set <- !training_set | all(training_set)
                        terms=c(1, c(attr(terms(as.formula(Aform[[v]])), "term.labels")))
                        nterms <- length(terms)

                        ncovar <- length(epsilon[[v]])

                        if(i==1){
                              ncandidates.temp[[v]] <- length(terms)
                              est.temp[[v]] <- likelihood.temp[[v]] <- varIC.temp[[v]] <- varDstar.temp[[v]] <- rep(Inf, ncandidates.temp[[v]])
                              Aform.evaluate[[v]] = "A~1"
                              nterms1[v] <- 1
                        }

                        if(i>1 & length(c(est.temp[[v]], rep(Inf, length(newterms)))) <= totalEst){
                              est.temp[[v]] <- c(est.temp[[v]], rep(Inf, length(newterms)))
                              likelihood.temp[[v]] <- c(likelihood.temp[[v]], rep(Inf, length(newterms)))
                              varIC.temp[[v]] <- c(varIC.temp[[v]], rep(Inf, length(newterms)))
                              varDstar.temp[[v]] <- c(varDstar.temp[[v]], rep(Inf, length(newterms)))
                        }

                        for (j in nterms1[v]:nterms){

                              Aform.evaluate[[v]] <- paste(Aform.evaluate[[v]], "+", terms[j], sep="")
                              g <- glm(Aform.evaluate[[v]], data=g.dataset, family=binomial) #fit g on all obs

                              Tx <- Tx_mech(g, g.dataset, gbound)
                              #Qstar <- Q + candidates$epsilon[j]  *  Tx$h
                              Qstar <- Q[[v]] + epsilon[[v]][j]  *  Tx$h

                              if(family[1]=="binomial"){
                                    Qstar <- plogis(Qstar)
                              }

                              est.temp[[v]][j]   <- (mean(Qstar[training_set,"Q1W"]) - mean(Qstar[training_set,"Q0W"]))*diff(ab)

                              if(like_type=="RSS") {
                                    likelihood.temp[[v]][j]   <- sum((Y[test_set]-Qstar[test_set,"QAW"])^2)
                              } else {
                                    likelihood.temp[[v]][j] <- -sum(Y[test_set]*log(Qstar[test_set,"QAW"]) + (1-Y[test_set])*log(1-Qstar[test_set,"QAW"]))
                              }

                              temp 	 <- calc_varIC(Y[test_set], Q=Qstar[test_set,], h=Tx$h[test_set,"hAW"], A=X[test_set,1],
                                                   W=X[test_set, which.cols(colnames(X), Aform[[v]])], g1W=Tx$g1W[test_set], ICg=TRUE)

                              varDstar.temp[[v]][j] <- temp[1]
                              varIC.temp[[v]][j]    <- temp[2]
                              if(is.nan(est.temp[[v]][j])|is.infinite(est.temp[[v]][j])){est.temp[[v]][j] <- NA}
                              if(is.nan(likelihood.temp[[v]][j])|is.infinite(likelihood.temp[[v]][j])){likelihood.temp[[v]][j] <- NA}
                              if(is.nan(varIC.temp[[v]][j])|is.infinite(varIC.temp[[v]][j])){varIC.temp[[v]][j] <- NA}
                              if(is.nan(varDstar.temp[[v]][j])|is.infinite(varDstar.temp[[v]][j])){varDstar.temp[[v]][j] <- NA}

                        } ##Richie comment: ends for (j in nterms1:nterms) loop

                  } ##ends if(nterms1[v]<nterms) statement


                  end <- length(est.temp[[v]])

                  est[1:end,v] <- est.temp[[v]]  ##Richie comment: est is a matrix with nrows= total # of estimates and ncol=#cross V
                  likelihood[[i]]   <- c(likelihood[[i]][1:end]+likelihood.temp[[v]], rep(Inf, nconstructed[1]+1-end))
                  varDstar[[i]]   <- c(varDstar[[i]][1:end]+varDstar.temp[[v]], rep(Inf, nconstructed[1]+1-end))
                  #bias[[i]]   <- c(bias[[i]][1:end]+(est.temp[[v]]-est.all[1:end]),rep(Inf, nconstructed[1]+1-end))  #need to comment this out since running cv before calculating actual estimates (est.all)

                  ##need to update Q after each clever covariate
                  Q[[v]]  <- Q[[v]] + epsilon[[v]][length(epsilon[[v]])] * Tx$h

                  ##updating nterms1
                  if(nterms1[v]<=nterms) nterms1[v]<- nterms+1

                  if(verbose & nterms1[v]<=totalEst) {cat("\tfold: ", v, "\n")}
                  if(verbose & nterms1[v]>totalEst) {cat("\tfold: ", v, "is done \n")}
            } ##Richie comment: end of CV loop


            Aform.evaluate<- Aform

            # multiple clever covars used in upate
            pen <- varDstar[[i]]*(diff(ab))^2/V  #+ n*(bias[[i]]/V)^2	#cannot calculate bias when running cv first (see comment above)
            if(identical(family[1], gaussian) | identical(family[1], "gaussian")){
                  pen <- pen * log(n)
            }
            if(PEN){
                  score <- likelihood[[i]]*(diff(ab)^2) + pen
            } else {
                  score <- likelihood[[i]]*(diff(ab)^2)
            }
            score[is.infinite(score)|is.nan(score)] <- Inf


            DONE<- max(likelihood[[i]]) < Inf


            ###checking patience parameter to stop early
            if(patience){
                  pterms<- sum(!is.infinite(score))
                  if(i==1){
                        pointer<- score[1]
                        pterms1<- 1
                  }
                  for(j in pterms1:pterms){
                        if(score[j]>pointer){
                              patience.count<- patience.count+1
                        }

                        if(patience.count>=patience) DONE<- TRUE
                        if(score[j]<=pointer & DONE==FALSE){
                              pointer<- score[j]
                              patience.count<- 0
                        }
                        if(j==pterms) DONE2<- TRUE
                  } 						##ends for(j in pterms1:pterms) loop
                  pterms1<- pterms+1
            } 							##ends if(patience) conditional

      } ##Richie comment: closees while(!DONE) in construct_candidates function

      best_k<- max(which(score[!is.infinite(score)] == min(score[!is.infinite(score)], na.rm=TRUE)))

      if(verbose){
            cat("\n\n")
            cat("\t terms evaluated by cross validation: ", terms, "\n")
            cat("\t best_k selected by cross validation: ", best_k, "\n\n")
      }
      return(list(best_k=best_k, likelihood=likelihood[[i]], like_type=like_type, penlikelihood=score, varIC=varDstar[[i]]/V, bias=bias[[i]]/V, pen=pen))

}

