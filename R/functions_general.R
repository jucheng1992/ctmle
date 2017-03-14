construct_candidates_general <- function(Y, X, Q,
                                         gn_candidates,
                                         ctmletype,
                                         family,Qbounds,training_set,
                                         gbound, like_type, verbose,stopFactor){
      g.dataset <- X
      Epsilon <- scorebws<-NULL
      # To keep track of when we need fluctuation
      best_gns <- NULL

      # Is the best score in all previous iterations
      minScore <- Inf

      gn_remaining <- gn_candidates
      n_candidates <- ncol(gn_candidates)

      i <- 0
      DONE <- (n_candidates<= 0)

      # browser()
      while(!DONE) {
            i <- i+1

            if(verbose){
                  cat("\tBeginning construction of clever covariate", i,"\n")
            }

            n_candidates <- ncol(gn_remaining)
            epsilon <- NULL


            score <- rep(NA, n_candidates)

            if(ctmletype == 1){
                  eps_try <- matrix(ncol = 2,nrow = n_candidates)
            }

            if(verbose){ cat("\n\t  Selecting best Nuisance Prameter Estimator to add to model...\n") }

            for (j in 1:n_candidates){
                  if(verbose) { cat("\t\tTrying", j, "th candidate estimator","\n") }

                  # Get the current g
                  g <- gn_remaining[training_set,j , drop=FALSE]


                  g1W.total <- bound(g, c(gbound, 1-gbound))
                  g0W.total <- 1 - bound(g, c(gbound, 1-gbound))
                  H1W <- X[training_set,1]/g1W.total
                  H0W <- (1 - X[training_set,1])/g0W.total
                  # browser()
                  if(ctmletype == 1){

                        suppressWarnings(eps_try[j,]<- coef(glm(Y[training_set] ~ -1 +
                                                                      offset(Q[training_set, "QAW"])
                                                                + H0W + H1W, family = family)))
                        Qstar <- Q[training_set,] + c((eps_try[j,1] * H0W +
                                                             eps_try[j,2] * H1W),
                                                      eps_try[j,1]/g0W.total,
                                                      eps_try[j,2]/g1W.total)
                  }

                  if(family=="binomial"){ Qstar <- plogis(Qstar)}

                  varIC <- calc_varIC(Y[training_set], Qstar, A=X[training_set,1],g1W=g1W.total)[1] * sum(training_set)/length(training_set)

                  # Compute current empirical loss
                  if(like_type == "RSS"){
                        score[j] <-  sum((Y[training_set] - Qstar[,"QAW"])^2) + varIC[1]
                  }else {
                        score[j] <- -sum(Y[training_set]*log(Qstar[,"QAW"]) + (1-Y[training_set])*log(1-Qstar[,"QAW"])) + varIC[1]
                  }
            }

            # browser()
            if(verbose) { cat("\t\t",paste("penalized", like_type,":"), round(score,5), "\n") }

            score[is.nan(score)] <- Inf

            # Find best score
            best <- which.min(abs(score))

            if(verbose) { cat("\t\tbest choice:", best, "th Estimator", "\n") }

            bestScore <- score[best]
            epsilon   <- rbind(epsilon, eps_try[best,])

            if(verbose) { cat("\t\tbest score:", bestScore, "\n") }

            earlyStop <- minScore*stopFactor < bestScore
            if(earlyStop & verbose){
                  cat("Stopping early because loss function of current candidate >", stopFactor, "times the best candidate seen so far\n")
                  cat("The ratio of best to current (",i,"), is ", round(bestScore/minScore, 2), "\n")
            }
            minScore <- ifelse(minScore<bestScore, minScore, bestScore)

            # Score the smaller the better
            # This means minScore > bestScore, which means at least one CTMLE get improved
            if(bestScore == minScore){
                  Epsilon  <- rbind(Epsilon, epsilon)
                  # Keep track the gn which mostly improve the empirical risk
                  # Then we need to fluctuate it before continue
                  best_gns  <- cbind(best_gns, gn_remaining[,best])
                  g <- gn_remaining[,best]

                  g1W.total <- bound(g, c(gbound, 1-gbound))
                  g0W.total <- 1 - bound(g, c(gbound, 1-gbound))
                  H1W <- X[,1]/g1W.total
                  H0W <- (1 - X[,1])/g0W.total

                  if(ctmletype==1){
                        Q[training_set,]  <- Q[training_set,] +
                              c(( Epsilon[dim(Epsilon)[1],1] * H0W[training_set] +
                                        Epsilon[dim(Epsilon)[1],2] * H1W[training_set]),
                                Epsilon[dim(Epsilon)[1],1]/g0W.total[training_set],
                                Epsilon[dim(Epsilon)[1],2]/g1W.total[training_set])
                        Q[training_set,]<-qlogis(bound(plogis(Q[training_set,]), Qbounds))
                  }

                  if( (best +1) > ncol(gn_remaining)){
                        gn_remaining <- integer(0)
                        n_candidates <-0
                  }     else{
                        gn_remaining <- gn_remaining[, (best + 1):ncol(gn_remaining), drop=FALSE]
                        n_candidates <- ncol(gn_remaining)
                  }

            }else{
                  n_candidates <- 0
            }
            DONE  <-  (n_candidates <= 0 | earlyStop)
            if(verbose){
                  if(bestScore==minScore){
                        cat("\tThe model for clever covariate", i, "is complete.\n")
                        cat("\t\t...Calculating h(A,W) based on candidate g-estimator", i,"\n")
                        cat("\t\t...Running a logistic regression to fit epsilon\n")
                        cat("\t\t...Updating estimate of Q(A,W) = Q(A,W) + epsilon * h(A,W)\n")
                  }
                  if(DONE){
                        cat("\tAll candidate TMLE estimators have been constructed\n\n")
                  } else {
                        cat("\n\tReady to use the updated estimate of Q(A,W) to construct the next clever covariate.\n")

                  }
            }
      }

      return(list(gns=best_gns, epsilon=Epsilon, earlyStop=earlyStop))
}


# Given a sequence of fluctuation (computed by construct_candidates_general function), this function returns C-TMLE estimates
evaluate_candidates_general <- function(Y,X,Q,
                                        gn_candidates,
                                        ctmletype, family,Qbounds,
                                        g.dataset,
                                        ab, candidates,
                                        like_type, gbound,
                                        training_set=1:length(Y)) {

      best_gns <- candidates$gns
      test_set <- !training_set | all(training_set)
      n_candidates <- ncol(gn_candidates)

      j <- 1
      nextbest_gn <- best_gns[,j]

      est <- likelihood <- varIC <- varDstar <-rep(Inf, n_candidates)

      if(ctmletype==1){ epsilon<-matrix(rep(Inf,n_candidates*2),ncol=2) }

      for (i in 1:n_candidates){

            g <- gn_candidates[,i]

            g1W.total <- bound(g, c(gbound, 1-gbound))
            g0W.total <- 1 - bound(g, c(gbound, 1-gbound))

            H1W <- X[,1]/g1W.total
            H0W <- (1 - X[,1])/g0W.total

            if(ctmletype==1){
                  suppressWarnings(epsilon[i,]<- coef(glm(Y ~ -1 + offset(Q[, "QAW"]) + H0W+ H1W, family = family)))
                  Qstar <- Q + c((epsilon[i,1] * H0W + epsilon[i,2] * H1W), epsilon[i,1]/g0W.total, epsilon[i,2]/g1W.total)
            }

            if(family=="binomial"){ Qstar <- plogis(Qstar) }

            est[i]   <- (mean(Qstar[training_set,"Q1W"]) - mean(Qstar[training_set,"Q0W"]))*diff(ab)

            # Compute the testing loss
            if(like_type=="RSS"){
                  likelihood[i]   <- sum((Y[test_set]-Qstar[test_set,"QAW"])^2)
            }else{
                  likelihood[i] <- -sum(Y[test_set]*log(Qstar[test_set,"QAW"]) + (1-Y[test_set])*log(1-Qstar[test_set,"QAW"]))
            }
            temp 	 <- calc_varIC(Y[test_set], Q=Qstar[test_set,],  A=X[test_set,1], W=X[test_set, 2], g1W=g1W.total[test_set], ICg=TRUE)
            varDstar[i] <- temp[1]
            varIC[i]    <- temp[2]
            if(is.nan(est[i])|is.infinite(est[i])){est[i] <- NA}
            if(is.nan(likelihood[i])|is.infinite(likelihood[i])){likelihood[i] <- NA}
            if(is.nan(varIC[i])|is.infinite(varIC[i])){varIC[i] <- NA}
            if(is.nan(varDstar[i])|is.infinite(varDstar[i])){varDstar[i] <- NA}


            if(  sum((gn_candidates[,i] - nextbest_gn) ^ 2)  < 1e-5){

                  # If current index is the one need fluctuation, than update Q
                  if(ctmletype==1){
                        Q <- Q +
                              c((epsilon[i,1] * H0W + epsilon[i,2] * H1W),
                                epsilon[i,1]/g0W.total, epsilon[i,2]/g1W.total)
                  }

                  Q <- qlogis(bound(plogis(Q), Qbounds))

                  # Update the gn for the next fluctuation
                  j <- j+1

                  if(ncol(best_gns) >= j){
                        nextbest_gn <- best_gns[, j]
                  }else{
                        # ???
                        nextbest_gn <- 0
                  }
            }
      }
      return(list(est=est, likelihood=likelihood, varDstar=varDstar, varIC=varIC))
}




# This is the core function for cv. This construct and evaluate ctmle estimator given training set index
stage2_general <- function(Y, X, Q,
                           gn_candidates,
                           ctmletype,
                           family, Qbounds,
                           ab,
                           training_set=rep(T,length(Y)),
                           like_type,  gbound,verbose,
                           stopFactor=10^6) {

      if(verbose) { cat ("\n\t\t-----Stage 2: Constructing candidate TMLE estimators-----\n\n") }


      candidates  <- construct_candidates_general(Y, X, Q, gn_candidates,ctmletype,
                                                  family=family[1] ,
                                                  Qbounds,
                                                  training_set=training_set,
                                                  gbound=gbound, like_type=like_type,
                                                  verbose=verbose,stopFactor=stopFactor)

      results.all <- evaluate_candidates_general(Y, X, Q, gn_candidates,ctmletype,
                                                 family[1],Qbounds, g.dataset=X,
                                                 candidates, ab=ab,
                                                 like_type=like_type, gbound=gbound,
                                                 training_set=training_set)

      return(list(candidates=candidates, results.all=results.all))
}


# This is core function of ctmle: v-fold CV to construct and select ctmle estimators.
cv_general <- function(Y,X, est.all, Q,
                       gn_candidates,
                       ctmletype, family, Qbounds, ab, like_type, gbound,
                       verbose=FALSE, PEN, V=5, folds = NULL) {
      # V <- 5

      # Adding user-specified folds optin
      n <- length(Y)
      nconstructed <- length(est.all)
      if(is.null(folds)){
            folds <- by(sample(1:n,n), rep(1:V, length=n), list)
      }else{
            if(length(folds) != V){
                  stop("The number of user-specified folds information does not match V")
            }else if(mean(sort(unlist(folds)) - 1:n)!= 0){
                  stop("Error in the indices of the user-specified folds")
            }
      }

      folds <- by(sample(1:n,n), rep(1:V, length=n), list)
      likelihood <- varDstar <- bias <- c(rep(0,nconstructed))
      est <- matrix(data=NA, nrow=nconstructed, ncol=V)

      if(verbose) {cat("\tfold: ")}
      for (v in 1:V){
            if(verbose) {cat(v," ")}
            test_set <- folds[[v]]
            training_set <- rep(TRUE, n)
            training_set[folds[[v]]] <- FALSE


            # candidates[[1]] is candidates, the fianl candidates
            # candidates[[2]] is results.all,  the information from CV
            candidates <- stage2_general(Y,X, Q, gn_candidates,ctmletype, family=family,
                                         Qbounds = Qbounds, ab=ab,
                                         training_set=training_set,
                                         like_type=like_type,
                                         gbound=gbound,
                                         verbose=FALSE)

            # Save estimators
            est[,v] <- candidates$results.all$est
            likelihood   <- c(likelihood+candidates$results.all$likelihood)
            varDstar   <- c(varDstar+candidates$results.all$varDstar)
            bias    <- c(bias+(candidates$results.all$est-est.all))

      }
      bias <- bias
      pen <- varDstar*(diff(ab))^2/V + n*(bias/V)^2
      if(identical(family[1], gaussian) | identical(family[1], "gaussian")){
            pen <- pen * log(n)
      }

      # If PEN, add penalization to negative log-likelihood
      if(PEN){
            score <- likelihood*(diff(ab)^2) + pen
      } else {
            score <- likelihood*(diff(ab)^2)
      }
      score[is.infinite(score)|is.nan(score)] <- Inf
      # Select by CV
      best_k <- which.min(abs(score))


      if(verbose){
            cat("\n\n")
            cat("\t terms: ", candidates$candidates$terms, "\n")
            cat("\t all estimates: ", est.all, "\n")
            cat("\t all RSS: ", score, "\n")
            cat("\t best_k: ", best_k, "\n\n")
      }
      return(list(best_k=best_k, est=est, likelihood=likelihood, like_type=like_type,
                  penlikelihood=score, varIC=varDstar/V, bias=bias/V, pen=pen))
}


