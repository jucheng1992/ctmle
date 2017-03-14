
# Given a sequence of fluctuation, this function returns C-TMLE estimates
evaluate_candidates_glmnet <- function(Y,X,Q,lambdas, ctmletype,
                                family,Qbounds, g.dataset, ab, candidates,
                                like_type, gbound,
                                training_set=1:length(Y),
                                best_k = NULL) {

      Qstar_best <- 0
      bestlambdas<-candidates$lambdas
      test_set <- !training_set | all(training_set)
      ncandidates <- length(lambdas)
      # Qstar_all <- list()
      j<-1
      nextbestlambdas<-bestlambdas[j]
      est <- likelihood <- varIC <- varDstar <-rep(Inf, ncandidates)
      if(ctmletype==1){ epsilon<-matrix(rep(Inf,ncandidates*2),ncol=2) }
      if(ctmletype==2){ epsilon<-matrix(rep(Inf,ncandidates*4),ncol=4) }
      for (i in 1:ncandidates){
            f1 <- eval(paste("A ~ ", paste(paste(names(X[,-1]), sep=''), collapse=" + ")))
            f2 <- as.formula(f1)
            lassox<- model.matrix(f2, X)[,-1]
            lassoy<-as.factor(X[,1])
            gbw <- glmnet(x=lassox,y=lassoy,family="binomial")
            g <- predict(gbw,newx=lassox,s=lambdas[i],type="response")
            g1W.total <- bound(g, c(gbound, 1-gbound))
            g0W.total <- 1 - bound(g, c(gbound, 1-gbound))
            H1W <- X[,1]/g1W.total
            H0W <- (1 - X[,1])/g0W.total

            if(ctmletype==1){
                  suppressWarnings(epsilon[i,]<- coef(glm(Y ~ -1 + offset(Q[, "QAW"]) + H0W+ H1W, family = family)))
                  Qstar <- Q + c((epsilon[i,1] * H0W + epsilon[i,2] * H1W), epsilon[i,1]/g0W.total, epsilon[i,2]/g1W.total)
            }
            if(ctmletype==2){
                  if(i==1){
                        g1<-predict(gbw,newx=lassox,s=(lambdas[i]+0.005),type="response")
                  }else{
                        g1<-predict(gbw,newx=lassox,s=(lambdas[i-1]),type="response")
                  }
                  g1W.total1 <- bound(g1, c(gbound, 1-gbound))
                  g0W.total1 <- 1 - bound(g1, c(gbound, 1-gbound))
                  ddg1<--(g1W.total-g1W.total1)
                  ddg1[which(ddg1==0)]<-1e-10
                  ddg0<--(g0W.total-g0W.total1)
                  ddg0[which(ddg0==0)]<-1e-10
                  ddH1W <- (X[,1]/(g1W.total^2))*ddg1
                  ddH0W <- ((1 - X[,1])/(g0W.total^2))*ddg0
                  suppressWarnings(epsilon[i,]<- coef(glm(Y ~ -1 +offset(Q[, "QAW"])+ H0W+H1W+ ddH0W + ddH1W, family = family)))
                  Qstar <- Q + cbind((epsilon[i,1] * H0W+ epsilon[i,2] * H1W),
                                     epsilon[i,1]/g0W.total, epsilon[i,2]/g1W.total)+ cbind((epsilon[i,3] * ddH0W+ epsilon[i,4] * ddH1W),
                                                                                            (epsilon[i,3]/(g0W.total^2))*ddg0, (epsilon[i,4]/(g1W.total^2))*ddg1)
            }




            if(family=="binomial"){ Qstar <- plogis(Qstar) }

            est[i]   <- (mean(Qstar[training_set,"Q1W"]) - mean(Qstar[training_set,"Q0W"]))*diff(ab)
            # Qstar_all[[i]] <- cbind(Qstar_all, Qstar)

            if(!is.null(best_k)){
                  if(i == best_k){
                        Qstar_best <- Qstar
                  }
            }


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


            if(lambdas[i]==nextbestlambdas){
                  if(ctmletype==1){
                        Q <- Q +  c((epsilon[i,1] * H0W + epsilon[i,2] * H1W),
                                    epsilon[i,1]/g0W.total, epsilon[i,2]/g1W.total)
                  }
                  if(ctmletype==2){
                        Q <- Q + cbind((epsilon[i,1] * H0W+ epsilon[i,2] * H1W),
                                       epsilon[i,1]/g0W.total, epsilon[i,2]/g1W.total)+
                              cbind((epsilon[i,3] * ddH0W+ epsilon[i,4] * ddH1W),
                                    (epsilon[i,3]/(g0W.total^2))*ddg0,
                                    (epsilon[i,4]/(g1W.total^2))*ddg1)
                  }

                  Q <- qlogis(bound(plogis(Q), Qbounds))
                  j <- j+1
                  if(length(bestlambdas)>=j){ nextbestlambdas<-bestlambdas[j] }else{ nextbestlambdas<-0 }
            }
      }
      return(list(est=est, likelihood=likelihood, varDstar=varDstar, varIC=varIC, Qstar_best = Qstar_best))
}

construct_candidates_glmnet <- function(Y, X, Q, lambdas, ctmletype,
                                 family,Qbounds,training_set,
                                 gbound, like_type, verbose,
                                 stopFactor){
      g.dataset <- X
      Epsilon <- scorelambdas<-NULL
      bestlambdas<-NULL
      minScore <- Inf
      lambdas_remaining<-lambdas
      ncandidates<-length(lambdas)
      i <- 0
      DONE <- ncandidates<= 0
      while(!DONE) {
            i <- i+1
            if(verbose) {cat("\tBeginning construction of clever covariate", i,"\n")}

            nlambdas<-length(lambdas_remaining)
            epsilon <- NULL
            score <- rep(NA, nlambdas)
            if(ctmletype==1){ eps_try <- matrix(ncol=2,nrow=nlambdas) }
            if(ctmletype==2){ eps_try <- matrix(ncol=4,nrow=nlambdas) }
            if(verbose){ cat("\n\t  Selecting best bandwidth to add to model...\n") }

            for (j in 1:nlambdas){
                  if(verbose) { cat("\t\tTrying", lambdas_remaining[j],"\n") }
                  f1<-eval(paste("A ~ ", paste(paste(names(X[,-1]), sep=''), collapse=" + ")))
                  f2<-as.formula(f1)
                  lassox<- model.matrix(f2, X[training_set,])[,-1]
                  lassoy<-as.factor(X[training_set,1])
                  gbw<-glmnet(x=lassox,y=lassoy,family="binomial")
                  g<-predict(gbw,newx=lassox,s=lambdas_remaining[j],type="response")
                  g1W.total <- bound(g, c(gbound, 1-gbound))
                  g0W.total <- 1 - bound(g, c(gbound, 1-gbound))
                  H1W <- X[training_set,1]/g1W.total
                  H0W <- (1 - X[training_set,1])/g0W.total

                  if(ctmletype==1){
                        suppressWarnings(eps_try[j,]<- coef(glm(Y[training_set] ~ -1 +
                                                                      offset(Q[training_set, "QAW"])
                                                                + H0W + H1W, family = family)))
                        Qstar <- Q[training_set,] + c((eps_try[j,1] * H0W +
                                                             eps_try[j,2] * H1W),
                                                      eps_try[j,1]/g0W.total,
                                                      eps_try[j,2]/g1W.total)
                  }
                  if(ctmletype==2){
                        if(j==1){
                              g1<-predict(gbw,newx=lassox,s=(lambdas_remaining[j]+0.005),type="response")
                        }else{
                              g1<-predict(gbw,newx=lassox,s=(lambdas_remaining[j-1]),type="response")
                        }
                        g1W.total1 <- bound(g1, c(gbound, 1-gbound))
                        g0W.total1 <- 1 - bound(g1, c(gbound, 1-gbound))
                        ddg1<--(g1W.total-g1W.total1)
                        ddg1[which(ddg1==0)]<-1e-10
                        ddg0<--(g0W.total-g0W.total1)
                        ddg0[which(ddg0==0)]<-1e-10
                        ddH1W <- (X[training_set,1]/(g1W.total^2))*ddg1
                        ddH0W <- ((1 - X[training_set,1])/(g0W.total^2))*ddg0
                        suppressWarnings(eps_try[j,]<- coef(glm(Y[training_set] ~ -1 + offset(Q[training_set, "QAW"])+ H0W+H1W+ ddH0W + ddH1W, family = family)))
                        Qstar <- Q[training_set,] + cbind((eps_try[j,1] * H0W+ eps_try[j,2] * H1W),
                                                          eps_try[j,1]/g0W.total, eps_try[j,2]/g1W.total)+ cbind((eps_try[j,3] * ddH0W+ eps_try[j,4] * ddH1W),
                                                                                                                 (eps_try[j,3]/(g0W.total^2))*ddg0, (eps_try[j,4]/(g1W.total^2))*ddg1)
                  }

                  if(family=="binomial"){ Qstar <- plogis(Qstar)}

                  varIC <- calc_varIC(Y[training_set], Qstar, A=X[training_set,1],g1W=g1W.total)[1] * sum(training_set)/length(training_set)

                  if(like_type == "RSS"){
                        score[j] <-  sum((Y[training_set] - Qstar[,"QAW"])^2) + varIC[1]
                  }else {
                        score[j] <- -sum(Y[training_set]*log(Qstar[,"QAW"]) + (1-Y[training_set])*log(1-Qstar[,"QAW"])) + varIC[1]
                  }
            }
            # browser()
            if(verbose) { cat("\t\t",paste("penalized", like_type,":"), round(score,5), "\n") }

            score[is.nan(score)] <- Inf
            best <- which.min(abs(score))

            if(verbose) { cat("\t\tbest choice:", lambdas_remaining[best], "\n") }

            bestScore <- score[best]
            epsilon   <- rbind(epsilon, eps_try[best,])

            if(verbose) { cat("\t\tbest score:", bestScore, "\n") }

            earlyStop <- minScore*stopFactor < bestScore
            if(earlyStop & verbose){
                  cat("Stopping early because loss function of current candidate >", stopFactor, "times the best candidate seen so far\n")
                  cat("The ratio of best to current (",i,"), is ", round(bestScore/minScore, 2), "\n")
            }
            minScore <- ifelse(minScore<bestScore, minScore, bestScore)
            if(bestScore==minScore){
                  Epsilon  <- rbind(Epsilon, epsilon)
                  bestlambdas  <- c(bestlambdas, lambdas_remaining[best])

                  f1<-eval(paste("A ~ ", paste(paste(names(X[,-1]), sep=''), collapse=" + ")))
                  f2<-as.formula(f1)
                  lassox<- model.matrix(f2, X)[,-1]
                  lassoy<-as.factor(X[,1])
                  gbw<-glmnet(x=lassox,y=lassoy,family="binomial")
                  g<-predict(gbw,newx=lassox,s=lambdas_remaining[best],type="response")
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

                  if(ctmletype==2){

                        if(best==1){
                              g1<-predict(gbw,newx=lassox,s=(lambdas_remaining[best]+0.005),
                                          type="response")
                        }else{
                              g1<-predict(gbw,newx=lassox,s=(lambdas_remaining[best-1]),
                                          type="response")
                        }


                        g1W.total1 <- bound(g1, c(gbound, 1-gbound))
                        g0W.total1 <- 1 - bound(g1, c(gbound, 1-gbound))
                        ddg1<--(g1W.total-g1W.total1)
                        ddg1[which(ddg1==0)]<-1e-10
                        ddg0<--(g0W.total-g0W.total1)
                        ddg0[which(ddg0==0)]<-1e-10
                        ddH1W <- (X[,1]/(g1W.total^2))*ddg1
                        ddH0W <- ((1 - X[,1])/(g0W.total^2))*ddg0
                        Q[training_set,]  <- Q[training_set,] +
                              cbind(( Epsilon[dim(Epsilon)[1],1] * H0W[training_set] +
                                            Epsilon[dim(Epsilon)[1],2] * H1W[training_set]),
                                    Epsilon[dim(Epsilon)[1],1]/g0W.total[training_set],  Epsilon[dim(Epsilon)[1],2]/g1W.total[training_set])+
                              cbind(( Epsilon[dim(Epsilon)[1],3] * ddH0W[training_set] +
                                            Epsilon[dim(Epsilon)[1],4] * ddH1W[training_set]),
                                    ((Epsilon[dim(Epsilon)[1],3]/(g0W.total^2))*ddg0)[training_set],  ((Epsilon[dim(Epsilon)[1],4]/(g1W.total^2))*ddg1)[training_set])

                        Q[training_set,]<-qlogis(bound(plogis(Q[training_set,]), Qbounds))
                  }

                  lambdas_remaining <- lambdas_remaining[which(lambdas_remaining<lambdas_remaining[best])]
                  ncandidates<-length(lambdas_remaining)
            }else{
                  ncandidates<-0
            }
            DONE  <-  ncandidates <= 0 | earlyStop
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

      return(list(lambdas=bestlambdas, epsilon=Epsilon, earlyStop=earlyStop))
}


stage2_glmnet <- function(Y, X, Q, lambdas, ctmletype, family, Qbounds,
                   ab,training_set=rep(T,length(Y)),
                   like_type,  gbound,verbose,
                   stopFactor=10^6, best_k = NULL) {

      if(verbose) { cat ("\n\t\t-----Stage 2: Constructing candidate TMLE estimators-----\n\n") }

      candidates  <- construct_candidates_glmnet(Y, X, Q,lambdas,ctmletype,
                                          family=family[1] ,
                                          Qbounds,training_set=training_set,
                                          gbound=gbound, like_type=like_type,
                                          verbose=verbose,stopFactor=stopFactor)

      results.all <- evaluate_candidates_glmnet(Y, X, Q,lambdas,ctmletype,
                                         family[1],Qbounds, g.dataset=X,
                                         candidates, ab=ab,
                                         like_type=like_type, gbound=gbound,
                                         training_set=training_set, best_k = best_k)

      return(list(candidates=candidates, results.all=results.all))
}

cv_glmnet <- function(Y,X, est.all, Q, lambdas,ctmletype, family, Qbounds, ab, like_type, gbound,
               verbose=FALSE, PEN, V=5, folds = NULL) {

      n <- length(Y)
      nconstructed <- length(est.all)

      # Adding user-specified folds optin
      if(is.null(folds)){
            folds <- by(sample(1:n,n), rep(1:V, length=n), list)
      }else{
            if(length(folds) != V){
                  stop("The number of user-specified folds information does not match V")
            }else if(mean(sort(unlist(folds)) - 1:n)!= 0){
                  stop("Error in the indices of the user-specified folds")
            }
      }

      likelihood <- varDstar <- bias <- c(rep(0,nconstructed))
      est <- matrix(data=NA, nrow=nconstructed, ncol=V)

      if(verbose) {cat("\tfold: ")}
      for (v in 1:V){
            if(verbose) {cat(v," ")}
            test_set <- folds[[v]]
            training_set <- rep(TRUE, n)
            training_set[folds[[v]]] <- FALSE

            candidates <- stage2_glmnet(Y,X, Q,lambdas,ctmletype, family=family,
                                 Qbounds = Qbounds, ab=ab,
                                 training_set=training_set, like_type=like_type, gbound=gbound,
                                 verbose=FALSE)
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
      if(PEN){
            score <- likelihood*(diff(ab)^2) + pen
      } else {
            score <- likelihood*(diff(ab)^2)
      }
      score[is.infinite(score)|is.nan(score)] <- Inf

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
