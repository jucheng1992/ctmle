#' Summarise a ctmle object
#' @param object a ctmle object
#' @param ... other parameter
#' @examples
#'\dontrun{
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
#'ctmle_summary = summary(ctmle_discrete_fit1)
#'ctmle_summary
#'ctmle_discrete_fit1
#'}
#' @export
summary.ctmle <- function(object, ...){
      if(identical(class(object), "ctmle")){
            d <- ncand <- terms <- covar <- NULL
            if(!is.null(object$candidates)){
                  ncand <- sum(!is.na(object$candidates$candidates$terms))
                  terms <- c(object$candidates$candidates$terms[-1])
                  covar <- object$candidates$candidates$covar[1:ncand]
                  d <- data.frame(c("(intercept)",object$candidates$candidates$terms[-1])[1:ncand],
                                  object$candidates$candidates$covar[1:ncand],
                                  object$candidates$results.all$est[1:ncand],
                                  object$cv.res$likelihood[1:ncand],    # this is how the thing was selected.
                                  object$cv.res$varIC[1:ncand],
                                  object$cv.res$penlikelihood[1:ncand])
                  dimnames(d) <- list(paste("cand", 1:ncand),
                                      c("term added", "cleverCovar","estimate", "cv-RSS", "cv-varIC", "cv-penRSS"))
            }
            names(object$best_k) <- names(object$est) <- NULL
            summary.ctmle <- list(est=object$est, var=object$var.psi,
                                  CI=object$CI,pvalue=object$pvalue, d=d,
                                  selected=object$best_k, ncand=ncand,
                                  terms=terms,
                                  covar=covar)

      } else {
            summary.ctmle <- NULL
      }
      class(summary.ctmle) <- "summary.ctmle"
      return(summary.ctmle)
}
