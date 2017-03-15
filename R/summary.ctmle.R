#' Summarise a ctmle object
#' @param object a ctmle object
#' @param ... other parameter
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
