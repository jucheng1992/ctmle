#' print a ctmle object
#' @param x a ctmle object
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
print.ctmle <- function(x, ...){
      if(identical(class(x), "ctmle")){
            cat("C-TMLE result:\n")
            cat("\tparameter estimate: ", round(x$est,5), "\n")
            cat("\testimated variance: ", round(x$var.psi,5), "\n")
            cat("\t           p-value: ", ifelse(x$pvalue <= 2*10^-16, "<2e-16",signif(x$pvalue,5)), "\n")
            cat("\t 95% conf interval:", paste("(", round(x$CI[1],5), ", ", round(x$CI[2],5), ")", sep=""),"\n")
      } else {
            stop("Error calling print.ctmle. 'x' needs to have class 'ctmle'")
      }
}
