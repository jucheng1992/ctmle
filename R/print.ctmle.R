#' print a ctmle object
#' @param x a ctmle object
#' @param ... other parameter
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
