#' set outliers to min/max allowable values. It assumes x contains only numerical data
#' @param x input data
#' @param bounds a vector with length 2, contains the min and max of the bound
#' @return x truncated input x by min/max in bounds
#' @examples
#' x <- rnorm(1000)
#' x <- bound(x, c(-1, 1))
#' @export
#-----------------------------------------
bound <- function(x, bounds){
      x[x<min(bounds)] <- min(bounds)
      x[x>max(bounds)] <- max(bounds)
      return(x)
}
