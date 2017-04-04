#' Help function to build the sequence of gn candidates in ctmleGeneral
#'
#' This function helps building gn candidates for ctmleGeneral. It returns gn_candidates_cv,
#' gn_candidates, and folds, which could be directly applied to ctmleGeneral.
#'
#' @export
#' @param A binary treatment indicator, 1 - treatment, 0 - control
#' @param W vector, matrix, or dataframe containing baseline covariates for Q bar
#' @param SL.library a vector of the names of the estimators for ctmle
#' (need to be prepared in the format for SL, see more details in SuperLearner package),
#' The theory of ctmle requires the estimators are ordered by the model complexity,
#' with the last one be a consistent estimator.
#' @param folds The list of indices for the ctmle cross-validation step
#' @return gn_candidates_cv matrix or dataframe, each column stand for a estimate of
#' propensity score. Estimate in the column with larger index should have smaller empirical loss
#' @return gn_candidates matrix or dataframe, each column stand for a the cross-validated estimate.
#' For example, the (i,j)-th element is the predicted propensity score by j-th estimator,
#' for i-th observation, when it is in the validation set
#' @return folds The list of indices for the ctmle cross-validation step
#' @examples
#'N <- 1000
#'p = 100
#'V = 5
#'Wmat <- matrix(rnorm(N * p), ncol = p)
#'gcoef <- matrix(c(-1,-1,rep(-(3/((p)-2)),(p)-2)),ncol=1)
#'
#'W <- as.data.frame(Wmat)
#'g <- 1/(1+exp(Wmat%*%gcoef / 3))
#'A <- rbinom(N, 1, prob = g)
#'
#'folds <-by(sample(1:N,N), rep(1:V, length=N), list)
#'
#'lasso_fit <- cv.glmnet(x = as.matrix(W), y = A, alpha = 1, nlambda = 100, nfolds = 10)
#'lasso_lambdas <- lasso_fit$lambda[lasso_fit$lambda <= lasso_fit$lambda.min][1:5]

#'# Build template for glmnet
#'SL.glmnet_new <- function (Y, X, newX, family, obsWeights, id, alpha = 1,
#'                           nlambda = 100, lambda = 0,...)
#'{
#'     # browser()
#'     if (!is.matrix(X)) {
#'           X <- model.matrix(~-1 + ., X)
#'          newX <- model.matrix(~-1 + ., newX)
#'    }
#'    fit <- glmnet::glmnet(x = X, y = Y,
#'                          lambda = lambda,
#'                          family = family$family, alpha = alpha)
#'    pred <- predict(fit, newx = newX, type = "response")
#'      fit <- list(object = fit)
#'    class(fit) <- "SL.glmnet"
#'    out <- list(pred = pred, fit = fit)
#'    return(out)
#'}
#'
#'# Use a sequence of estimator to build gn sequence:
#'SL.cv1lasso <- function (... , alpha = 1, lambda = lasso_lambdas[1]){
#'    SL.glmnet_new(... , alpha = alpha, lambda = lambda)
#'}
#'
#'SL.cv2lasso <- function (... , alpha = 1, lambda = lasso_lambdas[2]){
#'     SL.glmnet_new(... , alpha = alpha, lambda = lambda)
#'}
#'
#'SL.cv3lasso <- function (... , alpha = 1, lambda = lasso_lambdas[3]){
#'     SL.glmnet_new(... , alpha = alpha, lambda = lambda)
#'}
#'
#'SL.cv4lasso <- function (... , alpha = 1, lambda = lasso_lambdas[4]){
#'      SL.glmnet_new(... , alpha = alpha, lambda = lambda)
#'}
#'
#'SL.library = c('SL.cv1lasso', 'SL.cv2lasso', 'SL.cv3lasso', 'SL.cv4lasso', 'SL.glm')
#'
#'gn_seq <- build_gn_seq(A = A, W = W, SL.library = SL.library, folds = folds)
#'
#'gn_seq$gn_candidates_cv
#'gn_seq$gn_candidates
build_gn_seq <- function(A, W, SL.library, folds, verbose = TRUE){

      SL_fit <- SuperLearner(Y = A, X = as.data.frame(W),
                             SL.library = SL.library, family = 'binomial', verbose = verbose,
                             cvControl = list(V = length(folds), validRows = folds))

      return(list(gn_candidates_cv = SL_fit$Z, gn_candidates = SL_fit$library.predict, folds = folds))

}
