#' Calculate confidence/credible intervals of arbitrary functions of the drift matrix
#'
#' @param posterior array of dimension \eqn{s \times p \times p} giving s posterior samples
#' of a \eqn{p \times p} drift matrix
#' @param simplify argument to be passed to sapply. Defaults to TRUE
#' @param probs vector of quantiles to be computed. By default returns the 2.5, 50 and 97.5 percentile
#' @param const any number with which to multiply the drift matrix. Default to 1. Should only be used
#' in combination with the expm function, to allow for computation of \eqn{\Phi(\Delta t = const)}
#' @param FUN any function that takes a drift matrix as input and returns a vector
#' @param ... arguments to be passed to FUN
#' 
#' @return a matrix or vector containing lower bound, point estimate and upper bound quantiles
#' @seealso \code{\link{DE}}, \code{\link{TE}}, \code{\link{IE}}
#' @export

getCIs <- function(posterior, simplify = TRUE, probs = c(.025,.5,.975), const = 1, FUN,  ...){
  samps <- dim(posterior)[1]
  po <- sapply(1:samps,function(s){
    FUN(posterior[s,,]*const, ...)
  }, simplify = simplify)
  if(!is.null(nrow(po))){
    out <- apply(po,1,stats::quantile,probs = probs)
  } else {
    out <- stats::quantile(po,probs = probs)
  }
  out
}
