#' Calculate CT Centrality measures with CIs
#'
#' @param drift input \eqn{p \times p} drift matrix
#' @param dt the time-interval \eqn{\Delta t}
#' @param listout if TRUE output is a list, else a vector
#' @param posterior optional. Default is NULL. Users can supply an array of dimension
#' \eqn{s \times p \times p} of samples from the posterior of the drift matrix, 
#' as output from `ctsem::ctExtract()$pop_DRIFT`. If this is supplied
#' then the function returns quantiles of centrality metrics from the posterior, as defined by `probs`
#' @param probs vector of probabilities to be passed to `quantile()`. By default set to
#'  return mean (50th percentile) and 95 percent CIs (2.5 and 97.5 percentiles)
#' 
#' @return a list or vector containing the centrality measures for each variable. If posterior is supplied,
#' then returns an array of \eqn{1 \times 2p \times 3}, with slices point estimate, lower and upper bound
#' @seealso \code{\link{DE}}, \code{\link{TE}}, \code{\link{IE}}, \code{\link{getCIs}}
#' @export
#' @examples
#' drift <- matrix(c(-.357, 0, 0,
#'                  .771, -.511, 0,
#'                  -.450, .729, -.693), 3, 3, byrow = TRUE)
#' dt = 2
#' ctCentrality(drift = drift, dt = dt)

ctCentrality <- function(drift = NULL, dt, listout = FALSE, posterior = NULL, probs = c(.025,.5,.975)){
  # Get centrality using internal function call
  if(isTRUE(listout) & !is.null(posterior)){
    listout <- FALSE
    warning("option for list output is not available when posterior supplied")
  }
  if(is.null(posterior)){
    point <- ctCentrality_int(drift = drift,
                              dt = dt,
                              listout = listout)
    return(point)
  } else {
    return(getCIs(posterior = posterior, simplify = TRUE, probs = probs, FUN = ctCentrality_int, dt = dt, listout= FALSE))
  }
 
}


