#' Calculate CT Indirect Effect \eqn{IE(\Delta t)}
#'
#' @param drift input \eqn{p \times p} drift matrix
#' @param dt the time-interval \eqn{\Delta t}
#' @param IV index corresponding to the independent-variable
#' @param DV index corresponding to the dependent-variable
#' @param M index or indices corresponding to the mediators
#'
#' @return a scalar direct effect corresponding to a continuous intervention on M
#' @seealso \code{\link{DE}}, \code{\link{TE}}
#' @export
#' @examples
#' drift<- matrix(c(-.357, 0, 0,
#'                  .771, -.511, 0,
#'                  -.450, .729, -.693), 3, 3, byrow = TRUE)
#' dt = 2
#' IV = 1
#' DV = 3
#' M = 2
#' occ <- 1
#' IE(drift=drift,dt=dt,IV=IV,DV=DV,M=M)

IE <- function(drift,dt,IV,DV,M){
  TE(drift=drift,dt=dt,IV=IV,DV=DV) - DE(drift=drift,dt=dt,IV=IV,DV=DV,M=M)
}

