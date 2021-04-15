#' Calculate CT Direct Effect \eqn{DE(\Delta t)}
#'
#' @param drift input \eqn{p \times p} drift matrix
#' @param dt the time-interval \eqn{\Delta t}
#' @param IV index corresponding to the independent-variable
#' @param DV index corresponding to the dependent-variable
#' @param M index or indices corresponding to the mediators to be intervened on
#'
#' @return a scalar direct effect corresponding to a continuous intervention on M
#' @export
#' @examples
#' drift<- matrix(c(-.357, 0, 0,
#'                  .771, -.511, 0,
#'                  -.450, .729, -.693), 3, 3, byrow = TRUE)
#' dt = 2
#' IV = 1
#' DV = 3
#' M = 2
#' DE(drift = drift, dt = dt, IV = IV, DV = DV, M = M)

DE <- function(drift,dt,IV,DV,M){
  S <- diag(nrow(drift))
  S[M,M] <- 0
  drift_tilde <- S%*%drift%*%S # Remove drift parameters relating to M

  outmat <- expm::expm(drift_tilde * dt)
  outmat[DV,IV]
}
