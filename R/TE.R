#' CT Total Effect
#'
#' @param drift input \eqn{p \times p} drift matrix
#' @param dt the time-interval \eqn{\Delta t}
#' @param IV index corresponding to the independent-variable
#' @param DV index corresponding to the dependent-variable
#'
#' @return a scalar total effect
#' @export
#' @examples
#' drift<- matrix(c(-.357, 0, 0,
#'                  .771, -.511, 0,
#'                  -.450, .729, -.693), 3, 3, byrow = TRUE)
#' dt = 2
#' IV = 1
#' DV = 3
#' TE(drift=drift,dt=dt,IV=IV,DV=DV)

TE <- function(drift,dt,IV,DV){
      expm::expm(drift*dt)[DV,IV]
  }
