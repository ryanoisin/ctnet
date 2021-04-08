#' Calculate CT Centrality measures
#'
#' @param drift input \eqn{p \times p} drift matrix
#' @param dt the time-interval \eqn{\Delta t}
#'
#' @return a list containing the centrality measures for each variable
#' @seealso \code{\link{DE}}, \code{\link{TE}}, \code{\link{IE}}
#' @export
#' @examples
#' drift <- matrix(c(-.357, 0, 0,
#'                  .771, -.511, 0,
#'                  -.450, .729, -.693), 3, 3, byrow = TRUE)
#' dt = 2
#' ct_centrality(drift=drift, dt=dt)

ctCentrality <- function(drift, dt){

  # Total expected influence
  indices <- seq(1:nrow(drift))

  TEC <- numeric(nrow(drift))
  IEC <- numeric(nrow(drift))

  for(i in 1:nrow(drift)){
  DVs <- indices[indices != i]
  TEC[i] <- sum(sapply(DVs,function(j){
    TE(drift = drift,IV = i,DV = j,dt = dt) # Calculate all total effects from i to p\i
  }))


  # Indirect/Mediated Expected Influence
  # Generate all IV-DV pairs not including M
  pairs <- combinat::permn(DVs)

  IEC[i] <- sum(sapply(pairs,function(l){
    IE(drift = drift,dt = dt,IV = l[1], DV = l[2], M = i)
  } ))

  }

  out <- list()
  out$TotalEffectCentrality <- TEC
  out$IndirectEffectCentrality <- IEC
  out
}


