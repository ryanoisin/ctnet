#' Calculate CT Centrality measures (internal function)
#'
#' @param drift input \eqn{p \times p} drift matrix
#' @param dt the time-interval \eqn{\Delta t}
#' @param listout if TRUE output is a list, else a vector
#' 
#' @return a list or vector containing the centrality measures for each variable.
#' @seealso \code{\link{DE}}, \code{\link{TE}}, \code{\link{IE}}
#' @keywords Internal

ctCentrality_int <- function(drift, dt, listout = TRUE){
  
  if(is.null(dimnames(drift))){
    labs <- paste0("Y",1:nrow(drift))
  } else{ labs <- colnames(drift)}
  
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
    
    names(TEC) <- names(IEC) <- labs
    
  }
  if(isTRUE(listout)){
    out <- list()
    out$TotalEffectCentrality <- TEC
    out$IndirectEffectCentrality <- IEC
  } else{
    out <- c(TEC,IEC)
    names(out) <- c(paste0("TEC_",labs),paste0("IEC_",labs))
  }
  out
}


