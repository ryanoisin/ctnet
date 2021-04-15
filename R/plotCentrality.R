#' Function to plot centrality metrics across a range of time-intervals
#'
#' @param drift a \eqn{p \times p} drift matrix. Can be a matrix of point estimates. Optional input, not necessary if 
#' `posterior` or `CI_obj` supplied. Defaults to NULL
#' @param posterior array of dimension \eqn{s \times p \times p} giving s posterior samples
#' of a \eqn{p \times p} drift matrix. Must be supplied if `CI_obj` absent. Default is NULL
#' @param CI_obj a \eqn{3 \times 8 \times t} array containing lower bound, upper bound, and point estimates of each centrality metric
#' across a range of \eqn{t} time-intervals, as output by the `getCIs()` function or an earlier function call. See `plot`and details.
#' @param dts a vector of \eqn{t} time-intervals. If `CI_obj` is supplied, must have some dimensions as used in generation
#' @param plot logical return a plot or array object of type `CI_obj`. See details
#' @param type character string. "TEC" plots the total effect centrality, "IEC" the indirect effect centrality and "both" both.
#' @param x scalar or vector giving the index of the variable whose centrality you wish to plot. Defaults to 1
#'
#' @return  Returns a plot by default. If  `plot` is set to FALSE and `posterior` is supplied, returns an object of type `CI_obj`. This
#' can be used to do all necessary (slow) operations first, then used as input for a future plotting call calculations
#' @seealso \code{\link{ctCentrality}}, \code{\link{plotPhi}}
#' @import graphics
#' @export

plotCentrality <- function(drift = NULL, posterior = NULL, 
                           CI_obj = NULL, dts = seq(0,1,.5), plot = TRUE, 
                           type = "TEC", x = 1){
  
  if(is.null(CI_obj)){
    if(is.null(drift) & is.null(posterior)) stop("you must supply either CI_obj or a drift matrix and/or posterior object")
    CI_obj <- sapply(dts,function(dt){
      ctCentrality(drift,dt,listout=F, posterior = posterior)
    }, simplify = "array")
  } 
  

  if(plot == FALSE){ return(CI_obj) } else {
    p <- dim(posterior)[2]
    if(type == "TEC") inds <- x
    if(type == "IEC") inds <- x + p
    if(type == "both") inds <- c(x,x + p)
    
    
    for(i in 1:length(inds)){
      ind <- inds[i]
      ylab <- dimnames(CI_obj)[[2]][ind]
      

    if(!is.null(posterior)){
      plot(dts,CI_obj[ind,], col = "red", ylab = ylab, type = "l", xlab = "Time-Interval")
    } else { 
    
    ylim = c(min(CI_obj[,ind,]), max(CI_obj[,ind,]))
    
    plot(dts,CI_obj[2,ind,], col = "red", ylab = ylab, type = "l",
         ylim = ylim)
    lines(dts, CI_obj[1,ind,], col = "red", lty = 2)
    lines(dts, CI_obj[3,ind,], col = "red", lty = 2)
    abline(h = 0)
    
    } 
      }
  
}} # eoF