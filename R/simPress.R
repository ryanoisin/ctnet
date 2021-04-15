#' Calculate the effect of a press intervention
#'
#' @param drift input \eqn{p \times p} drift matrix
#' @param dt the time-interval considered \eqn{\Delta t}
#' @param M index or indices corresponding to the intervened-on variables
#' @param press scalar or vector giving the value of the press intervention(s)
#' @param start (optional) full vector of starting values.
#'
#' @return a vector representing the value of the CT process after a press intervention to set M to value press
#' @seealso \code{\link{DE}}
#' @export
#' @examples
#' n <- 4
#' drift <- matrix(c(-1.2,  .25,     0,  1,
#'                  1.1, -.5,     0,   0,
#'                    0,  1.18, -1.2,   0 ,
#'                  0, -1.46,     .5, -1.2),n,n, byrow = TRUE)
#'
#' M = 2; press = 1; start = c(0,M,0,0)
#'  # Set up the time vector
#'    dts2 <- seq(0,5,.001)
#'  # Calculate effect over a range of time intervals
#'    y_ct <- matrix(0,length(dts2),n+1)
#'
#' # Get Trajectories for CT and total (DT2) effects
#'
#'  for(i in 1:length(dts2)){
#'  # Total Effect
#'  y_ct[i,1:n] <- simPress(drift,dt=dts2[i],M,press, start)
#'  y_ct[i,n+1] <- dts2[i]
#' }
#'
#' # Create plot
#' cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")
#' plot.new()
#' plot.window(xlim=c(0, 5), ylim=c(-1,1))
#' for(c in 1:4) lines(y = y_ct[,c], x = y_ct[,5], type = "l", col = cols[c], lwd = 2)
#' abline(h = 0); axis(1); axis(2); title(xlab = "Time (t)", ylab = "Variable Values")
#'
#' # To get the new equilibrium plug in a large value for Delta t
#' simPress(drift,dt=100,M,press, start)
#'

simPress <- function(drift,dt,M,press, start = NULL){
  n <- nrow(drift)

  # Prepare starting values vector

  if(is.null(start)){
    start <- rep(0,n)
  }

  start0 <- start
  start[M] <- press

  # Transformation matrix
  D <- diag(nrow(drift))
  D[M,M] <- 0

  # Under press intervention on m, m is "exogenous" to all other variables
  # We delete the m row (incoming effects) but retain the m column (outgoing effects)
  drift_tilde <- D%*%drift

  if(any(Re(eigen(drift_tilde)$values) > 0)){
    stop("Intervening on this variable makes the system unstable!")
  }

  expm::expm(drift_tilde*dt)%*%start

}
