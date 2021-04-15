#' Visualize the model-implied lagged effects across a range of time-intervals
#'
#' @param drift a \eqn{p \times p} drift matrix. Can be a matrix of point estimates. Optional input, not necessary if 
#' `posterior` or `CI_obj` supplied. Defaults to NULL
#' @param posterior array of dimension \eqn{s \times p \times p} giving s posterior samples
#' of a \eqn{p \times p} drift matrix. Must be supplied if `CI_obj` absent. Default is NULL
#' @param CI_obj a \eqn{3 \times p^2 \times t} array containing lower bound, upper bound, and point estimates of each lagged parameter
#' across a range of \eqn{t} time-intervals, as output by the `getCIs()` function or an earlier function call. See `plot`and details.
#' @param dts a vector of \eqn{t} time-intervals. If `CI_obj` is supplied, must have some dimensions as used in generation
#' @param plot logical return a plot or array object of type `CI_obj`. See details
#' @param index matrix or character vector describing which parameters to plot. The matrix must be of dimension \eqn{q \times 2}
#' giving the row and column of the `q` parameters to be plotted. Character strings accepted are: \cr
#' "all": plot all lagged parameters \cr
#' "AR": plot only auto-regressive (diagonal) parameters \cr
#' "CL": plot only cross-lagged (off-diagonal) parameters
#' @param colvec vector of length `q` of colors to be used for plotting. Defaults to NULL. If not supplied, function will plot
#' a spectrum of `q` colors from blue to red using `colorRampPallete()`
#' @param leg logical. Include a legend or not. Defaults to TRUE
#'
#' @return  Returns a plot by default. If  `plot` is set to FALSE and `posterior` is supplied, returns an object of type `CI_obj`. This
#' can be used to do all necessary (slow) operations first, then used as input for a future plotting call calculations
#' @seealso \code{\link{plotCentrality}}
#' @import graphics grDevices
#' @export

plotPhi <-
  function(drift = NULL,
           posterior = NULL,
           CI_obj = NULL,
           dts = seq(0, 1, .5),
           plot = TRUE,
           index = "all",
           colvec = NULL,
           leg = TRUE) {
    if (is.null(CI_obj)) {
      if (is.null(drift) &
          is.null(posterior))
        stop("you must supply either CI_obj or a drift matrix and/or posterior object")
      CI_obj <- sapply(dts, function(dt) {
        getCIs(posterior,
               simplify = TRUE,
               FUN = expm::expm,
               const = dt)
      }, simplify = "array")
      
    }
    
    if (plot == FALSE) {
      return(CI_obj)
    } else {
      # extract number of dimensions
      p <- sqrt(dim(CI_obj)[2])
      
      # get matrix of indices
      if (is.matrix(index) ||
          is.data.frame(index)) {
        m <- index
      } else {
        if (index == "all")
          m <- expand.grid(1:p, 1:p)
        if (index == "AR")
          m <- cbind(1:p, 1:p)
        if (index == "CL")
          m <- expand.grid(1:p, 1:p)[-seq(1, p ^ 2, p + 1), ]
      }
      # So that we don't have to transpose CI_obj, need to map matrix indices to column numbers
      mcomp <- cbind(expand.grid(1:p, 1:p)[, c(1, 2)], 1:p ^ 2)
      
      ind <- rep(NA, nrow(m))
      for (i in 1:nrow(m)) {
        for (j in 1:nrow(mcomp)) {
          if (all(m[i, c(1, 2)] == mcomp[j, c(1, 2)])) {
            ind[i] <- mcomp[j, 3]
          } else
            next
        }
      }
      
      # create colors automatically if not supplied
      if (is.null(colvec))
        colvec <- colorRampPalette(c("blue", "red"))(nrow(m))
      
      # Now start plotting
      ylim = c(min(CI_obj[, ind, ]), max(CI_obj[, ind, ]))
      
      plot.new()
      plot.window(xlim = c(dts[1], max(dts)), ylim = ylim)
      axis(1)
      axis(2)
      title(main = "Lagged effects as a function of Time-Interval",
            xlab = "Time-Interval",
            ylab = "Effect size")
      abline(h = 0, col = "grey")
      
      for (i in 1:length(ind)) {
        lines(dts, CI_obj[2, ind[i], ], col = colvec[i])
        lines(dts, CI_obj[1, ind[i], ], col = colvec[i], lty = 2)
        lines(dts, CI_obj[3, ind[i], ], col = colvec[i], lty = 2)
      }
      
      if (leg == TRUE) {
        legtext <- as.vector(unlist(apply(m, 1, function(row) {
          bquote(Phi[.(paste0(row[1], row[2]))])
        })))
        
        
        legend(
          "topright",
          as.expression(legtext),
          lty = 1,
          col = colvec,
          cex = 1,
          xpd = TRUE,
          bty = "n"
        )
      }
    }
    
  }