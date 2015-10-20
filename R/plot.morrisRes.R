##### auxiliary function: plotting mu and sigma seperately ###########
plotSep <- function(res, pars, legendPos, ...) {
  t.vec <- res[1, ]
  k     <- (nrow(res) - 1) / 3
  my.cols <- rainbow(k)
  # Teile den Plot in mu.star und sigma:
  par(mfrow = c(1, 2))
  # mu.star:
  plot(t.vec, y = res[k + 2, ], type = "l", col = my.cols[1], lwd = 1,
       main = "mu.star(t)",
       ylim = c(min(res[(k+2):(2*k+1), ], na.rm = TRUE),
                max(res[(k+2):(2*k+1), ], na.rm = TRUE)),
       xlab = "time t", ylab = "mu.star")
  if(k > 1) {
    j <- 1
    for(i in (k+2+1):(2*k+1)) {
      lines(t.vec, y = res[i, ], col = my.cols[j], lwd = 1, type = "l")
      j <- j + 1
    }
  }
  legend("topleft", legend = pars, lty = 1,
         col = my.cols)
  # sigma:
  plot(t.vec, y = res[2+2*k, ], type = "l", col = my.cols[1], lwd = 1,
       main = "sigma(t)",
       ylim = c(min(res[(2*k+2):(3*k+1), ], na.rm = TRUE),
                max(res[(2*k+2):(3*k+1), ], na.rm = TRUE)),
       xlab = "time t", ylab = "sigma")
  if(k > 1) {
    j <- 1
    for(i in (2*k+2):(3*k+1)) {
      lines(t.vec, y = res[i, ], col = my.cols[j], lwd = 1, type = "l")
      j <- j + 1
    }
  }
  legend(legendPos, legend = pars, lty = 1, col = my.cols)
  par(mfrow = c(1, 1))
}

##### auxiliary function: plotting trajectories ######################
plotTrajectories <- function(res, pars, legendPos, ...) {
  t.vec <- res[1, ]
  k     <- (nrow(res) - 1) / 3
  my.cols <- rainbow(k)
  # Zeichne Trajektoren:
  plot(x = res[k+2, ], y = res[2+2*k, ], type = "b", col = my.cols[1], lwd = 1,
       main = "Trajectories: mu.star vs. sigma",
       xlab = "mu.star", ylab = "sigma")
  if(k > 1) {
    j <- 1
    for(i in (k+2):(2*k+1)) {
      lines(x = res[i, ], y = res[i+k, ], col = my.cols[j], lwd = 1, type = "b")
      j <- j + 1
    }
  }
  legend(legendPos, legend = pars, lty = 1, col = my.cols)
}

#' @title
#' Plotting the Results of Morris SA
#'
#' @description
#' \code{plot} plots the results of Morris SA.
#'
#' @details
#' \code{plot} with \code{type = "sep"} plots mu.star and
#'   sigma seperately versus time.
#'
#' \code{plot} with \code{type = "trajec"} plots mu.star versus
#'   sigma for every point of time.
#'
#' @param res [\code{morrisRes}]\cr
#'   resulting output of \code{\link{ODEmorris}}, of class \code{morrisRes}.
#' @param type [\code{character(1)}]\cr
#'   plot type, choose between \code{"sep"} and \code{"trajec"}.
#' @param legendPos [\code{character(1)}]\cr
#'   legend position, default is \code{"topleft"}.
#' @param ... additional arguments.
#'
#' @return NULL
#'
#' @method plot morrisRes
#'
#' @examples
#' ##### FitzHugh-Nagumo equations (Ramsay et al, 2007)
#' # definition of the model itself, parameters, initial values
#' # and the times vector:
#' FHNmod <- function(Time, State, Pars) {
#'   with(as.list(c(State, Pars)), {
#'
#'     dVoltage <- s * (Voltage - Voltage^3 / 3 + Current)
#'     dCurrent <- - 1 / s *(Voltage - a + b * Current)
#'
#'     return(list(c(dVoltage, dCurrent)))
#'   })
#' }
#'
#' FHNpars  <- c(a = 0.2,     # paramter a
#'               b = 0.3,     # paramter b
#'               s = 3)       # paramter s (= c in the original notation)
#'
#' FHNyini  <- c(Voltage = -1, Current = 1)
#' FHNtimes <- seq(0.1, 100, by = 10)
#'
#' FHNres <- ODEmorris(mod = FHNmod,
#'                     pars = names(FHNpars),
#'                     yini = FHNyini,
#'                     times = FHNtimes,
#'                     seed = 2015,
#'                     binf = c(0.18, 0.18, 2.8),
#'                     bsup = c(0.22, 0.22, 3.2),
#'                     r = 25,
#'                     design =
#'                         list(type = "oat", levels = 100, grid.jump = 1),
#'                     trafo = function(Y)  Y[, 1])    # voltage only
#'
#' # Plots:
#' plot(FHNres, type = "sep")
#' plot(FHNres, type = "trajec")
#'
#' @author Stefan Theers
#' @seealso \code{\link{ODEmorris}},
#'   \code{\link[sensitivity]{morris}}
#'
#' @export
#' @import
#'   checkmate
#'
plot.morrisRes <- function(res, type = "sep", legendPos = "topleft",
                           ...) {

  ##### Plausibilitaet #################################################
  assertClass(res, "morrisRes")
  assertCharacter(type, len = 1)
  notOk <- !any(rep(type, 2) == c("sep", "trajec"))
  if(notOk)
    stop("type must be one of \"sep\" or \"trajec\"!")
  assertCharacter(legendPos, len = 1)
  notOk <- !any(rep(legendPos, 9) == c("bottomright", "bottom",
    "bottomleft", "left", "topleft", "top", "topright", "right",
    "center"))
  if(notOk)
    stop("legendPos must be one of \"bottomright\", \"bottom\",
      \"bottomleft\", \"left\", \"topleft\", \"top\", \"topright\",
      \"right\", \"center\"!")

  ##### Plot ###########################################################
  if(type == "sep")    plotSep(res[[1]], res[[2]], legendPos, ...)
  if(type == "trajec") plotTrajectories(res[[1]], res[[2]], legendPos,
                                        ...)
  # for testing purposes:
  return(invisible(TRUE))
}
