#' @title
#' Plotting the results of Morris SA
#'
#' @description
#' \code{plot} plots the results of Morris SA.
#'
#' @details
#' \code{plot} with \code{type = "sep"} plots mu.star and
#'   sigma separately versus time.
#'
#' \code{plot} with \code{type = "trajec"} plots mu.star versus
#'   sigma for every point of time.
#'
#' @param x [\code{morrisRes}]\cr
#'   resulting output of \code{\link{ODEmorris}}, of class \code{morrisRes}.
#' @param type [\code{character(1)}]\cr
#'   plot type, choose between \code{"sep"} and \code{"trajec"}.
#' @param legendPos [\code{character(1)}]\cr
#'   legend position, default is \code{"topleft"}.
#' @param ... additional arguments passed to \code{\link{plot}}.
#'
#' @return TRUE (invisible; for testing purposes).
#'
#' @method plot morrisRes
#'
#' @examples
#' ##### FitzHugh-Nagumo equations (Ramsay et al., 2007)
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
#' FHNpars  <- c(a = 0.2,     # parameter a
#'               b = 0.3,     # parameter b
#'               s = 3)       # parameter s (= c in the original notation)
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
#' @note Unfortunately, the passing of arguments (e.g. "main") does not work
#'   correctly.
#'
#' @author Stefan Theers
#' @seealso \code{\link{ODEmorris}},
#'   \code{\link[sensitivity]{morris}}
#'
#' @export
#' @import
#'   checkmate
#'

plot.morrisRes <- function(x, type = "sep", legendPos = "topleft", ...) {

  ##### Check input #################################################
  assertClass(x, "morrisRes")
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
  if(type == "sep")    plotSep(x$res, x$pars, legendPos, ...)
  if(type == "trajec") plotTrajectories(x$res, x$pars, legendPos, ...)
  
  # For testing purposes:
  return(invisible(TRUE))
}
