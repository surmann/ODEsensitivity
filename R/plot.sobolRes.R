#' @title
#' Plotting the results of Sobol SA
#'
#' @description
#' \code{plot.sobolRes} plots the results of Sobol SA.
#'
#' @details
#' 1st order and total Sobol SA indices are plotted for each input
#' parameter against time.
#'
#' @param res [\code{sobolRes}]\cr
#'   resulting output of \code{\link{ODEsobol}}, of class \code{sobolRes}.
#' @param type [\code{character(1)}]\cr
#'   plot type, i.e. \code{"p", "l", "b", "c"} or \code{"n"}.
#' @param legendPos [\code{character(1)}]\cr
#'   legend position, default is \code{"topleft"}.
#' @param ... additional arguments.
#'
#' @return NULL
#'
#' @method plot sobolRes
#'
#' @examples
#'
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
#' FHNtimes <- seq(0.1, 20, by = 0.5)
#'
#' FHNres <- ODEsobol(mod = FHNmod,
#'                    pars = names(FHNpars),
#'                    yini = FHNyini,
#'                    times = FHNtimes,
#'                    seed = 2015,
#'                    n = 10,                        # use n >> 10!
#'                    trafo = function(Y) Y[, 1],    # voltage only
#'                    ncores = 4)
#'
#' # Plot:
#' plot(FHNres, type = "l", legendPos = "topright")
#'
#' @author Stefan Theers
#' @seealso \code{\link{ODEsobol}},
#'   \code{\link[sensitivity]{sobol}},
#'   \code{\link[sensitivity]{sobol2007}}
#'
#' @examples
#'
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
#' FHNtimes <- seq(0.2, 30, by = 0.2)
#'
#' FHNres <- ODEsobol(mod = FHNmod,
#'                    pars = names(FHNpars),
#'                    yini = FHNyini,
#'                    times = FHNtimes,
#'                    seed = 2015,
#'                    n = 10,                        # use n >> 10!
#'                    trafo = function(Y) Y[, 1],    # voltage only
#'                    ncores = 4)
#'
#' plot(FHNres, type = "l", legendPos = "topright")
#'
#' @export
#'
#' @import checkmate
#'
plot.sobolRes <- function(res, type = "b",
                          legendPos = "topleft", ...) {

  ##### Plausibilitaet #################################################
  assertClass(res, "sobolRes")
  assertCharacter(type, len = 1)
  notOk <- !any(rep(type, 5) == c("p", "l", "b", "c", "n"))
  if(notOk)
    stop("type must be one of \"p\", \"l\", \"b\", \"c\" or \"n\"!")
  assertCharacter(legendPos, len = 1)
  notOk <- !any(rep(legendPos, 9) == c("bottomright", "bottom",
    "bottomleft", "left", "topleft", "top", "topright", "right",
    "center"))
  if(notOk)
    stop("legendPos must be one of \"bottomright\", \"bottom\",
      \"bottomleft\", \"left\", \"topleft\", \"top\", \"topright\",
      \"right\", \"center\"!")

  ##### Vorbereitungen #################################################
  k <- nrow(res$S) - 1
  pars <- rownames(res$S)[- 1]
  parsCols <- rainbow(k)
  # Extrema SA Indizes:
  minMaxS <- c(0.95 * min(res$S[-1, ]), 1.05 * max(res$S[-1, ]))
  minMaxT <- c(0.95 * min(res$T[-1, ]), 1.05 * max(res$T[-1, ]))

  par(mfrow = c(1, 2))

  ##### 1st order SA indices ###########################################
  # Plot fuer ersten Parameter:
  plot(x = res$S[1, ], y = res$S[2, ],
       main = "1st order Sobol SA indices", xlab = "time",
       ylab = "1st order Sobol SA indices",
       type = type, col = parsCols[1], ylim = minMaxS, ...)
  # Plots fuer alle weiteren Parameter:
  for(i in 2:k) {
    lines(x = res$S[1, ], y = res$S[i + 1, ],
          type = type, col = parsCols[i], ...)
  }
  # Legende:
  legend(legendPos, legend = pars, col = parsCols, bg = "white",
         lty = 1, pch = 1)

  ##### total SA indices ###############################################
  # Plot fuer ersten Parameter:
  plot(x = res$T[1, ], y = res$T[2, ],
       main = "Total Sobol SA indices", xlab = "time",
       ylab = "Total Sobol SA indices",
       type = type, col = parsCols[1], ylim = minMaxT, ...)
  # Plots fuer alle weiteren Parameter:
  for(i in 2:k) {
    lines(x = res$T[1, ], y = res$T[i + 1, ],
          type = type, col = parsCols[i], ...)
  }
  # Legende:
  legend(legendPos, legend = pars, col = parsCols, bg = "white",
         lty = 1, pch = 1)

  par(mfrow = c(1, 1))
  return(invisible(TRUE))
}
