#' @title
#' Plotting the results of Sobol SA
#'
#' @description
#' \code{plot.sobolRes_ats} plots the results of Sobol SA for objects of class 
#' \code{sobolRes_ats}.
#'
#' @details
#' 1st order and total Sobol SA indices are plotted for each input
#' parameter against time.
#'
#' @param x [\code{sobolRes_ats}]\cr
#'   resulting output of \code{\link{ODEsobol_ats}}, of class 
#'   \code{sobolRes_ats}.
#' @param type [\code{character(1)}]\cr
#'   plot type, i.e. \code{"p", "l", "b", "c", "o", "s", "h"} or \code{"n"}. 
#'   Defaults to \code{"b"}.
#' @param legendPos [\code{character(1)}]\cr
#'   legend position, default is \code{"topleft"}.
#' @param overall_main [\code{character(1)}]\cr
#'   the common title for the two graphics. Default is \code{NULL}, which means
#'   an automatic title is generated.
#' @param ... additional arguments passed to \code{\link{plot}}.
#'
#' @return TRUE (invisible; for testing purposes).
#'
#' @method plot sobolRes_ats
#'
#' @note Unfortunately, the passing of arguments (e.g. "main") does not work
#'   correctly.
#'
#' @examples
#'
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
#' FHNyini  <- c(Voltage = -1, Current = 1)
#' FHNtimes <- seq(0.1, 50, by = 5)
#'
#' FHNres_ats <- ODEsobol_ats(mod = FHNmod,
#'                            pars = c("a", "b", "s"),
#'                            yini = FHNyini,
#'                            times = FHNtimes,
#'                            y_idx = 1,            # only Voltage
#'                            seed = 2015,
#'                            n = 10,               # use n >> 10!
#'                            rfuncs = c("runif", "runif", "rnorm"),
#'                            rargs = c(rep("min = 0.18, max = 0.22", 2),
#'                                  "mean = 3, sd = 0.2 / 3"),
#'                            method = "martinez",
#'                            nboot = 0)
#'
#' # Plot:
#' plot(FHNres_ats, type = "l", legendPos = "topright")
#'
#' @author Frank Weber
#' @seealso \code{\link{ODEsobol_ats}, \link[sensitivity]{sobol}, 
#' \link[sensitivity]{soboljansen_matrix}, 
#' \link[sensitivity]{sobolmartinez_matrix}}
#'
#' @export
#' @import checkmate
#'

plot.sobolRes_ats <- function(x, type = "b", legendPos = "topleft",
                              overall_main = NULL, ...) {

  ##### Plausibilitaet #################################################
  assertClass(x, "sobolRes_ats")
  assertCharacter(type, len = 1)
  notOk <- !type %in% c("p", "l", "b", "c", "n", "o", "s", "h")
  if(notOk)
    stop(paste("type must be one of \"p\", \"l\", \"b\", \"c\", \"n\",",
               "\"o\", \"s\" or \"h\"!"))
  assertCharacter(legendPos, len = 1)
  notOk <- !any(rep(legendPos, 9) == c("bottomright", "bottom", "bottomleft", 
              "left", "topleft", "top", "topright", "right", "center"))
  if(notOk)
    stop("legendPos must be one of \"bottomright\", \"bottom\",
         \"bottomleft\", \"left\", \"topleft\", \"top\", \"topright\",
         \"right\", \"center\"!")
  stopifnot(is.character(overall_main) && length(overall_main) == 1 ||
              is.null(overall_main))

  ##### Vorbereitungen #################################################
  k <- nrow(x$S) - 1
  pars <- rownames(x$S)[- 1]
  parsCols <- rainbow(k)
  # Extrema SA Indizes:
  minMaxS <- c(0.95 * min(x$S[-1, ]), 1.05 * max(x$S[-1, ]))
  minMaxT <- c(0.95 * min(x$T[-1, ]), 1.05 * max(x$T[-1, ]))
  # Gemeinsamer Titel fuer beide Grafiken:
  if(is.null(overall_main)){
    overall_main <- paste0("Sobol sensitivity indices for y_idx = ", x$y_idx,
                           " and method = \"", x$method, "\"")
  }
  
  oldpar <- par(mfrow = c(1, 2), mar = c(4, 4, 1, 2) + 0.2,
                oma = c(0, 0, 2, 0))
  
  ##### 1st order SA indices ###########################################
  # Plot fuer ersten Parameter:
  plot(x = x$S[1, ], y = x$S[2, ],
       xlab = "Time", ylab = "1st order Sobol SA indices",
       type = type, col = parsCols[1], ylim = minMaxS, ...)
  # Plots fuer alle weiteren Parameter:
  for(i in 2:k) {
    lines(x = x$S[1, ], y = x$S[i + 1, ], type = type, col = parsCols[i], ...)
  }
  # Legende:
  if(type %in% c("b", "o")){
    legend(legendPos, legend = pars, col = parsCols, bg = "white",
           lty = 1, pch = 1)
  } else if(type %in% c("l", "c", "s", "h")){
    legend(legendPos, legend = pars, col = parsCols, bg = "white", lty = 1)
  } else if(type == "p"){
    legend(legendPos, legend = pars, col = parsCols, bg = "white", pch = 1)
  }
  
  ##### total SA indices ###############################################
  # Plot fuer ersten Parameter:
  plot(x = x$T[1, ], y = x$T[2, ],
       xlab = "Time", ylab = "Total Sobol SA indices",
       type = type, col = parsCols[1], ylim = minMaxT, ...)
  # Plots fuer alle weiteren Parameter:
  for(i in 2:k) {
    lines(x = x$T[1, ], y = x$T[i + 1, ], type = type, col = parsCols[i], ...)
  }
  # Legende:
  if(type %in% c("b", "c", "o")){
    legend(legendPos, legend = pars, col = parsCols, bg = "white",
           lty = 1, pch = 1)
  } else if(type %in% c("l", "s", "h")){
    legend(legendPos, legend = pars, col = parsCols, bg = "white", lty = 1)
  } else if(type == "p"){
    legend(legendPos, legend = pars, col = parsCols, bg = "white", pch = 1)
  }
  
  # Gemeinsamer Titel:
  mtext(overall_main, side = 3, line = 0, outer = TRUE, cex = 1.2, font = 2)
  
  par(oldpar)
  return(invisible(TRUE))
}
