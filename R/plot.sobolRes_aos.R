#' @title
#' Plotting the results of Sobol SA
#'
#' @description
#' \code{plot.sobolRes_aos} plots the results of Sobol SA for objects of class 
#' \code{sobolRes_aos}.
#'
#' @param x [\code{sobolRes_aos}-object]\cr
#'   resulting output of \code{\link{ODEsobol_aos}}, of class 
#'   \code{sobolRes_aos}.
#' @param y_idx [\code{integer(1)}]\cr
#'   index of the output variable to be plotted. Defaults to 1.
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
#' @details
#' 1st order and total Sobol SA indices are plotted for one Y-variable (chosen
#' by argument \code{y_idx}) and for each input parameter against time.
#'
#' @method plot sobolRes_aos
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
#' FHNres_aos <- ODEsobol_aos(mod = FHNmod,
#'                            pars = c("a", "b", "s"),
#'                            yini = FHNyini,
#'                            times = FHNtimes,
#'                            seed = 2015,
#'                            n = 10,                      # use n >> 10!
#'                            rfuncs = c("runif", "runif", "rnorm"),
#'                            rargs = c(rep("min = 0.18, max = 0.22", 2),
#'                                  "mean = 3, sd = 0.2 / 3"),
#'                            method = "martinez",
#'                            nboot = 0)
#'
#' # Plot:
#' plot(FHNres_aos, y_idx = 1, type = "l", legendPos = "topright")
#'
#' @author Frank Weber
#' @seealso \code{\link{ODEsobol_aos}, \link[sensitivity]{sobol}, 
#' \link[sensitivity]{sobolmartinez_list}}
#'
#' @export
#' @import checkmate
#'

plot.sobolRes_aos <- function(x, y_idx = 1, type = "b", legendPos = "topleft",
                              overall_main = NULL, ...) {

  ##### Plausibilitaet #################################################
  assertClass(x, "sobolRes_aos")
  assertIntegerish(y_idx, lower = 1, upper = length(x$ST_by_y))
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
  k <- nrow((x$ST_by_y[[y_idx]])$S) - 1
  pars <- rownames((x$ST_by_y[[y_idx]])$S)[-1]
  parsCols <- rainbow(k)
  # Extrema SA Indizes:
  minMaxS <- c(0.95 * min((x$ST_by_y[[y_idx]])$S[-1, ]), 1.05 * 
                 max((x$ST_by_y[[y_idx]])$S[-1, ]))
  minMaxT <- c(0.95 * min((x$ST_by_y[[y_idx]])$T[-1, ]), 1.05 * 
                 max((x$ST_by_y[[y_idx]])$T[-1, ]))
  # Gemeinsamer Titel fuer beide Grafiken:
  if(is.null(overall_main)){
    overall_main <- paste0("Sobol sensitivity indices for y_idx = ", y_idx,
                           " and method = \"", x$method, "\"")
  }
  
  oldpar <- par(mfrow = c(1, 2), mar = c(4, 4, 1, 2) + 0.2,
                oma = c(0, 0, 2, 0))
  
  ##### 1st order SA indices ###########################################
  # Plot fuer ersten Parameter:
  plot(x = (x$ST_by_y[[y_idx]])$S[1, ], y = (x$ST_by_y[[y_idx]])$S[2, ],
       xlab = "Time", ylab = "1st order Sobol SA indices",
       type = type, col = parsCols[1], ylim = minMaxS, ...)
  # Plots fuer alle weiteren Parameter:
  for(i in 2:k) {
    lines(x = (x$ST_by_y[[y_idx]])$S[1, ], y = (x$ST_by_y[[y_idx]])$S[i + 1, ],
          type = type, col = parsCols[i], ...)
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
  plot(x = (x$ST_by_y[[y_idx]])$T[1, ], y = (x$ST_by_y[[y_idx]])$T[2, ],
       xlab = "Time", ylab = "Total Sobol SA indices",
       type = type, col = parsCols[1], ylim = minMaxT, ...)
  # Plots fuer alle weiteren Parameter:
  for(i in 2:k) {
    lines(x = (x$ST_by_y[[y_idx]])$T[1, ], y = (x$ST_by_y[[y_idx]])$T[i + 1, ],
          type = type, col = parsCols[i], ...)
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
