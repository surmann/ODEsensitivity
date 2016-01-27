#' @title
#' Plotting the results of Sobol' SA for objects of class \code{sobolRes_aos}
#'
#' @description
#' \code{plot.sobolRes_aos} plots the results of Sobol' SA for objects of class 
#' \code{sobolRes_aos}.
#'
#' @param x [\code{sobolRes_aos}-object]\cr
#'   resulting output of \code{\link{ODEsobol_aos}}, of class 
#'   \code{sobolRes_aos}.
#' @param y_plot [\code{character(1)}]\cr
#'   name of the \code{yini}-variable to be plotted. Defaults to the name of 
#'   the first \code{yini}-variable.
#' @param type [\code{character(1)}]\cr
#'   plot type, i.e. \code{"p", "l", "b", "c", "o", "s", "h"} or \code{"n"}. 
#'   Defaults to \code{"b"}.
#' @param main_title [\code{character(1)}]\cr
#'   the common title for the two graphics. Default is \code{NULL}, which means
#'   an automatic title is generated.
#' @param legendPos [\code{character(1)}]\cr
#'   legend position, default is \code{"topleft"}.
#' @param ... additional arguments passed to \code{\link{plot}}.
#'
#' @return TRUE (invisible; for testing purposes).
#'
#' @details
#' 1st order and total Sobol' SA indices are plotted for one 
#' \code{yini}-variable (chosen by argument \code{y_plot}) and for each input 
#' parameter against time.
#'
#' @note Unfortunately, the passing of arguments (e.g. "main") does not work
#'   correctly.
#'
#' @author Frank Weber
#' @seealso \code{\link{ODEsobol_aos}, \link[sensitivity]{soboljansen_list}, 
#' \link[sensitivity]{sobolmartinez_list}}
#'
#' @examples
#' ##### FitzHugh-Nagumo equations (Ramsay et al., 2007)
#' # definition of the model itself, parameters, initial values
#' # and the times vector:
#' 
#' @import checkmate
#' @method plot sobolRes_aos
#' @export
#'

plot.sobolRes_aos <- function(x, y_plot = names(x$ST_by_y)[1], type = "b",
                              main_title = NULL, 
                              legendPos = "topleft", ...) {

  ##### Plausibilitaet #################################################
  assertClass(x, "sobolRes_aos")
  assertCharacter(y_plot, len = 1)
  stopifnot(y_plot %in% names(x$ST_by_y))
  assertCharacter(type, len = 1)
  notOk <- !type %in% c("p", "l", "b", "c", "n", "o", "s", "h")
  if(notOk)
    stop(paste("type must be one of \"p\", \"l\", \"b\", \"c\", \"n\",",
               "\"o\", \"s\" or \"h\"!"))
  stopifnot(is.character(main_title) && length(main_title) == 1 ||
              is.null(main_title))
  assertCharacter(legendPos, len = 1)
  notOk <- !any(rep(legendPos, 9) == c("bottomright", "bottom", "bottomleft", 
              "left", "topleft", "top", "topright", "right", "center"))
  if(notOk)
    stop("legendPos must be one of \"bottomright\", \"bottom\",
         \"bottomleft\", \"left\", \"topleft\", \"top\", \"topright\",
         \"right\", \"center\"!")

  ##### Vorbereitungen #################################################
  
  # Index der zu plottenden yini-Variable:
  y_idx <- which(names(x$ST_by_y) == y_plot)
  # Extrahiere Anzahl Parameter:
  k <- nrow((x$ST_by_y[[y_idx]])$S) - 1
  pars <- rownames((x$ST_by_y[[y_idx]])$S)[-1]
  parsCols <- rainbow(k)
  # Extrema SA Indizes:
  minMaxS <- c(0.95 * min((x$ST_by_y[[y_idx]])$S[-1, ]), 1.05 * 
                 max((x$ST_by_y[[y_idx]])$S[-1, ]))
  minMaxT <- c(0.95 * min((x$ST_by_y[[y_idx]])$T[-1, ]), 1.05 * 
                 max((x$ST_by_y[[y_idx]])$T[-1, ]))
  # Gemeinsamer Titel fuer beide Grafiken:
  if(is.null(main_title)){
    main_title <- paste0("Sobol' sensitivity indices for State Variable \"", 
                         y_plot, "\" and method = \"", x$method, "\"")
  }
  
  oldpar <- par(mfrow = c(1, 2), mar = c(4, 4, 1, 2) + 0.2,
                oma = c(0, 0, 2, 0))
  
  ##### 1st order SA indices ###########################################
  
  # Plot fuer ersten Parameter:
  plot(x = (x$ST_by_y[[y_idx]])$S[1, ], y = (x$ST_by_y[[y_idx]])$S[2, ],
       xlab = "Time", ylab = "1st order Sobol' SA indices",
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
  
  ##### Total SA indices ###############################################
  
  # Plot fuer ersten Parameter:
  plot(x = (x$ST_by_y[[y_idx]])$T[1, ], y = (x$ST_by_y[[y_idx]])$T[2, ],
       xlab = "Time", ylab = "Total Sobol' SA indices",
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
  mtext(main_title, side = 3, line = 0, outer = TRUE, cex = 1.2, font = 2)
  
  par(oldpar)
  return(invisible(TRUE))
}
