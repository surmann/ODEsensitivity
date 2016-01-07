#' @title
#' Plot of the results of Morris SA for objects of class \code{morrisRes}
#'
#' @description
#' \code{plot.morrisRes} plots the results of Morris SA for objects of class 
#' \code{morrisRes}.
#'
#' @param x [\code{morrisRes}]\cr
#'   resulting output of \code{\link{ODEmorris}}, of class \code{morrisRes}.
#' @param y_idx [\code{integer(1)}]\cr
#'   index of the output variable to be plotted. Defaults to 1.
#' @param type [\code{character(1)}]\cr
#'   plot type, choose between \code{"sep"} and \code{"trajec"}.
#' @param colors_pars [\code{character(>= k)}]\cr
#'   vector of the colors to be used for the \code{k} different parameters 
#'   (where \code{k} is the length of the vector returned by 
#'   \code{\link[ODEnetwork]{createParamVecs}} applied to the \code{ODEnetwork}-
#'   object for which the sensitivity analysis was done). Must be at least of
#'   length \code{k}.
#' @param main_title [\code{character(1)}]\cr
#'   title for the plot. If \code{type = "sep"}, this is the overall title for
#'   the two separate plots. Defaults to NULL, so a standard title is generated.
#' @param legendPos [\code{character(1)}]\cr
#'   legend position, default is \code{"topleft"}.
#' @param ... additional arguments passed to \code{\link{plot}}.
#'
#' @return TRUE (invisible; for testing purposes).
#'
#' @details
#' \code{plot} with \code{type = "sep"} plots mu.star and
#'   sigma separately versus time.
#'
#' \code{plot} with \code{type = "trajec"} plots mu.star versus
#'   sigma for every point of time.
#'
#' @note Unfortunately, the passing of arguments (e.g. "main") does not work
#'   correctly.
#'
#' @author Frank Weber
#' @seealso \code{\link{ODEmorris}},
#'   \code{\link[sensitivity]{morris_list}}
#'
#' @examples
#' library(ODEnetwork)
#' masses <- c(1, 1)
#' dampers <- diag(c(1, 1))
#' springs <- diag(c(1, 1))
#' springs[1, 2] <- 1
#' distances <- diag(c(0, 2))
#' distances[1, 2] <- 1
#' odenet <- ODEnetwork(masses, dampers, springs, 
#'                      cartesian = TRUE, distances = distances)
#' odenet <- setState(odenet, c(0.5, 1), c(0, 0))
#' 
#' ODEtimes <- seq(0.01, 20, by = 0.1)
#' ODEbinf <- c(rep(0.001, 9), -4, 0.001, -10)
#' ODEbsup <- c(2, 1.5, 6, 4, 6, 10, 2, 1.5, 6, -0.001, 6, -0.001)
#' 
#' ODEres <- ODEmorris(odenet, ODEtimes, ode_method = "adams", seed = 2015, 
#'                     binf = ODEbinf, bsup = ODEbsup, r = 20)
#' 
#' # Standard (separate plots):
#' plot(ODEres, y_idx = 1, type = "sep", legendPos = "topleft")
#' plot(ODEres, y_idx = 1, type = "sep", legendPos = "outside")
#' # Palette "Dark2" from the package "RColorBrewer" with some 
#' # additional colors:
#' my_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
#'              "#E6AB02", "#A6761D", "#666666", "black", "firebrick",
#'              "darkblue", "darkgreen")
#' plot(ODEres, y_idx = 1, type = "sep", colors_pars = my_cols, 
#'      legendPos = "outside")
#' 
#' # Trajectories:
#' plot(ODEres, y_idx = 1, type = "trajec", legendPos = "topleft")
#' plot(ODEres, y_idx = 1, type = "trajec", legendPos = "outside")
#' plot(ODEres, y_idx = 1, type = "trajec", colors_pars = my_cols, 
#'      legendPos = "outside")
#'
#' @import
#'   checkmate
#' @method plot morrisRes
#' @export
#'

plot.morrisRes <- function(x, y_idx = 1, type = "sep", colors_pars = NULL,
                           main_title = NULL, legendPos = "outside", ...) {

  ##### Check input #################################################
  assertClass(x, "morrisRes")
  assertIntegerish(y_idx, lower = 1, upper = length(x))
  assertCharacter(type, len = 1)
  notOk <- !type %in% c("sep", "trajec")
  if(notOk)
    stop("type must be one of \"sep\" or \"trajec\"!")
  stopifnot((is.character(colors_pars) && 
              length(colors_pars) >= (nrow(x[[y_idx]]) - 1) / 3) || 
              is.null(colors_pars))
  assertCharacter(legendPos, len = 1)
  notOk <- !legendPos %in% c("outside", "bottomright", "bottom",
    "bottomleft", "left", "topleft", "top", "topright", "right",
    "center")
  if(notOk)
    stop("legendPos must be one of \"outside\", \"bottomright\", \"bottom\",
      \"bottomleft\", \"left\", \"topleft\", \"top\", \"topright\",
      \"right\", \"center\"!")

  ##### Plot ###########################################################
  
  # Extrahiere die Parameter-Namen:
  k <- (nrow(x[[y_idx]]) - 1) / 3
  pars_tmp <- rownames(x[[y_idx]])[2:(k + 1)]
  pars <- substr(pars_tmp, start = 4, stop = nchar(pars_tmp))
  state_name <- names(x)[y_idx]
  
  # Separate Plots fuer mu.star und sigma:
  if(type == "sep"){
    plotSep(x[[y_idx]], pars, colors_pars, 
            state_name, common_title = main_title, legendPos, ...)
  }
  
  # Trajectories:
  if(type == "trajec"){
    plotTrajectories(x[[y_idx]], pars, colors_pars, 
                     state_name, main_title, legendPos, ...)
  }
  
  # For testing purposes:
  return(invisible(TRUE))
}
