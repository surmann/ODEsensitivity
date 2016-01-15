#' @title
#' Plot of the results of Morris SA for objects of class \code{morrisRes}
#'
#' @description
#' \code{plot.morrisRes} plots the results of Morris SA for objects of class 
#' \code{morrisRes}.
#'
#' @param x [\code{morrisRes}]\cr
#'   resulting output of \code{\link{ODEmorris}}, of class \code{morrisRes}.
#' @param state_plot [\code{character(1)}]\cr
#'   name of the state variable to be plotted. Defaults to the name of the
#'   first state variable.
#' @param type [\code{character(1)}]\cr
#'   plot type, choose between \code{"sep"} and \code{"trajec"}.
#' @param colors_pars [\code{character(>= k)}]\cr
#'   vector of the colors to be used for the \code{k} different parameters. Must
#'   be at least of length \code{k}. If \code{NULL} (the default), 
#'   \code{rainbow(k)} is used.
#' @param main_title [\code{character(1)}]\cr
#'   title for the plot. If \code{type = "sep"}, this is the overall title for
#'   the two separate plots. Defaults to \code{NULL}, so a standard title is 
#'   generated.
#' @param legendPos [\code{character(1)}]\cr
#'   keyword for the legend position, either one of those specified in
#'   \code{\link{legend}} or \code{"outside"} (the default), which means the 
#'   legend is placed under the plot (useful, if there are many parameters in 
#'   the model).
#' @param ... additional arguments passed to \code{\link{plot}}.
#'
#' @return \code{TRUE} (invisible; for testing purposes).
#'
#' @details
#' \code{plot} with \code{type = "sep"} plots mu.star and
#'   sigma separately versus time.
#'
#' \code{plot} with \code{type = "trajec"} plots mu.star versus
#'   sigma for every point of time.
#'
#' @note Not all plotting arguments can be passed by \code{...}, for example
#' \code{xlab} and \code{ylab} are fixed.
#'
#' @author Frank Weber
#' @seealso \code{\link{ODEmorris}, \link[sensitivity]{morris_list}}
#'
#' @examples
#' library(ODEnetwork)
#' 
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
#' ODEpars <- c("m.1", "d.1", "k.1", "k.1.2", "m.2", "d.2", "k.2")
#' ODEtimes <- seq(0.01, 20, by = 0.1)
#' ODEbinf <- rep(0.001, length(ODEpars))
#' ODEbsup <- c(2, 1.5, 6, 6, 2, 1.5, 6)
#' 
#' ODEres <- ODEmorris(odenet, ODEpars, ODEtimes, ode_method = "adams", 
#'                     seed = 2015, binf = ODEbinf, bsup = ODEbsup, r = 20)
#' 
#' # Plots for state variable "x.2":
#' # Palette "Dark2" from the package "RColorBrewer" with some 
#' # additional colors:
#' my_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
#'              "#E6AB02", "#A6761D", "#666666", "black", "firebrick",
#'              "darkblue", "darkgreen")
#' # Standard (separate plots):
#' plot(ODEres, state_plot = "x.2", type = "sep", colors_pars = my_cols, 
#'      legendPos = "outside")
#' 
#' # Trajectories:
#' plot(ODEres, state_plot = "x.2", type = "trajec", colors_pars = my_cols, 
#'      legendPos = "outside")
#'
#' @import checkmate
#' @method plot morrisRes
#' @export
#'

plot.morrisRes <- function(x, state_plot = names(x)[1], type = "sep", 
                           colors_pars = NULL, main_title = NULL, 
                           legendPos = "outside", ...) {
  
  ##### Input checks ###################################################
  
  assertClass(x, "morrisRes")
  assertCharacter(state_plot, len = 1)
  stopifnot(state_plot %in% names(x))
  assertCharacter(type, len = 1)
  if(!type %in% c("sep", "trajec")){
    stop("type must be one of \"sep\" or \"trajec\"!")
  }
  stopifnot(is.null(colors_pars) || (is.character(colors_pars) && 
    length(colors_pars) >= (nrow(x[[which(names(x) == state_plot)]]) - 1) / 3))
  stopifnot(is.null(main_title) || (is.character(main_title) && 
    length(main_title) == 1))
  assertCharacter(legendPos, len = 1)
  if(!legendPos %in% c("outside", "bottomright", "bottom", "bottomleft", 
                       "left", "topleft", "top", "topright", "right", 
                       "center")){
    stop("legendPos must be one of \"outside\", \"bottomright\", \"bottom\",
         \"bottomleft\", \"left\", \"topleft\", \"top\", \"topright\",
         \"right\", \"center\"!")
  }
  
  ##### Preparation ####################################################
  
  # Index of the state-Variable which shall be plotted:
  y_idx <- which(names(x) == state_plot)
  # Extract the parameter names:
  k <- (nrow(x[[y_idx]]) - 1) / 3
  pars_tmp <- rownames(x[[y_idx]])[2:(k + 1)]
  pars <- substr(pars_tmp, start = 4, stop = nchar(pars_tmp))
  
  ##### Plot ###########################################################
  
  # Separate plots for mu.star und sigma:
  if(type == "sep"){
    plotSep(x[[y_idx]], pars, y_name = state_plot, colors_pars = colors_pars, 
            common_title = main_title, legendPos = legendPos, ...)
  }
  
  # Trajectories:
  if(type == "trajec"){
    plotTrajectories(x[[y_idx]], pars, y_name = state_plot, 
                     colors_pars = colors_pars,
                     main_title, legendPos = legendPos, ...)
  }
  
  # For testing purposes:
  return(invisible(TRUE))
}
