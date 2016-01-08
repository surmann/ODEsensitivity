#' @title
#' Plotting the results of Morris SA for objects of class \code{morrisRes_aos}
#'
#' @description
#' \code{plot.morrisRes_aos} plots the results of Morris SA for objects of class 
#' \code{morrisRes_aos}.
#'
#' @details
#' \code{plot} with \code{type = "sep"} plots mu.star and
#'   sigma separately versus time.
#'
#' \code{plot} with \code{type = "trajec"} plots mu.star versus
#'   sigma for every point of time.
#'
#' @param x [\code{morrisRes_aos}]\cr
#'   resulting output of \code{\link{ODEmorris_aos}}, of class 
#'   \code{morrisRes_aos}.
#' @param y_plot [\code{character(1)}]\cr
#'   name of the \code{yini}-variable to be plotted. Defaults to the name of 
#'   the first \code{yini}-variable.
#' @param type [\code{character(1)}]\cr
#'   plot type, choose between \code{"sep"} and \code{"trajec"}.
#' @param colors_pars [\code{character(>= k)}]\cr
#'   vector of the colors to be used for the \code{k} different parameters 
#'   (where \code{k} is the of parameters in the model for which SA was done). 
#'   Must be at least of length \code{k}. If \code{NULL} (the default), 
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
#' @method plot morrisRes_aos
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
#' FHNyini  <- c(Voltage = -1, Current = 1)
#' FHNtimes <- seq(0.1, 50, by = 5)
#'
#' FHNres_aos <- ODEmorris_aos(mod = FHNmod,
#'                             pars = c("a", "b", "s"),
#'                             yini = FHNyini,
#'                             times = FHNtimes,
#'                             ode_method = "adams",
#'                             seed = 2015,
#'                             binf = c(0.18, 0.18, 2.8),
#'                             bsup = c(0.22, 0.22, 3.2),
#'                             r = 10,
#'                             design =
#'                                list(type = "oat", levels = 30, 
#'                                     grid.jump = 1),
#'                             scale = TRUE)
#' 
#' # Palette "Dark2" from the package "RColorBrewer" with some 
#' # additional colors:
#' my_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
#'              "#E6AB02", "#A6761D", "#666666", "black", "firebrick",
#'              "darkblue", "darkgreen")
#' # Voltage:
#' plot(FHNres_aos, y_plot = "Voltage", type = "sep")
#' plot(FHNres_aos, y_plot = "Voltage", type = "sep", legendPos = "topleft")
#' plot(FHNres_aos, y_plot = "Voltage", type = "sep", colors_pars = my_cols)
#' plot(FHNres_aos, y_plot = "Voltage", type = "trajec")
#' plot(FHNres_aos, y_plot = "Voltage", type = "trajec", legendPos = "topleft")
#' plot(FHNres_aos, y_plot = "Voltage", type = "trajec", colors_pars = my_cols)
#' # Current:
#' plot(FHNres_aos, y_plot = "Current", type = "sep")
#' plot(FHNres_aos, y_plot = "Current", type = "sep", legendPos = "topleft")
#' plot(FHNres_aos, y_plot = "Current", type = "sep", colors_pars = my_cols)
#' plot(FHNres_aos, y_plot = "Current", type = "trajec")
#' plot(FHNres_aos, y_plot = "Current", type = "trajec", legendPos = "topleft")
#' plot(FHNres_aos, y_plot = "Current", type = "trajec", colors_pars = my_cols)
#'
#' @note Unfortunately, the passing of arguments (e.g. "main") does not work
#'   correctly.
#'
#' @author Frank Weber
#' @seealso \code{\link{ODEmorris_aos}},
#'   \code{\link[sensitivity]{morris_list}}
#'
#' @export
#' @import
#'   checkmate
#'

plot.morrisRes_aos <- function(x, y_plot = names(x)[1], type = "sep", 
                               colors_pars = NULL, main_title = NULL, 
                               legendPos = "outside", ...) {
  
  ##### Check input #################################################
  assertClass(x, "morrisRes_aos")
  assertCharacter(y_plot, len = 1)
  stopifnot(y_plot %in% names(x))
  assertCharacter(type, len = 1)
  notOk <- !type %in% c("sep", "trajec")
  if(notOk)
    stop("type must be one of \"sep\" or \"trajec\"!")
  stopifnot(is.null(colors_pars) || (is.character(colors_pars) && 
            length(colors_pars) >= (nrow(x[[which(names(x) == y_plot)]]) - 1) / 
              3))
  stopifnot(is.null(main_title) || (is.character(main_title) && 
            length(main_title) == 1))
  assertCharacter(legendPos, len = 1)
  notOk <- !legendPos %in% c("outside", "bottomright", "bottom",
    "bottomleft", "left", "topleft", "top", "topright", "right", "center")
  if(notOk)
    stop("legendPos must be one of \"outside\", \"bottomright\", \"bottom\",
         \"bottomleft\", \"left\", \"topleft\", \"top\", \"topright\",
         \"right\", \"center\"!")

  ##### Plot ###########################################################
  
  # Index der zu plottenden yini-Variable:
  y_idx <- which(names(x) == y_plot)
  # Extrahiere die Parameter-Namen:
  k <- (nrow(x[[y_idx]]) - 1) / 3
  pars_tmp <- rownames(x[[y_idx]])[2:(k + 1)]
  pars <- substr(pars_tmp, start = 4, stop = nchar(pars_tmp))
  
  # Separate Plots fuer mu.star und sigma:
  if(type == "sep"){
    plotSep(x[[y_idx]], pars, y_name = y_plot, colors_pars,
            common_title = main_title, legendPos = legendPos, ...)
  }
  
  # Trajectories:
  if(type == "trajec"){
    plotTrajectories(x[[y_idx]], pars, y_name = y_plot, 
                     colors_pars = colors_pars, 
                     main_title = main_title, legendPos = legendPos, ...)
  }
  
  # For testing purposes:
  return(invisible(TRUE))
}
