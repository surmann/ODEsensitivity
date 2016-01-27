#' @title
#' Plotting the results of Morris SA for objects of class \code{morrisRes_ats}
#'
#' @description
#' \code{plot.morrisRes_ats} plots the results of Morris SA for objects of class 
#' \code{morrisRes_ats}.
#'
#' @param x [\code{morrisRes_ats}]\cr
#'   resulting output of \code{\link{ODEmorris_ats}}, of class 
#'   \code{morrisRes_ats}.
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
#' @details
#' \code{plot} with \code{type = "sep"} plots mu.star and
#'   sigma separately versus time.
#'
#' \code{plot} with \code{type = "trajec"} plots mu.star versus
#'   sigma for every point of time.
#'
#' @examples
#' ##### FitzHugh-Nagumo equations (Ramsay et al., 2007)
#' # definition of the model itself, parameters, initial values
#' # and the times vector:
#' FHNmod <- function(Time, State, Pars) {
#' with(as.list(c(State, Pars)), {
#'   
#'   dVoltage <- s * (Voltage - Voltage^3 / 3 + Current)
#'   dCurrent <- - 1 / s *(Voltage - a + b * Current)
#'   
#'   return(list(c(dVoltage, dCurrent)))
#' })
#' }
#' 
#' FHNyini  <- c(Voltage = -1, Current = 1)
#' FHNtimes <- seq(0.1, 50, by = 5)
#' FHNres_ats <- ODEmorris_ats(mod = FHNmod,
#'                             pars = c("a", "b", "s"),
#'                             yini = FHNyini,
#'                             times = FHNtimes,
#'                             y_analyzed = "Voltage",
#'                             ode_method = "adams",
#'                             seed = 2015,
#'                             binf = c(0.18, 0.18, 2.8),
#'                             bsup = c(0.22, 0.22, 3.2),
#'                             r = 10,
#'                             design =
#'                               list(type = "oat", levels = 30, 
#'                                    grid.jump = 1),
#'                             scale = TRUE)
#' 
#' # Palette "Dark2" from the package "RColorBrewer" with some 
#' # additional colors:
#' my_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
#'              "#E6AB02", "#A6761D", "#666666", "black", "firebrick",
#'              "darkblue", "darkgreen")
#' plot(FHNres_ats, type = "sep")
#' plot(FHNres_ats, type = "sep", legendPos = "topleft")
#' plot(FHNres_ats, type = "sep", colors_pars = my_cols)
#' plot(FHNres_ats, type = "trajec")
#' plot(FHNres_ats, type = "trajec", legendPos = "topleft")
#' plot(FHNres_ats, type = "trajec", colors_pars = my_cols)
#'
#' @note Unfortunately, the passing of arguments (e.g. "main") does not work
#'   correctly.
#'
#' @author Frank Weber
#' @seealso \code{\link{ODEmorris_ats}},
#'   \code{\link[sensitivity]{morris_matrix}}
#' 
#' @import
#'   checkmate
#' @method plot morrisRes_ats
#' @export
#'

plot.morrisRes_ats <- function(x, type = "sep", colors_pars = NULL, 
                               main_title = NULL, legendPos = "outside", ...) {

  ##### Check input #################################################
  assertClass(x, "morrisRes_ats")
  assertCharacter(type, len = 1)
  notOk <- !type %in% c("sep", "trajec")
  if(notOk)
    stop("type must be one of \"sep\" or \"trajec\"!")
  stopifnot(is.null(colors_pars) || (is.character(colors_pars) && 
            length(colors_pars) >= (nrow(x$res) - 1) / 3))
  stopifnot(is.null(main_title) || (is.character(main_title) && 
            length(main_title) == 1))
  assertCharacter(legendPos, len = 1)
  notOk <- !legendPos %in% c("outside", "bottomright", "bottom",
    "bottomleft", "left", "topleft", "top", "topright", "right", "center")
  if(notOk)
    stop("legendPos must be one of \"outside\", \"bottomright\", \"bottom\",
         \"bottomleft\", \"left\", \"topleft\", \"top\", \"topright\",
         \"right\", \"center\"!")

  ##### Plots ###########################################################
  
  # Extrahiere die Parameter-Namen:
  k <- (nrow(x$res) - 1) / 3
  pars_tmp <- rownames(x$res)[2:(k + 1)]
  pars <- substr(pars_tmp, start = 4, stop = nchar(pars_tmp))
  
  # Separate Plots fuer mu.star und sigma:
  if(type == "sep"){
    plotSep(x$res, pars, y_name = x$y_analyzed, colors_pars = colors_pars,
            common_title = main_title, legendPos = legendPos, ...)
  }
  
  # Trajectories:
  if(type == "trajec"){
    plotTrajectories(x$res, pars, y_name = x$y_analyzed, 
                     colors_pars = colors_pars, 
                     main_title = main_title, legendPos = legendPos, ...)
  }
  
  # For testing purposes:
  return(invisible(TRUE))
}
