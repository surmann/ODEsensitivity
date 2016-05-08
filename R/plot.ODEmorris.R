#' @title
#' Plot of the Results of Morris Screening for Objects of Class \code{ODEmorris}
#'
#' @description
#'   \code{plot.ODEmorris} plots the results of Morris screening for objects of 
#'   class \code{ODEmorris}.
#'
#' @param x [\code{ODEmorris}]\cr
#'   output of \code{\link{ODEmorris}} (of class \code{ODEmorris}).
#' @param pars_plot [\code{character(k)}]\cr
#'   names of the \code{k} parameters to be plotted. If \code{NULL} (the 
#'   default), all parameters are plotted.
#' @param state_plot [\code{character(1)}]\cr
#'   name of the state variable to be plotted. Defaults to the name of the
#'   first state variable.
#' @param kind [\code{character(1)}]\cr
#'   kind of the plot, choose between \code{"sep"} and \code{"trajec"}.
#' @param colors_pars [\code{character(>= k)}]\cr
#'   vector of the colors to be used for the \code{k} different parameters. Must
#'   be at least of length \code{k} (only the first \code{k} elements will be
#'   used, though). If \code{NULL} (the default), \code{rainbow(k)} is used.
#' @param main_title [\code{character(1)}]\cr
#'   title for the plot. If \code{kind = "sep"}, this is the overall title for
#'   the two separate plots. If \code{NULL} (the default), a standard title is 
#'   generated.
#' @param legendPos [\code{character(1)}]\cr
#'   keyword for the legend position, either one of those specified in
#'   \code{\link{legend}} or \code{"outside"} (the default), which means the 
#'   legend is placed under the plot (useful, if there are many parameters in 
#'   the model).
#' @param type [\code{character(1)}]\cr
#'   plot type, i.e. \code{"p", "l", "b", "c", "o", "s", "h"} or \code{"n"}. 
#'   Defaults to \code{"l"}.
#' @param ... additional arguments passed to \code{\link{plot.default}}.
#'
#' @return \code{TRUE} (invisible; for testing purposes).
#'
#' @details
#'   \code{plot.ODEmorris} with \code{kind = "sep"} plots mu.star and
#'   sigma separately versus time.
#'
#'   \code{plot.ODEmorris} with \code{kind = "trajec"} plots mu.star versus
#'   sigma for every point of time.
#'
#' @note 
#'   Not all plotting arguments can be passed by \code{...}, for example
#'   \code{xlab} and \code{ylab} are fixed.
#'
#' @author Stefan Theers, Frank Weber
#' @seealso \code{\link{ODEmorris}, \link[sensitivity]{morris}}
#'
#' @examples
#' ##### FitzHugh-Nagumo equations (Ramsay et al., 2007) #####
#' # Definition of the model itself, parameters, initial state values
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
#' FHNstate  <- c(Voltage = -1, Current = 1)
#' FHNtimes <- seq(0.1, 50, by = 5)
#' 
#' set.seed(4628)
#' FHNres <- ODEmorris(mod = FHNmod,
#'                     pars = c("a", "b", "s"),
#'                     state_init = FHNstate,
#'                     times = FHNtimes,
#'                     binf = c(0.18, 0.18, 2.8),
#'                     bsup = c(0.22, 0.22, 3.2),
#'                     r = 50,
#'                     design = list(type = "oat", levels = 100, grid.jump = 1),
#'                     scale = TRUE,
#'                     ode_method = "adams",
#'                     parallel_eval = TRUE,
#'                     parallel_eval_ncores = 2)
#' 
#' # Define custom colors for the plot:
#' my_cols <- c("firebrick", "chartreuse3", "dodgerblue")
#' plot(FHNres, state_plot = "Current", kind = "sep", colors_pars = my_cols)
#' plot(FHNres, state_plot = "Current", kind = "trajec", colors_pars = my_cols)
#' 
#' ##### A network of ordinary differential equations #####
#' # Definition of the network using the package "ODEnetwork":
#' library(ODEnetwork)
#' masses <- c(1, 1)
#' dampers <- diag(c(1, 1))
#' springs <- diag(c(1, 1))
#' springs[1, 2] <- 1
#' distances <- diag(c(0, 2))
#' distances[1, 2] <- 1
#' lfonet <- ODEnetwork(masses, dampers, springs, 
#'                      cartesian = TRUE, distances = distances)
#' lfonet <- setState(lfonet, c(0.5, 1), c(0, 0))
#' LFOpars <- c("k.1", "k.2", "k.1.2")
#' LFOtimes <- seq(0.01, 20, by = 0.1)
#' LFObinf <- rep(0.001, 3)
#' LFObsup <- c(6, 6, 3)
#' 
#' set.seed(4628)
#' LFOres <- ODEmorris(lfonet,
#'                     LFOpars,
#'                     LFOtimes,
#'                     binf = LFObinf,
#'                     bsup = LFObsup,
#'                     r = 50,
#'                     design = list(type = "oat", levels = 100, grid.jump = 1),
#'                     scale = TRUE,
#'                     ode_method = "adams",
#'                     parallel_eval = TRUE,
#'                     parallel_eval_ncores = 2)
#' 
#' plot(LFOres, state_plot = "x.2", kind = "sep", colors_pars = my_cols)
#' plot(LFOres, state_plot = "x.2", kind = "trajec", colors_pars = my_cols)
#'
#' @import checkmate
#' @method plot ODEmorris
#' @export
#'

plot.ODEmorris <- function(x, pars_plot = NULL, state_plot = names(x)[1], 
                           kind = "sep", colors_pars = NULL, main_title = NULL, 
                           legendPos = "outside", type = "l", ...) {
  
  ##### Input checks ###################################################
  
  assertClass(x, "ODEmorris")
  stopifnot(is.null(pars_plot) || is.character(pars_plot))
  assertCharacter(state_plot, len = 1)
  stopifnot(state_plot %in% names(x))
  assertCharacter(kind, len = 1)
  if(!kind %in% c("sep", "trajec")){
    stop("kind must be one of \"sep\" or \"trajec\"!")
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
  assertCharacter(type, len = 1)
  if(!type %in% c("p", "l", "b", "c", "n", "o", "s", "h")){
    stop(paste("type must be one of \"p\", \"l\", \"b\", \"c\", \"n\",",
               "\"o\", \"s\" or \"h\"!"))
  }
  
  ##### Preparation ####################################################
  
  # Index of the state variable to be plotted:
  state_idx <- which(names(x) == state_plot)
  
  # Extract the parameter names, set "pars_plot" if not specified and perform
  # a validity check if parameter names are user-specified:
  k <- (nrow(x[[state_idx]]) - 1) / 3
  pars_tmp <- rownames(x[[state_idx]])[2:(k + 1)]
  all_pars <- substr(pars_tmp, start = 4, stop = nchar(pars_tmp))
  if(is.null(pars_plot)){
    pars_plot <- all_pars
  } else{
    stopifnot(all(pars_plot %in% all_pars))
    if(any(duplicated(pars_plot))){
      pars_plot <- unique(pars_plot)
    }
  }
  
  ##### Plot ###########################################################
  
  # Separate plots for mu.star und sigma:
  if(kind == "sep"){
    plotSep(x[[state_idx]], pars = pars_plot, state_name = state_plot, 
            colors_pars = colors_pars, common_title = main_title, 
            legendPos = legendPos, type = type, ...)
  }
  
  # Trajectories:
  if(kind == "trajec"){
    plotTrajectories(x[[state_idx]], pars = pars_plot, state_name = state_plot, 
                     colors_pars = colors_pars, main_title = main_title, 
                     legendPos = legendPos, type = type, ...)
  }
  
  # For testing purposes:
  return(invisible(TRUE))
}
