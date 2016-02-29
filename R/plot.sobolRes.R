#' @title
#' Plotting the results of Sobol' SA for objects of class \code{sobolRes}
#'
#' @description
#'   \code{plot.sobolRes} plots the results of Sobol' SA for objects of class 
#'   \code{sobolRes}.
#'
#' @param x [\code{sobolRes}]\cr
#'   resulting output of \code{\link{ODEsobol}}, of class \code{sobolRes}.
#' @param state_plot [\code{character(1)}]\cr
#'   name of the state variable to be plotted. Defaults to the name of the
#'   first state variable.
#' @param colors_pars [\code{character(>= k)}]\cr
#'   vector of the colors to be used for the \code{k} different parameters. Must
#'   be at least of length \code{k}. If \code{NULL} (the default), 
#'   \code{rainbow(k)} is used.
#' @param main_title [\code{character(1)}]\cr
#'   common title for the two graphics. Default is \code{NULL}, which means
#'   an automatic title is generated.
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
#' @return TRUE (invisible; for testing purposes).
#'
#' @details
#'   First order and total Sobol' SA indices are plotted for one state variable 
#'   (chosen by argument \code{state_plot}) and all of the parameters against 
#'   time.
#'
#' @note 
#'   Not all arguments of \code{\link{plot.default}} can be passed by 
#'   \code{...}, for example \code{xlab} and \code{ylab} are fixed.
#'
#' @author Frank Weber
#' @seealso \code{\link{ODEsobol}, \link[sensitivity]{soboljansen},
#' \link[sensitivity]{sobolmartinez}}
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
#' # Warning: The following code might take a long time!
#' FHNres <- ODEsobol(mod = FHNmod,
#'                    pars = c("a", "b", "s"),
#'                    state_init = FHNstate,
#'                    times = FHNtimes,
#'                    seed = 2015,
#'                    n = 1000,
#'                    rfuncs = c("runif", "rnorm", "rexp"),
#'                    rargs = c("min = 0.18, max = 0.22", 
#'                              "mean = 0.2, sd = 0.2 / 3",
#'                              "rate = 1 / 3"),
#'                    sobol_method = "martinez",
#'                    ode_method = "adams",
#'                    ode_parallel = TRUE,
#'                    ode_parallel_ncores = 2)
#' 
#' # Define custom colors for the plot:
#' my_cols <- c("firebrick", "chartreuse3", "dodgerblue")
#' plot(FHNres, state_plot = "Current", colors_pars = my_cols)
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
#' 
#' LFOres <- ODEsobol(lfonet, 
#'                    LFOpars, 
#'                    LFOtimes, 
#'                    seed = 2015, 
#'                    n = 1000,
#'                    rfuncs = c("runif", "rnorm", "rexp"),
#'                    rargs = c("min = 0.001, max = 6",
#'                              "mean = 3, sd = 0.5",
#'                              "rate = 1 / 3"),
#'                    sobol_method = "martinez",
#'                    ode_method = "adams",
#'                    ode_parallel = TRUE,
#'                    ode_parallel_ncores = 2)
#' # (A warning is thrown, concerning the state variables "v.1" and "v.2". 
#' # Assuming that we are only interested in a sensitivity analysis of "x.1" and
#' # "x.2", this warning is ignored.)
#' 
#' plot(LFOres, state_plot = "x.1", colors_pars = my_cols)
#' plot(LFOres, state_plot = "x.2", colors_pars = my_cols)
#' 
#' @import checkmate
#' @method plot sobolRes
#' @export
#'

plot.sobolRes <- function(x, state_plot = names(x$ST_by_state)[1],
                          colors_pars = NULL, main_title = NULL,
                          legendPos = "outside", type = "l", ...) {
  
  ##### Input checks ###################################################
  
  assertClass(x, "sobolRes")
  assertCharacter(state_plot, len = 1)
  stopifnot(state_plot %in% names(x$ST_by_state))
  assertCharacter(type, len = 1)
  if(!type %in% c("p", "l", "b", "c", "n", "o", "s", "h")){
    stop(paste("type must be one of \"p\", \"l\", \"b\", \"c\", \"n\",",
               "\"o\", \"s\" or \"h\"!"))
  }
  stopifnot(is.null(colors_pars) || (is.character(colors_pars) && 
    length(colors_pars) >= 
      nrow((x$ST_by_state[[which(names(x$ST_by_state) == state_plot)]])$S) - 1))
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
  
  # Index of the state variable to be plotted:
  state_idx <- which(names(x$ST_by_state) == state_plot)
  
  # Extract the matrices S and T:
  S <- (x$ST_by_state[[state_idx]])$S
  T <- (x$ST_by_state[[state_idx]])$T
  
  # Extract the parameter names:
  k <- nrow(S) - 1
  pars <- rownames(S)[-1]
  
  # Extreme values of the SA indices:
  minMaxS <- c(0.95 * min(S[-1, ]), 1.05 * 
                 max(S[-1, ]))
  minMaxT <- c(0.95 * min(T[-1, ]), 1.05 * 
                 max(T[-1, ]))
  
  # Set colors if not set by the user:
  if(is.null(colors_pars)){
    colors_pars <- rainbow(k)
  }
  
  # Set main title if not set by the user:
  if(x$sobol_method == "jansen"){
    method_title <- "Jansen"
  } else if(x$sobol_method == "martinez"){
    method_title <- "Martinez"
  }
  if(is.null(main_title)){
    main_title <- paste0("Sobol' sensitivity indices for \"", state_plot, 
                         "\" (Sobol'-", method_title, " method)")
  }
  
  ##### Plot ###########################################################
  
  if(legendPos == "outside"){
    if(k > 6){
      # Multirow legend:
      legend_ncol <- ceiling(k / 2)
      oldpar <- par(mfrow = c(1, 2),
                    oma = c(2.5, 0, 2, 0), mar = c(4, 4, 1, 2) + 0.2)
    } else{
      # One-row legend:
      legend_ncol <- k
      oldpar <- par(mfrow = c(1, 2),
                    oma = c(1.7, 0, 2, 0), mar = c(4, 4, 1, 2) + 0.2)
    }
  } else{
    oldpar <- par(mfrow = c(1, 2),
                  oma = c(0, 0, 2, 0), mar = c(4, 4, 1, 2) + 0.2)
  }
  
  # First order SA indices:
  # First parameter:
  plot(x = S[1, ], y = S[2, ],
       xlab = "Time", ylab = "First order Sobol' SA indices",
       type = type, col = colors_pars[1], ylim = minMaxS, ...)
  # All remaining parameters:
  for(i in 2:k) {
    lines(x = S[1, ], 
          y = S[i + 1, ],
          type = type, col = colors_pars[i], ...)
  }
  # Legend:
  if(legendPos != "outside"){
    if(type %in% c("b", "o")){
      legend(legendPos, legend = pars, col = colors_pars, bg = "white",
             lty = 1, pch = 1)
    } else if(type %in% c("l", "c", "s", "h")){
      legend(legendPos, legend = pars, col = colors_pars, bg = "white", lty = 1)
    } else if(type == "p"){
      legend(legendPos, legend = pars, col = colors_pars, bg = "white", pch = 1)
    }
  }
  
  # Total SA indices:
  # First parameter:
  plot(x = T[1, ], y = T[2, ],
       xlab = "Time", ylab = "Total Sobol' SA indices",
       type = type, col = colors_pars[1], ylim = minMaxT, ...)
  # All remaining parameters:
  for(i in 2:k) {
    lines(x = T[1, ], 
          y = T[i + 1, ],
          type = type, col = colors_pars[i], ...)
  }
  # Legend:
  if(legendPos != "outside"){
    if(type %in% c("b", "o")){
      legend(legendPos, legend = pars, col = colors_pars, bg = "white",
             lty = 1, pch = 1)
    } else if(type %in% c("l", "c", "s", "h")){
      legend(legendPos, legend = pars, col = colors_pars, bg = "white", lty = 1)
    } else if(type == "p"){
      legend(legendPos, legend = pars, col = colors_pars, bg = "white", pch = 1)
    }
  }
  
  # Create the big common title for the two plots:
  mtext(main_title, side = 3, line = 0, outer = TRUE, cex = 1.2, font = 2)
  
  # Legend outside of the plotting region:
  if(legendPos == "outside"){
    # Dummy plot:
    oldpar2 <- par(mfrow = c(1, 1), oma = rep(0, 4), mar = rep(0, 4), 
                   new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    # Draw the legend depending on the plot type:
    if(type %in% c("b", "o")){
      legend("bottom", legend = pars, col = colors_pars, lty = 1, pch = 1, 
             bty = "n", xpd = TRUE, ncol = legend_ncol, inset = c(0, 0))
    } else if(type %in% c("l", "c", "s", "h")){
      legend("bottom", legend = pars, col = colors_pars, lty = 1, 
             bty = "n", xpd = TRUE, ncol = legend_ncol, inset = c(0, 0))
    } else if(type == "p"){
      legend("bottom", legend = pars, col = colors_pars, pch = 1, 
             bty = "n", xpd = TRUE, ncol = legend_ncol, inset = c(0, 0))
    }
    par(oldpar2)
  }
  
  par(oldpar)
  # For testing purposes:
  return(invisible(TRUE))
}
