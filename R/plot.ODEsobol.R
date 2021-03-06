#' @title
#' Plot of the Results of Sobol' Sensitivity Analysis for Objects of Class 
#' \code{ODEsobol}
#'
#' @description
#'   \code{plot.ODEsobol} plots the results of Sobol' SA for objects of class 
#'   \code{ODEsobol}.
#'
#' @param x [\code{ODEsobol}]\cr
#'   output of \code{\link{ODEsobol}} (of class \code{ODEsobol}).
#' @param pars_plot [\code{character(k)}]\cr
#'   names of the \code{k} parameters to be plotted. If \code{NULL} (the 
#'   default), all parameters are plotted.
#' @param state_plot [\code{character(1)}]\cr
#'   name of the state variable to be plotted. Defaults to the name of the
#'   first state variable.
#' @param colors_pars [\code{character(>= k)}]\cr
#'   vector of the colors to be used for the \code{k} different parameters. Must
#'   be at least of length \code{k} (only the first \code{k} elements will be
#'   used, though). If \code{NULL} (the default), \code{rainbow(k)} is used.
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
#'   First order and total Sobol' sensitivity indices are plotted for one state 
#'   variable (chosen by argument \code{state_plot}) and the parameters named
#'   in \code{pars_plot} against time. If no parameters are named in 
#'   \code{pars_plot}, the sensitivity indices for all parameters are plotted.
#'
#' @note 
#'   Not all arguments of \code{\link{plot.default}} can be passed by 
#'   \code{...}, for example \code{xlab} and \code{ylab} are fixed.
#'
#' @author Stefan Theers, Frank Weber
#' @seealso \code{\link{ODEsobol}, \link[sensitivity]{soboljansen},
#' \link[sensitivity]{sobolmartinez}}
#'
#' @examples
#' ##### Lotka-Volterra equations #####
#' LVmod <- function(Time, State, Pars) {
#'   with(as.list(c(State, Pars)), {
#'     Ingestion    <- rIng  * Prey * Predator
#'     GrowthPrey   <- rGrow * Prey * (1 - Prey/K)
#'     MortPredator <- rMort * Predator
#'     
#'     dPrey        <- GrowthPrey - Ingestion
#'     dPredator    <- Ingestion * assEff - MortPredator
#'     
#'     return(list(c(dPrey, dPredator)))
#'   })
#' }
#' LVpars  <- c("rIng", "rGrow", "rMort", "assEff", "K")
#' LVbinf <- c(0.05, 0.05, 0.05, 0.05, 1)
#' LVbsup <- c(1.00, 3.00, 0.95, 0.95, 20)
#' LVinit  <- c(Prey = 1, Predator = 2)
#' LVtimes <- c(0.01, seq(1, 50, by = 1))
#' set.seed(59281)
#' # Warning: The following code might take very long!
#' \donttest{
#' LVres_sobol <- ODEsobol(mod = LVmod,
#'                         pars = LVpars,
#'                         state_init = LVinit,
#'                         times = LVtimes,
#'                         n = 500,
#'                         rfuncs = "runif",
#'                         rargs = paste0("min = ", LVbinf,
#'                                        ", max = ", LVbsup),
#'                         sobol_method = "Martinez",
#'                         ode_method = "lsoda",
#'                         parallel_eval = TRUE,
#'                         parallel_eval_ncores = 2)
#' my_cols <- c("firebrick", "orange2", "dodgerblue", 
#'              "forestgreen", "black")
#' plot(LVres_sobol, colors_pars = my_cols)
#' plot(LVres_sobol, pars_plot = c("rGrow", "rMort"), state_plot = "Predator", 
#'      colors_pars = my_cols[2:3])
#' }
#' 
#' ##### A network of 4 mechanical oscillators connected in a circle #####
#' M_mat <- rep(2, 4)
#' K_mat <- diag(rep(2 * (2*pi*0.17)^2, 4))
#' K_mat[1, 2] <- K_mat[2, 3] <- 
#'   K_mat[3, 4] <- K_mat[1, 4] <- 2 * (2*pi*0.17)^2 / 10
#' D_mat <- diag(rep(0.05, 4))
#' library("ODEnetwork")
#' lfonet <- ODEnetwork(masses = M_mat, dampers = D_mat, springs = K_mat)
#' LFOpars <- c("k.1", "k.2", "k.3", "k.4",
#'              "d.1", "d.2", "d.3", "d.4")
#' LFObinf <- c(rep(0.2, 4), rep(0.01, 4))
#' LFObsup <- c(rep(20, 4), rep(0.1, 4))
#' lfonet <- setState(lfonet, state1 = rep(2, 4), state2 = rep(0, 4))
#' LFOtimes <- seq(25, 150, by = 2.5)
#' set.seed(1739)
#' # Warning: The following code might take very long!
#' \donttest{
#' suppressWarnings(
#'   LFOres_sobol <- ODEsobol(mod = lfonet,
#'                            pars = LFOpars,
#'                            times = LFOtimes,
#'                            n = 500,
#'                            rfuncs = "runif",
#'                            rargs = paste0("min = ", LFObinf,
#'                                           ", max = ", LFObsup),
#'                            sobol_method = "Martinez",
#'                            parallel_eval = TRUE,
#'                            parallel_eval_ncores = 2)
#' )
#' plot(LFOres_sobol, pars_plot = paste0("k.", 1:4), state_plot = "x.2",
#'      colors_pars = my_cols)
#' }
#' 
#' @export
#'

plot.ODEsobol <- function(x, pars_plot = NULL, state_plot = names(x)[1],
                          colors_pars = NULL, main_title = NULL,
                          legendPos = "outside", type = "l", ...) {
  
  ##### Input checks ###################################################
  
  assertClass(x, "ODEsobol")
  stopifnot(is.null(pars_plot) || is.character(pars_plot))
  assertCharacter(state_plot, len = 1)
  stopifnot(state_plot %in% names(x))
  stopifnot(is.null(colors_pars) || is.character(colors_pars))
  if(is.character(colors_pars) && is.null(pars_plot) &&
     length(colors_pars) < nrow((x[[which(names(x) == state_plot)]])$S) - 1){
    stop("\"colors_pars\" has to be at least of length",
         (nrow(x[[which(names(x) == state_plot)]]) - 1) / 3)
  } else if(is.character(colors_pars) && !is.null(pars_plot) &&
            length(colors_pars) < length(pars_plot)){
    stop("\"colors_pars\" has to be at least of length", length(pars_plot))
  }
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
  
  # Extract the matrices S and T:
  S <- (x[[state_idx]])$S
  T <- (x[[state_idx]])$T
  
  # Extract the timepoints:
  t.vec <- S[1, ]
  
  # Extract the parameter names, set "pars_plot" if not specified and perform
  # a validity check if parameter names are user-specified:
  all_pars <- rownames(S)[-1]
  if(is.null(pars_plot)){
    pars_plot <- all_pars
  } else{
    stopifnot(all(pars_plot %in% all_pars))
    if(any(duplicated(pars_plot))){
      pars_plot <- unique(pars_plot)
    }
  }
  k <- length(pars_plot)
  
  # Extreme values of the SA indices:
  minMaxS <- c(0.95 * min(S[pars_plot, ]), 1.05 * max(S[pars_plot, ]))
  minMaxT <- c(0.95 * min(T[pars_plot, ]), 1.05 * max(T[pars_plot, ]))
  
  # Set colors if not set by the user:
  if(is.null(colors_pars)){
    colors_pars <- grDevices::rainbow(k)
  }
  
  # Set main title if not set by the user:
  method_title <- attr(x, "sobol_method")
  if(is.null(main_title)){
    main_title <- paste0("Sobol' sensitivity indices for \"", state_plot, 
                         "\" (Sobol'-", method_title, " method)")
  }
  
  ##### Plot ###########################################################
  
  if(legendPos == "outside"){
    if(k > 6){
      # Multirow legend:
      legend_ncol <- ceiling(k / 2)
      oldpar <- graphics::par(mfrow = c(1, 2), oma = c(2.5, 0, 2, 0), mar = c(4, 4, 1, 2) + 0.2)
    } else{
      # One-row legend:
      legend_ncol <- k
      oldpar <- graphics::par(mfrow = c(1, 2),
                    oma = c(1.7, 0, 2, 0), mar = c(4, 4, 1, 2) + 0.2)
    }
  } else{
    oldpar <- graphics::par(mfrow = c(1, 2), oma = c(0, 0, 2, 0), mar = c(4, 4, 1, 2) + 0.2)
  }
  
  # First order SA indices:
  # First parameter:
  graphics::plot(x = t.vec, y = S[pars_plot[1], ],
       xlab = "Time", ylab = "First order Sobol' indices",
       type = type, col = colors_pars[1], ylim = minMaxS, ...)
  # All remaining parameters:
  if(k >= 2){
    for(i in 2:k) {
      graphics::lines(x = t.vec, y = S[pars_plot[i], ], type = type, 
            col = colors_pars[i], ...)
    }
  }
  # Legend:
  if(legendPos != "outside"){
    if(type %in% c("b", "o")){
      graphics::legend(legendPos, legend = pars_plot, col = colors_pars, bg = "white",
             lty = 1, pch = 1)
    } else if(type %in% c("l", "c", "s", "h")){
      graphics::legend(legendPos, legend = pars_plot, col = colors_pars, bg = "white", 
             lty = 1)
    } else if(type == "p"){
      graphics::legend(legendPos, legend = pars_plot, col = colors_pars, bg = "white", 
             pch = 1)
    }
  }
  
  # Total SA indices:
  # First parameter:
  graphics::plot(x = t.vec, y = T[pars_plot[1], ],
       xlab = "Time", ylab = "Total Sobol' indices",
       type = type, col = colors_pars[1], ylim = minMaxT, ...)
  # All remaining parameters:
  if(k >= 2){
    for(i in 2:k) {
      graphics::lines(x = t.vec, y = T[pars_plot[i], ], type = type, 
            col = colors_pars[i], ...)
    }
  }
  # Legend:
  if(legendPos != "outside"){
    if(type %in% c("b", "o")){
      graphics::legend(legendPos, legend = pars_plot, col = colors_pars, bg = "white",
             lty = 1, pch = 1)
    } else if(type %in% c("l", "c", "s", "h")){
      graphics::legend(legendPos, legend = pars_plot, col = colors_pars, bg = "white", 
             lty = 1)
    } else if(type == "p"){
      graphics::legend(legendPos, legend = pars_plot, col = colors_pars, bg = "white", 
             pch = 1)
    }
  }
  
  # Create the big common title for the two plots:
  graphics::mtext(main_title, side = 3, line = 0, outer = TRUE, cex = 1.2, font = 2)
  
  # Legend outside of the plotting region:
  if(legendPos == "outside"){
    # Dummy plot:
    oldpar2 <- graphics::par(mfrow = c(1, 1), oma = rep(0, 4), mar = rep(0, 4), new = TRUE)
    graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    # Draw the legend depending on the plot type:
    if(type %in% c("b", "o")){
      graphics::legend("bottom", legend = pars_plot, col = colors_pars, lty = 1, pch = 1, 
             bty = "n", xpd = TRUE, ncol = legend_ncol, inset = c(0, 0))
    } else if(type %in% c("l", "c", "s", "h")){
      graphics::legend("bottom", legend = pars_plot, col = colors_pars, lty = 1, 
             bty = "n", xpd = TRUE, ncol = legend_ncol, inset = c(0, 0))
    } else if(type == "p"){
      graphics::legend("bottom", legend = pars_plot, col = colors_pars, pch = 1, 
             bty = "n", xpd = TRUE, ncol = legend_ncol, inset = c(0, 0))
    }
    graphics::par(oldpar2)
  }
  
  graphics::par(oldpar)
  # For testing purposes:
  return(invisible(TRUE))
}
