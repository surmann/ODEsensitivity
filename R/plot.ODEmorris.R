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
#'   kind of the plot, choose between \code{"sep"} and \code{"trajec"} (see 
#'   details).
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
#'   Morris sensitivity indices are plotted for one state variable (chosen by 
#'   argument \code{state_plot}) and the parameters named in \code{pars_plot}. 
#'   If no parameters are named in \code{pars_plot}, the sensitivity indices for
#'   all parameters are plotted. There are two kinds of plots:
#'   \itemize{
#'     \item{\code{kind = "sep"}: }{separate plots of the Morris sensitivity 
#'       indices \eqn{\mu^*}{\mu*} and \eqn{\sigma} against time}
#'     \item{\code{kind = "trajec"}: }{plot of \eqn{\mu^*}{\mu*} against 
#'       \eqn{\sigma}}
#'   }
#'
#' @note 
#'   Not all plotting arguments can be passed by \code{...}, for example
#'   \code{xlab} and \code{ylab} are fixed.
#'
#' @author Stefan Theers, Frank Weber
#' @seealso \code{\link{ODEmorris}, \link[sensitivity]{morris}}
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
#' set.seed(7292)
#' LVres_morris <- ODEmorris(mod = LVmod,
#'                           pars = LVpars,
#'                           state_init = LVinit,
#'                           times = LVtimes,
#'                           binf = LVbinf,
#'                           bsup = LVbsup,
#'                           r = 500,
#'                           design = list(type = "oat", 
#'                                         levels = 10, grid.jump = 1),
#'                           scale = TRUE,
#'                           ode_method = "lsoda",
#'                           parallel_eval = TRUE,
#'                           parallel_eval_ncores = 2)
#' my_cols <- c("firebrick", "orange2", "dodgerblue", 
#'              "forestgreen", "black")
#' plot(LVres_morris, kind = "sep", colors_pars = my_cols)
#' plot(LVres_morris, pars_plot = c("rGrow", "rMort"), state_plot = "Predator", 
#'      kind = "trajec", colors_pars = my_cols[2:3])
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
#' LFOres_morris <- ODEmorris(mod = lfonet,
#'                            pars = LFOpars,
#'                            times = LFOtimes,
#'                            binf = LFObinf,
#'                            bsup = LFObsup,
#'                            r = 500,
#'                            design = list(type = "oat", 
#'                                          levels = 10, grid.jump = 1),
#'                            scale = TRUE,
#'                            parallel_eval = TRUE,
#'                            parallel_eval_ncores = 2)
#' plot(LFOres_morris, pars_plot = paste0("k.", 1:4), state_plot = "x.2",
#'      kind = "sep", colors_pars = my_cols)
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
  k_all <- (nrow(x[[state_idx]]) - 1) / 3
  pars_tmp <- rownames(x[[state_idx]])[2:(k_all + 1)]
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
