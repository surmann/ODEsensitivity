#' @title
#' Plotting the results of Morris SA
#'
#' @description
#' \code{plot} plots the results of Morris SA.
#'
#' @details
#' \code{plot} with \code{type = "sep"} plots mu.star and
#'   sigma separately versus time.
#'
#' \code{plot} with \code{type = "trajec"} plots mu.star versus
#'   sigma for every point of time.
#'
#' @param x [\code{morrisRes_trafo}]\cr
#'   resulting output of \code{\link{ODEmorris_trafo}}, of class 
#'   \code{morrisRes_trafo}.
#' @param type [\code{character(1)}]\cr
#'   plot type, choose between \code{"sep"} and \code{"trajec"}.
#' @param main_title [\code{character(1)}]\cr
#'   title for the plot. If \code{type = "sep"}, this is the overall title for
#'   the two separate plots. Defaults to NULL, so a standard title is generated.
#' @param legendPos [\code{character(1)}]\cr
#'   legend position, default is \code{"topleft"}.
#' @param ... additional arguments passed to \code{\link{plot}}.
#'
#' @return TRUE (invisible; for testing purposes).
#'
#' @method plot morrisRes_trafo
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
#' FHNres_trafo <- ODEmorris_trafo(mod = FHNmod,
#'                                 pars = c("a", "b", "s"),
#'                                 yini = FHNyini,
#'                                 times = FHNtimes,
#'                                 ode_method = "adams",
#'                                 seed = 2015,
#'                                 binf = c(0.18, 0.18, 2.8),
#'                                 bsup = c(0.22, 0.22, 3.2),
#'                                 r = 25,
#'                                 design =
#'                                   list(type = "oat", levels = 100, 
#'                                        grid.jump = 1),
#'                                 scale = TRUE,
#'                                 trafo = function(Y) Y[, 1],    # voltage only
#'                                 ncores = 2)  
#'
#' # Plots:
#' plot(FHNres, type = "sep")
#' plot(FHNres, type = "trajec")
#'
#' @note Unfortunately, the passing of arguments (e.g. "main") does not work
#'   correctly.
#'
#' @author Stefan Theers
#' @seealso \code{\link{ODEmorris_trafo}},
#'   \code{\link[sensitivity]{morris}}
#'
#' @export
#' @import
#'   checkmate
#'

plot.morrisRes_trafo <- function(x, type = "sep", main_title = NULL, 
                                 legendPos = "topleft", ...) {

  ##### Check input #################################################
  assertClass(x, "morrisRes_trafo")
  assertCharacter(type, len = 1)
  notOk <- !any(rep(type, 2) == c("sep", "trajec"))
  if(notOk)
    stop("type must be one of \"sep\" or \"trajec\"!")
  assertCharacter(legendPos, len = 1)
  notOk <- !any(rep(legendPos, 9) == c("bottomright", "bottom",
    "bottomleft", "left", "topleft", "top", "topright", "right",
    "center"))
  if(notOk)
    stop("legendPos must be one of \"bottomright\", \"bottom\",
      \"bottomleft\", \"left\", \"topleft\", \"top\", \"topright\",
      \"right\", \"center\"!")
  stopifnot(is.character(main_title) && length(main_title) == 1 ||
              is.null(main_title))

  ##### Plot ###########################################################
  
  # Extrahiere die Parameter-Namen:
  k <- (nrow(x) - 1) / 3
  pars_tmp <- rownames(x)[2:(k + 1)]
  pars <- substr(pars_tmp, start = 4, stop = nchar(pars_tmp))
  
  # Separate Plots fuer mu.star und sigma:
  if(type == "sep"){
    oldpar <- par(mfrow = c(1, 2), mar = c(4, 4, 1, 2) + 0.2,
                  oma = c(0, 0, 2, 0))
    # Erstelle die separaten Plots:
    plotSep(x, pars, legendPos, ...)
    # Erstelle die Gesamtueberschrift:
    if(is.null(main_title)){
      main_title <- "Morris SA"
    }
    mtext(main_title, side = 3, line = 0, outer = TRUE, cex = 1.2, font = 2)
    par(oldpar)
  }
  
  # Trajectories:
  if(type == "trajec"){
    # Erstelle die Ueberschrift:
    if(is.null(main_title)){
      main_title <- "Morris SA: Trajectories"
    }
    plotTrajectories(x, pars, legendPos, main_title, ...)
  }
  
  # For testing purposes:
  return(invisible(TRUE))
}
