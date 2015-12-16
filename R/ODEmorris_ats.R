#' @title Morris SA for ODEs at All Timepoints Simultaneously
#'
#' @description
#' \code{ODEmorris_ats} performs a sensitivity analysis for
#' ordinary differential equations using Morris's elementary effects 
#' screening method. The analysis is done for one output variable at all
#' timepoints simultaneously using \code{\link[sensitivity]{morris_matrix}}.
#'
#' @param mod [\code{function(Time, State, Pars)}]\cr
#'   model to examine, cf. example below.
#' @param pars [\code{character(k)}]\cr
#'   vector of \code{k} input variable names.
#' @param yini [\code{numeric(z)}]\cr
#'   vector of \code{z} initial values.
#' @param times [\code{numeric}]\cr
#'   points of time at which the SA should be executed
#'   (vector of arbitrary length). Also the
#'   first point of time must be positive.
#' @param y_idx [\code{integer(1)}]\cr
#'   index of the output variable to be analyzed. Defaults to 1, so the first
#'   output variable is used.
#' @param seed [\code{numeric(1)}]\cr
#'   seed.
#' @param binf [\code{numeric(k)}]\cr
#'   vector of lower borders of possible input parameter values.
#'   If they are all equal, a single value can be set.
#' @param bsup [\code{numeric(k)}]\cr
#'   vector of upper borders of possible input parameter values.
#'   If they are all equal, a single value can be set.
#' @param r [\code{integer(1)}]\cr
#'   number of repetitions of the \code{design},
#'   cf. \code{\link[sensitivity]{morris}}.
#' @param design [\code{list}]\cr
#'   a list specifying the design type and its parameters,
#'   cf. \code{\link[sensitivity]{morris}}.
#' @param scale [\code{logical(1)}]\cr
#'   if \code{TRUE}, scaling is done for the input design of experiments after 
#'   building the design and before calculating the elementary effects,
#'   cf. \code{\link[sensitivity]{morris}}. Defaults to \code{FALSE}, but it is
#'   highly recommended to use \code{scale = TRUE} if the factors have different
#'   orders of magnitude.
#'
#' @return list of class \code{morrisRes_ats} with
#'   \itemize{
#'     \item \code{res}, a matrix containing the Morris SA results (i.e.
#'       \code{mu, mu.star} and \code{sigma} for every parameter) for one point
#'       of time per column
#'     \item \code{pars}, the parameter names.
#'   }
#' 
#' @details \code{ODEmorris_ats} is faster than \code{ODEmorris} since 
#' \code{\link[sensitivity]{morris_matrix}} can handle matrix output for its 
#' model function. Thus, one random matrix of parameter combinations can be 
#' generated and the adopted model function for 
#' \code{\link[sensitivity]{morris_matrix}} can 
#' return the values of the output variable for all parameter combinations 
#' (rows) and all timepoints (columns) together. This saves time since the
#' adopted model function for 
#' \code{\link[sensitivity]{morris_matrix}} mainly relies on executing 
#' \code{\link[deSolve]{ode}} from the
#' package \code{deSolve}. Using \code{\link[deSolve]{ode}} for all timepoints
#' simultaneously is a lot faster than \code{apply}-ing it on every timepoint
#' separately.
#'
#' @author Frank Weber
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
#' FHNpars  <- c(a = 0.2,     # parameter a
#'               b = 0.3,     # parameter b
#'               s = 3)       # parameter s (= c in the original notation)
#'
#' FHNyini  <- c(Voltage = -1, Current = 1)
#' FHNtimes <- seq(0.1, 100, by = 10)
#'
#' FHNres_ats <- ODEmorris_ats(mod = FHNmod,
#'                             pars = names(FHNpars),
#'                             yini = FHNyini,
#'                             times = FHNtimes,
#'                             y_idx = 1,        # voltage only
#'                             seed = 2015,
#'                             binf = c(0.18, 0.18, 2.8),
#'                             bsup = c(0.22, 0.22, 3.2),
#'                             r = 25,
#'                             design =
#'                               list(type = "oat", levels = 100, 
#'                                    grid.jump = 1),
#'                             scale = TRUE)
#'
#' @seealso \code{\link[sensitivity]{morris}},
#'   \code{\link[sensitivity]{morris_matrix}},
#'   \code{\link{ODEmorris_aos}},
#'   \code{\link{plot.morrisRes_ats}}
#'
#' @note \code{ODEmorris_ats} is only purposed for analyzing one output 
#' variable. If multiple output variables shall be analysed simultaneously,
#' please use \code{\link{ODEmorris_aos}}.
#' 
#' \code{\link[deSolve]{ode}} or rather its standard solver \code{lsoda}
#'   sometimes cannot solve an ODE system if unrealistic parameters
#'   are sampled by \code{\link[sensitivity]{morris}}. Hence
#'   \code{NA}s might occur in the Morris sensitivity results, such
#'   that \code{\link{ODEmorris_ats}} fails for one or many points of time!
#'   For this reason, if \code{NA}s occur, please make use of
#'   \code{\link{ODEsobol}} instead or
#'   restrict the input parameter value intervals usefully using
#'   \code{binf} and \code{bsup}!
#'   
#'   If \code{\link[sensitivity]{morris_matrix}} throws a warning message saying
#'   "In ... keeping ... repetitions out of ...", try using a bigger number of 
#'   \code{levels} in the \code{design} argument.
#'
#' @references J. O. Ramsay, G. Hooker, D. Campbell and J. Cao, 2007,
#'   \emph{Parameter estimation for differential equations: a generalized 
#'   smoothing approach}, Journal of the Royal Statistical Society, Series B, 
#'   69, Part 5, 741--796.
#'
#' @export
#' @import
#'   checkmate
#' @importFrom deSolve ode
#' @importFrom sensitivity morris_matrix
#'

ODEmorris_ats <- function(mod,
                          pars,
                          yini,
                          times,
                          y_idx = 1,
                          seed = 2015,
                          binf = 0,
                          bsup = 1,
                          r = 25,
                          design =
                            list(type = "oat", levels = 100, grid.jump = 1),
                          scale = FALSE) {
  
  ##### Plausibilitaet #################################################
  ## stopifnot(!missing(...))
  assertFunction(mod)
  assertCharacter(pars)
  assertNumeric(yini)
  assertNumeric(times, lower = 0, finite = TRUE, unique = TRUE)
  times <- sort(times)
  stopifnot(!any(times == 0))
  assertIntegerish(y_idx)
  assertNumeric(seed)
  assertNumeric(binf)
  notOk <- length(binf) != length(pars) & length(binf) != 1
  if(notOk)
    stop("binf must be of length 1 or of the same length as pars!")
  assertNumeric(bsup)
  notOk <- length(bsup) != length(pars) & length(bsup) != 1
  if(notOk)
    stop("bsup must be of length 1 or of the same length as pars!")
  assertIntegerish(r, len = 1)
  if(r < 1)
    stop("r must be greater or equal to 1.")
  assertList(design)
  
  ##### Vorarbeiten ####################################################
  set.seed(seed)
  # Anzahl Parameter:
  k <- length(pars)
  # Anzahl Outputgroessen:
  z <- length(yini)
  # Anzahl Zeitpunkte von Interesse:
  timesNum <- length(times)
  # Forme DGL-Modell um, sodass fuer morris_matrix()-Argument "model_matrix" 
  # passend:
  model_fit <- function(X){
    # X   - (nxk)-Matrix mit den n einzugebenden Parameter-Konstellationen
    #       als Zeilen
    colnames(X) <- pars
    res <- apply(X, 1, function(x){
      ode(yini, times = c(0, times), 
          mod, parms = x)[2:(timesNum + 1), y_idx + 1]
    })
    # (Jede Zeile in "res" steht fuer einen Zeitpunkt.)
    
    if(timesNum == 1){
      # Korrektur noetig, falls timesNum == 1:
      res <- matrix(res)
    } else{
      # Transponiere die Ergebnis-Matrix, sodass jede Spalte fuer einen
      # Zeitpunkt steht:
      res <- t(res)
    }
    return(res)
  }
  
  ##### Sensitivitaet ##################################################
  x <- morris_matrix(model_matrix = model_fit, factors = k, r = r, 
                     design = design, binf = binf, bsup = bsup, scale = scale)
  mu <- lapply(x$ee_by_col, colMeans)
  mu.star <- lapply(x$ee_by_col, abs)
  mu.star <- lapply(mu.star, colMeans)
  sigma <- lapply(x$ee_by_col, function(M){
    apply(M, 2, sd)
  })
  out_y_idx <- mapply(c, mu, mu.star, sigma, SIMPLIFY = TRUE)
  out_y_idx <- rbind(times, out_y_idx)
  rownames(out_y_idx) <- c("time", paste0("mu_", pars), 
                           paste0("mu.star_", pars),
                           paste0("sigma_", pars))
  
  # Warnungen, falls NAs auftreten (unrealistische Parameter => nicht
  # loesbare ODEs):
  if(any(is.na(out_y_idx[1:(1 + k*2), ]))){
    warning("deSolve/ lsoda cannot solve the ODE system!
            This might be due to arising unrealistic parameters by means of 
            Morris Screening. Use ODEsobol() instead or set binf and bsup 
            differently!")
  } else if(all(is.na(out_y_idx[(2 + k*2):(1 + k*3), ])) && r == 1){
    warning("Calculation of sigma requires r >= 2.")
  } else if(any(is.na(out_y_idx[(2 + k*2):(1 + k*3), ]))){
    warning("NAs for sigma. This might be due to r being too small.")
  }
  
  # Rueckgabe:
  res <- list(res = out_y_idx, pars = pars)
  class(res) <- "morrisRes_ats"
  return(res)
}
