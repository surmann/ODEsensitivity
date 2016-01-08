#' @title Morris SA for ODEs at All Timepoints Simultaneously
#'
#' @description
#' \code{ODEmorris_ats} performs a sensitivity analysis for
#' ordinary differential equations using Morris's elementary effects 
#' screening method. The analysis is done for one output variable at all
#' timepoints simultaneously using \code{\link[sensitivity]{morris_matrix}} from
#' the package \code{sensitivity}.
#'
#' @param mod [\code{function(Time, State, Pars)}]\cr
#'   model to examine, cf. example below.
#' @param pars [\code{character(k)}]\cr
#'   vector of \code{k} input variable names.
#' @param yini [\code{numeric(z)}]\cr
#'   vector of \code{z} initial values. Must be named and must not contain
#'   duplicated names.
#' @param times [\code{numeric}]\cr
#'   points of time at which the SA should be executed
#'   (vector of arbitrary length). Also the
#'   first point of time must be positive.
#' @param y_analyzed [\code{character(1)}]\cr
#'   name of the \code{yini}-variable to be analyzed. Defaults to the name of 
#'   the first \code{yini}-variable.
#' @param ode_method [\code{character(1)}]\cr
#'   method to be used for solving the differential equations, see 
#'   \code{\link[deSolve]{ode}}. Defaults to \code{"lsoda"}.
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
#'   cf. \code{\link[sensitivity]{morris}}. Defaults to \code{TRUE}, which is
#'   highly recommended if the factors have different orders of magnitude, see
#'   \code{\link[sensitivity]{morris}}.
#'
#' @return List of class \code{morrisRes_ats} with
#'   \itemize{
#'     \item \code{res}, a list of length \code{length(yini)} containing in each
#'       element a matrix for one yini-variable. The matrices itself contain in 
#'       their rows the Morris SA results (i.e. \code{mu, mu.star} and 
#'       \code{sigma} for every parameter) for all timepoints (columns).
#'     \item \code{y_analyzed}, the name of the analyzed \code{yini}-variable.
#'   }
#' 
#' @details \code{ODEmorris_ats} is faster than \code{ODEmorris} since 
#' \code{\link[sensitivity]{morris_matrix}} can handle matrix output for 
#' its model function. In \code{ODEmorris_ats}, a model function is
#' created returning a matrix with the \code{\link[deSolve]{ode}}-results for 
#' all timepoints (one per column). Thus, \code{\link[deSolve]{ode}} only needs
#' to be executed once.
#' 
#' @note \code{ODEmorris_ats} is only purposed for analyzing one output 
#' variable. If multiple output variables shall be analyzed simultaneously,
#' please use \code{\link{ODEmorris_aos}}.
#' 
#' \code{\link[deSolve]{ode}} sometimes cannot solve an ODE system if 
#'   unrealistic parameter combinations are sampled by 
#'   \code{\link[sensitivity]{morris_matrix}}. Hence \code{NA}s might occur in 
#'   the Morris sensitivity results, such that \code{\link{ODEmorris_ats}} fails
#'   for one or many points of time! For this reason, if \code{NA}s occur, 
#'   please make use of \code{\link{ODEsobol_ats}} instead or
#'   restrict the input parameter value intervals usefully using
#'   \code{binf}, \code{bsup} and \code{scale = TRUE}. It is also helpful to try
#'   another ODE-solver (argument \code{ode_method}). Problems are known for the
#'   \code{ode_method}s \code{"euler"}, \code{"rk4"} and \code{"ode45"}. 
#'   In contrast, the \code{ode_method}s \code{"vode"}, \code{"bdf"}, 
#'   \code{"bdf_d"}, \code{"adams"}, \code{"impAdams"} and \code{"impAdams_d"} 
#'   might be even faster than the standard \code{ode_method} \code{"lsoda"}.
#'   
#'   If \code{\link[sensitivity]{morris_matrix}} throws a warning message saying
#'   "In ... keeping ... repetitions out of ...", try using a bigger number of 
#'   \code{levels} in the \code{design} argument.
#' 
#' @author Frank Weber
#' @references J. O. Ramsay, G. Hooker, D. Campbell and J. Cao, 2007,
#'   \emph{Parameter estimation for differential equations: a generalized 
#'   smoothing approach}, Journal of the Royal Statistical Society, Series B, 
#'   69, Part 5, 741--796.
#' @seealso \code{\link[sensitivity]{morris}},
#'   \code{\link[sensitivity]{morris_matrix}},
#'   \code{\link{ODEmorris_aos}},
#'   \code{\link{plot.morrisRes_ats}}
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
#' @import
#'   checkmate
#' @importFrom deSolve ode
#' @importFrom sensitivity morris_matrix
#' @export
#'

ODEmorris_ats <- function(mod,
                          pars,
                          yini,
                          times,
                          y_analyzed = names(yini)[1],
                          ode_method = "lsoda",
                          seed = 2015,
                          binf = 0,
                          bsup = 1,
                          r = 25,
                          design =
                            list(type = "oat", levels = 100, grid.jump = 1),
                          scale = TRUE) {
  
  ##### Plausibilitaet #################################################
  assertFunction(mod)
  assertCharacter(pars)
  assertNumeric(yini)
  checkNamed(yini, type = "unique")
  assertNumeric(times, lower = 0, finite = TRUE, unique = TRUE)
  times <- sort(times)
  stopifnot(!any(times == 0))
  assertCharacter(y_analyzed, len = 1)
  stopifnot(y_analyzed %in% names(yini))
  stopifnot(ode_method %in% c("lsoda", "lsode", "lsodes","lsodar","vode", 
                              "daspk", "euler", "rk4", "ode23", "ode45", 
                              "radau", "bdf", "bdf_d", "adams", "impAdams", 
                              "impAdams_d" ,"iteration"))
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
  # Index der zu analysierenden y-Variable:
  y_idx <- which(y_analyzed == names(yini))
  # Forme DGL-Modell um, sodass fuer morris_matrix()-Argument "model_matrix" 
  # passend:
  model_fit <- function(X){
    # X   - (nxk)-Matrix mit den n einzugebenden Parameter-Konstellationen
    #       als Zeilen
    colnames(X) <- pars
    res <- apply(X, 1, function(x){
      ode(yini, times = c(0, times), mod, parms = x, 
          method = ode_method)[2:(timesNum + 1), y_idx + 1]
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
  x <- morris_matrix(model = model_fit, factors = k, r = r, 
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
    warning(paste("The ODE system can't be solved. This might be due to", 
      "arising unrealistic parameters by means of Morris Screening. Use",
      "ODEsobol() instead or set binf and bsup differently together with",
      "scale = TRUE. It might also be helpful to try another ODE-solver by",
      "using the \"ode_method\"-argument."))
  } else if(all(is.na(out_y_idx[(2 + k*2):(1 + k*3), ])) && r == 1){
    warning("Calculation of sigma requires r >= 2.")
  } else if(any(is.na(out_y_idx[(2 + k*2):(1 + k*3), ]))){
    warning("NAs for sigma. This might be due to r being too small.")
  }
  
  # Rueckgabe:
  res <- list(res = out_y_idx, y_analyzed = y_analyzed)
  class(res) <- "morrisRes_ats"
  return(res)
}
