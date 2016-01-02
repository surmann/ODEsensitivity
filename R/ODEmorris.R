#' @title Morris SA for ODEs
#'
#' @description
#' \code{ODEmorris} performs a sensitivity analysis for
#' ordinary differential equations using
#' Morris's elementary effects screening method.
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
#' @param trafo [\code{function}]\cr
#'   function to transform \code{z > 1} output variables to
#'   IR [only needed, if \code{z > 1}]. Must be able to deal with a
#'   matrix.
#' @param ncores [\code{integer(1)}]\cr
#'   number of processor cores to be used for calculating the sensitivity
#'   indices. Must be between 1 and 4.
#'
#' @return List of class \code{morrisRes} with
#'   \itemize{
#'     \item \code{res}, the matrix of Morris SA results
#'       (i.e. \code{mu, mu.star} and \code{sigma}) for every point of
#'       time in the \code{times} vector and
#'     \item \code{pars}, the parameter names.
#'   }
#'
#' @note Package \code{parallel} is needed for this function. However, this 
#'   function is outdated and the use of \code{\link{ODEmorris_ats}} (for one 
#'   output variable) or \code{\link{ODEmorris_aos}} (for multiple output 
#'   variables) is recommended if no transformation has to be made to the 
#'   \code{\link[deSolve]{ode}}-results.
#' 
#' \code{\link[deSolve]{ode}} sometimes cannot solve an ODE system if 
#'   unrealistic parameter combinations are sampled by 
#'   \code{\link[sensitivity]{morris}}. Hence \code{NA}s might occur in the 
#'   Morris sensitivity results, such that \code{\link{ODEmorris}} fails for 
#'   one or many points of time! For this reason, if \code{NA}s occur, please 
#'   make use of \code{\link{ODEsobol}} instead or
#'   restrict the input parameter value intervals usefully using
#'   \code{binf}, \code{bsup} and \code{scale = TRUE}. It is also helpful to try
#'   another ODE-solver (argument \code{ode_method}). Problems are known for the
#'   \code{ode_method}s \code{"euler"}, \code{"rk4"} and \code{"ode45"}. 
#'   In contrast, the \code{ode_method}s \code{"vode"}, \code{"bdf"}, 
#'   \code{"bdf_d"}, \code{"adams"}, \code{"impAdams"} and \code{"impAdams_d"} 
#'   might be even faster than the standard \code{ode_method} \code{"lsoda"}.
#'   
#'   If \code{\link[sensitivity]{morris}} throws a warning message saying
#'   "In ... keeping ... repetitions out of ...", try using a bigger number of 
#'   \code{levels} in the \code{design} argument.
#'
#' @author Stefan Theers
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
#' FHNres <- ODEmorris(mod = FHNmod,
#'                     pars = c("a", "b", "s"),
#'                     yini = FHNyini,
#'                     times = FHNtimes,
#'                     ode_method = "adams",
#'                     seed = 2015,
#'                     binf = c(0.18, 0.18, 2.8),
#'                     bsup = c(0.22, 0.22, 3.2),
#'                     r = 25,
#'                     design =
#'                         list(type = "oat", levels = 100, grid.jump = 1),
#'                     scale = TRUE,
#'                     trafo = function(Y) Y[, 1],    # voltage only
#'                     ncores = 4)
#'
#' @seealso \code{\link[sensitivity]{morris}, \link{ODEmorris_ats},
#' \link{ODEmorris_aos}, \link{plot.morrisRes}}
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
#' @importFrom sensitivity morris
#'

ODEmorris <- function(mod,
                      pars,
                      yini,
                      times,
                      ode_method = "lsoda",
                      seed = 2015,
                      binf = 0,
                      bsup = 1,
                      r = 25,
                      design =
                        list(type = "oat", levels = 100, grid.jump = 1),
                      scale = TRUE,
                      trafo = function(Y) rowSums(Y^2),
                      ncores = 1) {

  ##### Plausibilitaet #################################################
  ## stopifnot(!missing(...))
  assertFunction(mod)
  assertCharacter(pars)
  assertNumeric(yini)
  assertNumeric(times, lower = 0, finite = TRUE, unique = TRUE)
  times <- sort(times)
  stopifnot(!any(times == 0))
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
  assertIntegerish(r)
  assertList(design)
  assertFunction(trafo)
  notOk <- !testVector(trafo(matrix(1:30, nrow = 6)), len = 6)
  if(notOk)
    stop("Make sure that trafo() transforms matrices to suitable vectors!")
  assertIntegerish(ncores, lower = 1L, upper = 4L)
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop(paste("Package \"parallel\" needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }

  ##### Vorarbeiten ####################################################
  # set.seed(seed)
  # Anzahl Parameter:
  k <- length(pars)
  # Anzahl Outputgroessen:
  z <- length(yini)
  # Anzahl Zeitpunkte von Interesse:
  timesNum <- length(times)

  # Umformen DGL-Modell, sodass fuer morris()-Argument model passend
  modFun <- function(X, pot) {
    # X   - (nxk)-Matrix
    # pot - point of time
    colnames(X) <- pars
    res <-
      t(apply(X, 1, function(x)
              ode(yini, times = c(0, pot), mod, parms = x,
                  method = ode_method)[2, 2:(z+1)]))
    # Transformation der Output-Variablen nach IR:
    trafo(res)
  }

  ##### Sensitivitaet ##################################################
  
  # performing Morris SA for 1 point of time:
  oneRun <- function(pot) {
    # Setze hier erst den Seed, sodass immer die gleiche Zufallsmatrix
    # fuer die Parameter-Konstellationen verwendet wird:
    set.seed(seed)
    x <- morris(model = modFun, factors = k, r = r, pot = pot,
                design = design, binf = binf, bsup = bsup, scale = scale)
    # analog zur Hilfeseite von morris()/ weniger primitiv als Schleife:
    mu <- colMeans(x$ee)
    mu.star <- colMeans(abs(x$ee))
    sigma <- apply(x$ee, 2, sd)

    # Ergebnisse:
    res <- c(pot, mu, mu.star, sigma)
    names(res) <- c("time", paste0("mu_", pars), paste0("mu.star_", pars),
                    paste0("sigma_", pars))
    return(res)
  }

  cl <- parallel::makeCluster(rep("localhost", ncores), type = "SOCK")
  # clusterSetRNGStream(cl)
  parallel::clusterExport(cl, list("mod", "modFun", "times", "timesNum", "pars",
                                   "yini", "z", "r", "design", "oneRun",
                                   "morris", "ode", "trafo", "k"),
                          envir = environment())
  ## parallel::clusterExport(cl, list(ls(), "sobol2007", "ode"),
  ##                         envir = environment())
  res <- parallel::parSapply(cl, times, oneRun)
  parallel::stopCluster(cl)

  # Warnungen, falls NAs auftreten (unrealistische Paramter => nicht
  # loesbare ODEs):
  if(any(is.na(res[1:(1 + k*2), ]))){
    warning(paste("The ODE system can't be solved. This might be due to", 
      "arising unrealistic parameters by means of Morris Screening. Use",
      "ODEsobol() instead or set binf and bsup differently together with",
      "scale = TRUE. It might also be helpful to try another ODE-solver by",
      "using the \"ode_method\"-argument."))
  } else if(all(is.na(res[(2 + k*2):(1 + k*3), ])) && r == 1){
    warning("Calculation of sigma requires r >= 2.")
  } else if(any(is.na(res[(2 + k*2):(1 + k*3), ]))){
    warning("NAs for sigma. This might be due to r being too small.")
  }
  
  # Rueckgabe:
  class(res) <- "morrisRes"
  return(res)
}
