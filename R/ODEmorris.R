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
#' @param trafo [\code{function}]\cr
#'   function to transform \code{z > 1} output variables to
#'   IR [only needed, if \code{z > 1}]. Must be able to deal with a
#'   matrix.
#' @param ncores [\code{integer(1)}]\cr
#'   number of processor cores to be used for calculating the sensitivity
#'   indices. Must be between 1 and 4.
#'
#' @return list of class \code{morrisRes} with
#'   \itemize{
#'     \item \code{res}, the matrix of Morris SA results
#'       (i.e. \code{mu, mu.star} and \code{sigma}) for every point of
#'       time in the \code{times} vector and
#'     \item \code{pars}, the parameter names.
#'   }
#'
#'
#' @examples
#' ##### FitzHugh-Nagumo equations (Ramsay et al, 2007)
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
#' FHNres <- ODEmorris(mod = FHNmod,
#'                     pars = names(FHNpars),
#'                     yini = FHNyini,
#'                     times = FHNtimes,
#'                     seed = 2015,
#'                     binf = c(0.18, 0.18, 2.8),
#'                     bsup = c(0.22, 0.22, 3.2),
#'                     r = 25,
#'                     design =
#'                         list(type = "oat", levels = 100, grid.jump = 1),
#'                     trafo = function(Y) Y[, 1],    # voltage only
#'                     ncores = 4)
#'
#' @seealso \code{\link[sensitivity]{morris}},
#'   \code{\link{plot.morrisRes}}
#'
#' @note \code{\link[deSolve]{ode}} or rather its standard solver \code{lsoda}
#'   sometimes cannot solve an ODE system if unrealistic parameters
#'   are sampled by \code{\link[sensitivity]{morris}}. Hence
#'   \code{NA}s might occur in the Morris sensitivity results, such
#'   that \code{\link{ODEmorris}} fails for one or many points of time!
#'   For this reason, if \code{NA}s occur, please make use of
#'   \code{\link{ODEsobol}} instead or
#'   restrict the input parameter value intervals usefully using
#'   \code{binf} and \code{bsup}!
#'
#'
#' @export
#' @import
#'   checkmate
#'   deSolve
#'   sensitivity
#'   boot
#'   parallel
#'   BBmisc
#'

ODEmorris <- function(mod,
                      pars,
                      yini,
                      times,
                      seed = 2015,
                      binf = 0,
                      bsup = 1,
                      r = 25,
                      design =
                        list(type = "oat", levels = 100, grid.jump = 1),
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

  ##### Vorarbeiten ####################################################
  set.seed(seed)
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
              ode(yini, times = c(0, pot), mod, parms = x)[2, 2:(z+1)]))

    # Transformation der Output-Variablen nach IR:
    trafo(res)
  }

  ##### Sensitivitaet ##################################################
  
  # performing Morris SA for 1 point of time:
  oneRun <- function(pot) {
    set.seed(seed)
    x <- morris(model = modFun, factors = k, r = r, pot = pot,
                design = design, binf = binf, bsup = bsup)
    # analog zur Hilfeseite von morris()/ weniger primitiv als Schleife:
    mu <- colMeans(x$ee)
    mu.star <- colMeans(abs(x$ee))
    sigma <- apply(x$ee, 2, sd)

    # Ergebnisse:
    res <- c(pot, mu, mu.star, sigma)
    # k <- ncol(x$ee)
    names(res) <- c("time",
                    paste("mu", 1:k, sep = ""),
                    paste("mu.star", 1:k, sep = ""),
                    paste("sigma", 1:k, sep = ""))
    return(res)
  }

  cl <- makeCluster(rep("localhost", ncores), type = "SOCK")
  # clusterSetRNGStream(cl)
  clusterExport(cl, list("mod", "modFun", "times", "timesNum", "pars",
                         "yini", "z", "r", "design", "oneRun",
                         "morris", "ode", "trafo", "k"),
                envir = environment())
  ## clusterExport(cl, list(ls(), "sobol2007", "ode"), envir = environment())
  res <- parSapply(cl, times, oneRun)
  stopCluster(cl)

  # Warnungen, falls NAs auftreten (unrealistische Paramter => nicht
  # loesbare ODEs):
  if(any(is.na(res)))
    warning("deSolve/ lsoda cannot solve the ODE system!
This might be due to arising unrealistic parameters by means of Morris
Screening. Use ODEsobol() instead or set binf and bsup differently!")

  # Rueckgabe:
  res <- list(res = res, pars = pars)
  res <- setClasses(res, "morrisRes")
  return(res)
}
