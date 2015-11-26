#' @title Morris SA for ODEs at All Timepoints Simultaneously
#'
#' @description
#' \code{ODEmorris_ats} performs a sensitivity analysis for
#' ordinary differential equations at all timepoints 
#' simultaneously using Morris's elementary effects 
#' screening method.
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
#' @param ncores [\code{integer(1)}]\cr
#'   number of processor cores to be used for calculating the sensitivity
#'   indices. Must be between 1 and 4.
#'
#' @return list of class \code{morrisRes_ats} with
#'   \itemize{
#'     \item \code{res}, the list of Morris SA results for every output
#'       variable (each element of the list contains a matrix of
#'       \code{mu, mu.star} and \code{sigma} for every parameter and
#'       every point of time in the \code{times} vector)
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
#'   boot
#'   parallel
#'   BBmisc
#'

ODEmorris_ats <- function(mod,
                          pars,
                          yini,
                          times,
                          seed = 2015,
                          binf = 0,
                          bsup = 1,
                          r = 25,
                          design =
                            list(type = "oat", levels = 100, grid.jump = 1),
                          scale = FALSE,
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
  assertIntegerish(r, len = 1)
  assertList(design)
  assertIntegerish(ncores, lower = 1L, upper = 4L)
  
  ##### Vorarbeiten ####################################################
  set.seed(seed)
  # Anzahl Parameter:
  k <- length(pars)
  # Anzahl Outputgroessen:
  z <- length(yini)
  # Anzahl Zeitpunkte von Interesse:
  timesNum <- length(times)
  
  ##### Sensitivitaet ##################################################
  
  # Morris SA fuer eine Outputvariable:
  one_y <- function(y_idx){
    # Forme DGL-Modell um, sodass fuer morris_ats()-Argument "model_matrix" 
    # passend:
    model_fit <- function(X){
      # X   - (nxk)-Matrix mit den n einzugebenden Parameter-Konstellationen
      #       als Zeilen
      colnames(X) <- pars
      res <- t(apply(X, 1, function(x){
        ode(yini, times = c(0, times), FHNmod, parms = x)[2:(timesNum+1), y_idx+1]
      }))
      return(res)
    }
    
    x <- morris_ats(model_matrix = model_fit, pars = pars, p = k, r = r, 
                    design = design, binf = binf, bsup = bsup, scale = scale)
    mu <- lapply(x$ee, colMeans)
    mu.star <- lapply(x$ee, abs)
    mu.star <- lapply(mu.star, colMeans)
    # mu.star <- lapply(x$ee, function(M){
    #   apply(M, 2, function(x) mean(abs(x)))
    # })
    sigma <- lapply(x$ee, function(M){
      apply(M, 2, sd)
    })
    out_y_idx <- mapply(c, mu, mu.star, sigma, SIMPLIFY = TRUE)
    out_y_idx <- rbind(times, out_y_idx)
    rownames(out_y_idx) <- c("time", paste0("mu_", pars), 
                             paste0("mu.star_", pars),
                             paste0("sigma_", pars))
    
    # Warnungen, falls NAs auftreten (unrealistische Paramter => nicht
    # loesbare ODEs):
    if(any(is.na(out_y_idx)))
      warning("deSolve/ lsoda cannot solve the ODE system!
              This might be due to arising unrealistic parameters by means of Morris
              Screening. Use ODEsobol() instead or set binf and bsup differently!")
    
    return(out_y_idx)
  }
  
#   cl <- makeCluster(rep("localhost", ncores), type = "SOCK")
#   clusterSetRNGStream(cl)
#   clusterExport(cl, list("mod", "times", "timesNum", "pars",
#                          "yini", "z", "r", "design", "one_y",
#                          "morris_ats", "ode", "k", "out_ats", "response_ats",
#                          "tell_ats", "random.oat", "ind.rep", 
#                          "haussdorf.distance", "kennard.stone",
#                          "morris.maximin"), envir = environment())
#   res <- parLapply(cl, 1:z, one_y)
#   stopCluster(cl)
  
#   out_all <- lapply(1:z, one_y)
#   
#   out_all <- list(y1 = one_y(1))
  out_1 <- one_y(1)
  
  # Rueckgabe:
#   res <- list(res = out_all, pars = pars)
#   class(res) <- "morrisRes_ats"
  
  res <- list(res = out_1, pars = pars)
  class(res) <- "morrisRes"
  return(res)
}
