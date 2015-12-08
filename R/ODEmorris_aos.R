#' @title Morris SA for ODEs for All Output Variables and All Timepoints
#' Simultaneously
#'
#' @description
#' \code{ODEmorris_aos} performs a sensitivity analysis for
#' ordinary differential equations using Morris's elementary effects 
#' screening method. The analysis is done for all output variables and all
#' timepoints simultaneously using \code{\link[sensitivity]{morris_list}}.
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
#' @param scale [\code{logical(1)}]\cr
#'   if \code{TRUE}, scaling is done for the input design of experiments after 
#'   building the design and before calculating the elementary effects,
#'   cf. \code{\link[sensitivity]{morris}}. Defaults to \code{FALSE}, but it is
#'   highly recommended to use \code{scale = TRUE} if the factors have different
#'   orders of magnitude.
#'
#' @return list of class \code{morrisRes_aos} with
#'   \itemize{
#'     \item \code{res}, a list containing in each element a matrix for one
#'       output variable. The matrices itself contain the Morris SA results 
#'       (i.e. \code{mu, mu.star} and \code{sigma} for every parameter) for one 
#'       point of time per column
#'     \item \code{pars}, the parameter names.
#'   }
#'
#' @details \code{ODEmorris_aos} uses \code{\link[sensitivity]{morris_list}}
#' which can handle \code{list}s as output for its model function. Thus, each
#' element of the list can be used to contain the results for one output
#' variable. This saves time since \code{\link[deSolve]{ode}} from the
#' package \code{deSolve} does its calculations for all output variables anyway.
#' This way, \code{\link[deSolve]{ode}} only runs once and the rest is
#' reformatting the data to a \code{list} output.
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
#' FHNpars  <- c(a = 0.2,     # parameter a
#'               b = 0.3,     # parameter b
#'               s = 3)       # parameter s (= c in the original notation)
#'
#' FHNyini  <- c(Voltage = -1, Current = 1)
#' FHNtimes <- seq(0.1, 100, by = 10)
#'
#' FHNres_aos <- ODEmorris_aos(mod = FHNmod,
#'                             pars = names(FHNpars),
#'                             yini = FHNyini,
#'                             times = FHNtimes,
#'                             seed = 2015,
#'                             binf = c(0.18, 0.18, 2.8),
#'                             bsup = c(0.22, 0.22, 3.2),
#'                             r = 25,
#'                             design =
#'                                list(type = "oat", levels = 100, 
#'                                     grid.jump = 1),
#'                             scale = TRUE)
#'
#' @seealso \code{\link[sensitivity]{morris}},
#'   \code{\link[sensitivity]{morris_list}},
#'   \code{\link{plot.morrisRes_aos}}
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
#'   If \code{\link[sensitivity]{morris_list}} throws a warning message saying
#'   "In ... keeping ... repetitions out of ...", try using a bigger number of 
#'   \code{levels} in the \code{design} argument.
#'
#' @references J. O. Ramsay, G. Hooker, D. Campbell and J. Cao, 2007,
#'   \emph{Parameter estimation for differential equations: a generalized 
#'   smoothing approach}, Journal of the Royal Statistical Society, Series B, 
#'   69, Part 5, 741--796.
#'   
#'
#' @export
#' @import
#'   checkmate
#'   deSolve
#'   sensitivity
#'

ODEmorris_aos <- function(mod,
                          pars,
                          yini,
                          times,
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
  
  ##### Sensitivitaet ##################################################
  
  # Forme DGL-Modell um, sodass fuer morris_list()-Argument "model_list" 
  # passend:
  model_fit <- function(X){
    # X   - (nxk)-Matrix mit den n einzugebenden Parameter-Konstellationen
    #       als Zeilen
    colnames(X) <- pars
    res_per_par <- lapply(1:nrow(X), function(i){
      ode(yini, times = c(0, times), 
          mod, parms = X[i, ])[2:(timesNum + 1), 2:(z + 1)]
    })
    if(timesNum == 1){
      # Correction needed if timesNum == 1:
      res_vec <- unlist(res_per_par)
      res_matrix <- matrix(res_vec, ncol = 1)
    } else{
      # Transpose the results matrices, so columns represent timepoints:
      res_matrix <- t(do.call(cbind, res_per_par))
    }
    # Entferne verwirrende Zeilennamen:
    rownames(res_matrix) <- NULL
    nrow_res_matrix <- nrow(res_matrix)
    res_per_y <- lapply(1:z, function(i){
      res_matrix[seq(i, nrow_res_matrix, z), , drop = FALSE]
    })
    names(res_per_y) <- names(yini)
    return(res_per_y)
  }
  
  x <- morris_list(model_list = model_fit, factors = pars, r = r, 
                   design = design, binf = binf, bsup = bsup, scale = scale)
  
  one_y <- function(L){
    mu <- lapply(L, colMeans)
    mu.star <- lapply(L, abs)
    mu.star <- lapply(mu.star, colMeans)
    sigma <- lapply(L, function(M){
      apply(M, 2, sd)
    })
    out_y <- mapply(c, mu, mu.star, sigma, SIMPLIFY = TRUE)
    out_y <- rbind(times, out_y)
    rownames(out_y) <- c("time", paste0("mu_", pars), 
                         paste0("mu.star_", pars),
                         paste0("sigma_", pars))
    return(out_y)
  }
  
  out_all_y <- lapply(x$ee_by_y, one_y)
  
  # Warnungen, falls NAs auftreten (unrealistische Parameter => nicht
  # loesbare ODEs):
  NA_check_mu <- function(M){
    any(is.na(M[1:(1 + k*2), ]))
  }
  NA_check_sigma <- function(M){
    all(is.na(M[(2 + k*2):(1 + k*3), ]))
  }
  if(any(unlist(lapply(out_all_y, NA_check_mu)))){
    warning("deSolve/ lsoda cannot solve the ODE system!
            This might be due to arising unrealistic parameters by means of 
            Morris Screening. Use ODEsobol() instead or set binf and bsup 
            differently!")
  } else if(all(unlist(lapply(out_all_y, NA_check_sigma))) && r == 1){
    warning("Calculation of sigma requires r >= 2.")
  } else{
    NA_check_sigma_any <- function(M){
      any(is.na(M[(2 + k*2):(1 + k*3), ]))
    }
    if(any(unlist(lapply(out_all_y, NA_check_sigma_any)))){
      warning("NAs for sigma. This might be due to r being too small.")
    }
  }
  
  res <- list(res = out_all_y, pars = pars)
  class(res) <- "morrisRes_aos"
  return(res)
}
