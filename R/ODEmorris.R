#' @title Morris SA for objects of class \code{ODEnetwork}
#'
#' @description
#' \code{ODEmorris} performs a sensitivity analysis for objects of class
#' \code{ODEnetwork} using Morris's elementary effects screening method. Package
#' \code{ODEnetwork} is required for this function to work. 
#'
#' @param odenet [\code{ODEnetwork}]\cr
#'   list of class \code{ODEnetwork}.
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
#'   vector of lower borders of possible values for the \code{k} input 
#'   parameters. If they are all equal, a single value can be set.
#' @param bsup [\code{numeric(k)}]\cr
#'   vector of upper borders of possible values for the \code{k} input 
#'   parameters. If they are all equal, a single value can be set.
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
#' @return List of class \code{morrisRes} of length 
#'   \code{2 * nrow(odenet$state)} containing in each element a matrix for 
#'   one state variable. The matrices itself contain in their rows the Morris SA
#'   results (i.e. \code{mu, mu.star} and \code{sigma} for every parameter) for
#'   all timepoints (columns).
#'
#' @details The sensitivity analysis is done for all output variables, since 
#' \code{ODEmorris} is an adapted version of \code{\link{ODEmorris_aos}}.
#' 
#' @note \code{\link[deSolve]{ode}} sometimes cannot solve an ODE system if 
#'   unrealistic parameter combinations are sampled by 
#'   \code{\link[sensitivity]{morris_list}}. Hence \code{NA}s might occur in the 
#'   Morris sensitivity results, such that \code{\link{ODEmorris_aos}} fails for 
#'   one or many points of time! For this reason, if \code{NA}s occur, please 
#'   make use of \code{\link{ODEsobol_ats}} instead or
#'   restrict the input parameter value intervals usefully using
#'   \code{binf}, \code{bsup} and \code{scale = TRUE}. It is also helpful to try
#'   another ODE-solver (argument \code{ode_method}). Problems are known for the
#'   \code{ode_method}s \code{"euler"}, \code{"rk4"} and \code{"ode45"}. 
#'   In contrast, the \code{ode_method}s \code{"vode"}, \code{"bdf"}, 
#'   \code{"bdf_d"}, \code{"adams"}, \code{"impAdams"} and \code{"impAdams_d"} 
#'   might be even faster than the standard \code{ode_method} \code{"lsoda"}.
#'   
#'   If \code{\link[sensitivity]{morris_list}} throws a warning message saying
#'   "In ... keeping ... repetitions out of ...", try using a bigger number of 
#'   \code{levels} in the \code{design} argument.
#'
#' @author Frank Weber
#' @seealso \code{\link[sensitivity]{morris}},
#'   \code{\link[sensitivity]{morris_list}},
#'   \code{\link{plot.morrisRes}}
#' 
#' @examples
#' library(ODEnetwork)
#' masses <- c(1, 1)
#' dampers <- diag(c(1, 1))
#' springs <- diag(c(1, 1))
#' springs[1, 2] <- 1
#' distances <- diag(c(0, 2))
#' distances[1, 2] <- 1
#' odenet <- ODEnetwork(masses, dampers, springs, 
#'                      cartesian = TRUE, distances = distances)
#' odenet <- setState(odenet, c(0.5, 1), c(0, 0))
#' 
#' ODEtimes <- seq(0.01, 20, by = 0.1)
#' ODEbinf <- c(rep(0.001, 9), -4, 0.001, -10)
#' ODEbsup <- c(2, 1.5, 6, 4, 6, 10, 2, 1.5, 6, -0.001, 6, -0.001)
#' 
#' ODEres <- ODEmorris(odenet, ODEtimes, ode_method = "adams", seed = 2015, 
#'                     binf = ODEbinf, bsup = ODEbsup, r = 20)
#'
#' @import
#'   checkmate
#' @importFrom deSolve ode
#' @importFrom sensitivity morris_list
#' @export
#'

ODEmorris <- function(odenet, 
                      times, 
                      ode_method = "lsoda",
                      seed = 2015,
                      binf = 0,
                      bsup = 1,
                      r = 25,
                      design =
                        list(type = "oat", levels = 100, 
                             grid.jump = 1),
                      scale = TRUE){
  UseMethod("ODEmorris", odenet)
}

#' @method ODEmorris ODEnetwork
#' @export

ODEmorris.ODEnetwork <- function(odenet, 
                                 times, 
                                 ode_method = "lsoda",
                                 seed = 2015,
                                 binf = 0,
                                 bsup = 1,
                                 r = 25,
                                 design =
                                   list(type = "oat", levels = 100, 
                                 grid.jump = 1),
                                 scale = TRUE){
  ##### Package checks #################################################
  
  if (!requireNamespace("ODEnetwork", quietly = TRUE)) {
    stop(paste("Package \"ODEnetwork\" needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }
  
  ##### Input checks #################################################
  
  set.seed(seed)
  odenet_pars <- ODEnetwork::createParamVec(odenet)
  pars <- names(odenet_pars)
  yini <- ODEnetwork::createState(odenet)
  # Anzahl Parameter:
  k <- length(pars)
  # Anzahl Outputgroessen:
  z <- length(yini)
  # Anzahl Zeitpunkte von Interesse:
  timesNum <- length(times)
  
  ##### Sensitivity analysis #########################################
  
  # Forme DGL-Modell um, sodass fuer morris_list()-Argument "model_list" 
  # passend:
  model_fit <- function(X){
    # X   - (nxk)-Matrix mit den n einzugebenden Parameter-Konstellationen
    #       als Zeilen
    colnames(X) <- pars
    res_per_par <- lapply(1:nrow(X), function(i){
      pars_upd <- X[i, ]
      names(pars_upd) <- pars
      odenet_parmod <- ODEnetwork::updateOscillators(odenet, 
                                                     ParamVec = pars_upd)
      ODEnetwork::simuNetwork(odenet_parmod, c(0, times), 
                  method = ode_method)$simulation$results[2:(timesNum + 1), 
                                                          2:(z + 1)]
    })
    if(timesNum == 1){
      # Korrektur noetig, falls timesNum == 1:
      res_vec <- unlist(res_per_par)
      res_matrix <- matrix(res_vec, ncol = 1)
    } else{
      # Transponiere die Ergebnis-Matrix, sodass jede Spalte fuer einen
      # Zeitpunkt steht:
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
  
  x <- morris_list(model = model_fit, factors = pars, r = r, 
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
    warning(paste("The ODE system can't be solved. This might be due to", 
                  "arising unrealistic parameters by means of Morris Screening. Use",
                  "ODEsobol() instead or set binf and bsup differently together with",
                  "scale = TRUE. It might also be helpful to try another ODE-solver by",
                  "using the \"ode_method\"-argument."))
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
  
  class(out_all_y) <- "morrisRes"
  return(out_all_y)
}