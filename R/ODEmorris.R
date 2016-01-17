#' @title Morris SA for Objects of Class \code{ODEnetwork}
#'
#' @description
#' \code{ODEmorris} performs a sensitivity analysis for objects of class
#' \code{ODEnetwork} using Morris's elementary effects screening method. Package
#' \code{ODEnetwork} is required for this function to work.
#'
#' @param odenet [\code{ODEnetwork}]\cr
#'   list of class \code{ODEnetwork}.
#' @param pars [\code{character(k)}]\cr
#'   vector of \code{k} input variable names. All parameters must be 
#'   contained in \code{names(ODEnetwork::createParamVec(odenet))} and must not
#'   be derivable from other parameters supplied (e.g., \code{"k.2.1"} can be 
#'   derived from \code{"k.1.2"}, so supplying \code{"k.1.2"} suffices).
#' @param times [\code{numeric}]\cr
#'   points of time at which the SA should be executed (vector of arbitrary 
#'   length). The first point of time must be greater than zero.
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
#'   one state variable (all components of the 2 state variables are analyzed
#'   independently). The matrices themselves contain in their rows the Morris SA
#'   results (i.e. \code{mu, mu.star} and \code{sigma} for every parameter) for
#'   all timepoints (columns).
#'
#' @details The sensitivity analysis is done for all state variables, since 
#' \code{ODEmorris} is an adapted version of \code{\link{ODEmorris_aos}}.
#' 
#' @note \code{\link[deSolve]{ode}} sometimes cannot solve an ODE system if 
#'   unrealistic parameter combinations are sampled by 
#'   \code{\link[sensitivity]{morris_list}}. Hence \code{NA}s might occur in the 
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
#' 
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
#' ODEpars <- c("m.1", "d.1", "k.1", "k.1.2", "m.2", "d.2", "k.2")
#' ODEtimes <- seq(0.01, 20, by = 0.1)
#' ODEbinf <- rep(0.001, length(ODEpars))
#' ODEbsup <- c(2, 1.5, 6, 6, 2, 1.5, 6)
#' 
#' ODEres <- ODEmorris(odenet, ODEpars, ODEtimes, ode_method = "adams", 
#'                     seed = 2015, binf = ODEbinf, bsup = ODEbsup, r = 20)
#'
#' @import checkmate
#' @importFrom deSolve ode
#' @importFrom sensitivity morris_list
#' @export
#'

ODEmorris <- function(odenet,
                      pars,
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
                                 pars,
                                 times, 
                                 ode_method = "lsoda",
                                 seed = 2015,
                                 binf = 0,
                                 bsup = 1,
                                 r = 25,
                                 design = list(type = "oat", levels = 100, 
                                               grid.jump = 1),
                                 scale = TRUE){
  
  ##### Package checks #################################################
  
  if(!requireNamespace("ODEnetwork", quietly = TRUE)){
    stop(paste("Package \"ODEnetwork\" needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }
  
  ##### Input checks ###################################################
  
  assertClass(odenet, "ODEnetwork")
  assertCharacter(pars)
  stopifnot(all(pars %in% names(ODEnetwork::createParamVec(odenet))))
  # Check if there are duplicated parameters:
  if(any(duplicated(pars))){
    pars <- unique(pars)
    warning("Duplicated parameter names in \"pars\". Only taking unique names.")
  }
  # Check if there are parameters which can be derived from others (like 
  # "k.2.1" from "k.1.2"):
  pars_offdiag <- pars[nchar(pars) == 5]
  pars_offdiag_exchanged <- pars_offdiag
  # Exchange third and fifth position:
  substr(pars_offdiag_exchanged, 3, 3) <- substr(pars_offdiag, 5, 5)
  substr(pars_offdiag_exchanged, 5, 5) <- substr(pars_offdiag, 3, 3)
  if(any(pars_offdiag_exchanged %in% pars)){
    pars_deriv <- pars_offdiag[pars_offdiag %in% pars_offdiag_exchanged]
    pars_deriv_exchanged <- 
      pars_offdiag_exchanged[pars_offdiag_exchanged %in% pars_offdiag]
    pars_keep <- character(length(pars_deriv))
    for(i in seq_along(pars_deriv)){
      if(!pars_deriv_exchanged[i] %in% pars_keep){
        pars_keep[i] <- pars_deriv[i]
      } else{
        pars_keep[i] <- "drop"
      }
    }
    pars_keep <- c(pars[nchar(pars) != 5], 
                   pars_offdiag[!pars_offdiag %in% pars_deriv], 
                   pars_keep[pars_keep != "drop"])
    rfuncs <- rfuncs[pars %in% pars_keep]
    rargs <- rargs[pars %in% pars_keep]
    pars <- pars[pars %in% pars_keep]
    warning(paste("Derivable parameters in \"pars\". Keeping only one", 
                  "parameter of each derivable pair."))
  }
  assertNumeric(times, lower = 0, finite = TRUE, unique = TRUE)
  times <- sort(times)
  stopifnot(!any(times == 0))
  stopifnot(ode_method %in% c("lsoda", "lsode", "lsodes","lsodar","vode", 
                              "daspk", "euler", "rk4", "ode23", "ode45", 
                              "radau", "bdf", "bdf_d", "adams", "impAdams", 
                              "impAdams_d" ,"iteration"))
  assertNumeric(seed)
  assertNumeric(binf)
  if(length(binf) != length(pars) && length(binf) != 1)
    stop("binf must be of length 1 or of the same length as pars!")
  assertNumeric(bsup)
  if(length(bsup) != length(pars) & length(bsup) != 1)
    stop("bsup must be of length 1 or of the same length as pars!")
  assertIntegerish(r, len = 1)
  if(r < 1)
    stop("r must be greater or equal to 1.")
  assertList(design)
  assertLogical(scale, len = 1)
  
  ##### Preparation ####################################################
  
  set.seed(seed)
  yini <- ODEnetwork::createState(odenet)
  # Number of parameters:
  k <- length(pars)
  # Number of output variables (state variables):
  z <- length(yini)
  # Number of timepoints:
  timesNum <- length(times)
  
  # Adapt the ODE-model for argument "model" of soboljansen_list() resp.
  # sobolmartinez_list():
  model_fit <- function(X){
    # Input: matrix X with k columns
    colnames(X) <- pars
    res_per_par <- lapply(1:nrow(X), function(i){
      pars_upd <- X[i, ]
      names(pars_upd) <- pars
      odenet_parmod <- ODEnetwork::updateOscillators(odenet, 
                                                     ParamVec = pars_upd)
      ODEnetwork::simuNetwork(odenet_parmod, c(0, times), 
        method = ode_method)$simulation$results[2:(timesNum + 1), 2:(z + 1)]
    })
    if(timesNum == 1){
      # Correction needed, if timesNum == 1:
      res_vec <- unlist(res_per_par)
      res_matrix <- matrix(res_vec, ncol = 1)
    } else{
      # Transpose the matrix of the results, so that each column represents
      # one timepoint:
      res_matrix <- t(do.call(cbind, res_per_par))
    }
    rownames(res_matrix) <- NULL
    nrow_res_matrix <- nrow(res_matrix)
    res_per_y <- lapply(1:z, function(i){
      res_matrix[seq(i, nrow_res_matrix, z), , drop = FALSE]
    })
    names(res_per_y) <- names(yini)
    return(res_per_y)
  }
  
  ##### Sensitivity analysis #########################################
  
  # Sensitivity analysis with function morris_list() from package "sensitivity":
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
  
  # Throw a warning if NAs occur (probably not suitable parameters, so ODE
  # system can't be solved):
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
  
  # Return:
  class(out_all_y) <- "morrisRes"
  return(out_all_y)
}