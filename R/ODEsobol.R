#' @title Sobol' SA for Objects of Class \code{ODEnetwork}
#'
#' @description
#' \code{ODEsobol} performs a variance-based sensitivity analysis for objects of
#' class \code{ODEnetwork} according to either the Sobol'-Jansen- or the 
#' Sobol'-Martinez-method. Package \code{ODEnetwork} is required for this 
#' function to work.
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
#' @param n [\code{integer(1)}]\cr
#'   number of random parameter values (\code{n} per input factor) used to 
#'   estimate the variance-based sensitivity indices by Monte-Carlo-method.
#'   (Variance-based methods for sensitivity analysis rely on 
#'   Monte-Carlo-simulation to estimate the integrals needed for the calculation
#'   of the sensitivity indices.) Defaults to 1000.
#' @param rfuncs [\code{character(k)}]\cr
#'   names of the \code{k} functions used to generate the \code{n} random values
#'   for the \code{k} parameters (see details). This way, different 
#'   distributions can be supplied for the \code{k} parameters. Defaults to 
#'   \code{"runif"} for each of the \code{k} parameters.
#' @param rargs [\code{character(k)}]\cr
#'   arguments to be passed to the \code{k} functions of \code{rfuncs}. Each 
#'   element of \code{rargs} has to be a string of type \code{"tag1 = value1, 
#'   tag2 = value2, ..."}. By default, \code{min = 0} and \code{max = 1} are 
#'   used for each of the \code{k} \code{runif}'s, meaning a uniform 
#'   distribution of all parameters on [0, 1].
#' @param method [\code{character(1)}]\cr
#'   either \code{"jansen"} or \code{"martinez"}, specifying which modification
#'   of the variance-based Sobol' method shall be used. Defaults to 
#'   \code{"martinez"}, which is slightly faster than \code{"jansen"}.
#' @param nboot [\code{integer(1)}]\cr
#'   parameter \code{nboot} used in \code{\link{soboljansen_list}} resp.
#'   \code{\link{sobolmartinez_list}}, i.e. the number of bootstrap 
#'   replicates. Defaults to 0, so no bootstrapping is done.
#'
#' @return List of length \code{2 * nrow(odenet$state)} and of class 
#' \code{sobolRes} containing in each element a list of the Sobol' SA results 
#' for the corresponding state-variable (i.e. first order sensitivity indices
#' \code{S} and total sensitivity indices \code{T}) for every point of time of 
#' the \code{times} vector.
#'
#' @details The sensitivity analysis is done for all state-variables, since 
#' \code{ODEsobol} is an adapted version of \code{\link{ODEsobol_aos}}.
#'
#' @note 
#'   Sometimes, it is also helpful to try another ODE-solver (argument 
#'   \code{ode_method}). Problems are known for the
#'   \code{ode_method}s \code{"euler"}, \code{"rk4"} and \code{"ode45"}. 
#'   In contrast, the \code{ode_method}s \code{"vode"}, \code{"bdf"}, 
#'   \code{"bdf_d"}, \code{"adams"}, \code{"impAdams"} and \code{"impAdams_d"} 
#'   might be even faster than the standard \code{ode_method} \code{"lsoda"}.
#'
#' @author Frank Weber
#' @references J. O. Ramsay, G. Hooker, D. Campbell and J. Cao, 2007,
#'   \emph{Parameter estimation for differential equations: a generalized 
#'   smoothing approach}, Journal of the Royal Statistical Society, Series B, 
#'   69, Part 5, 741--796.
#' @seealso \code{\link[sensitivity]{soboljansen_list},
#' \link[sensitivity]{sobolmartinez_list},
#' \link{plot.sobolRes_aos}}
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
#' ODEbinf <- rep(0.001, length(ODEpars) - 1)
#' ODEbsup <- c(2, 1.5, 6, 6, 2, 1.5)
#' 
#' ODEres <- ODEsobol(odenet, ODEpars, ODEtimes, ode_method = "adams", 
#'                    seed = 2015, n = 10,
#'                    rfuncs = c(rep("runif", length(ODEbinf)), "rnorm"),
#'                    rargs = c(paste0("min = ", ODEbinf, ", max = ", ODEbsup),
#'                              "mean = 3, sd = 0.8"),
#'                    method = "martinez",
#'                    nboot = 0)
#'
#' @import checkmate
#' @importFrom deSolve ode
#' @importFrom sensitivity soboljansen_list
#' @importFrom sensitivity sobolmartinez_list
#' @export
#'

ODEsobol <- function(odenet,
                     pars,
                     times,
                     ode_method = "lsoda",
                     seed = 2015,
                     n = 1000,
                     rfuncs = rep("runif", length(pars)),
                     rargs = rep("min = 0, max = 1", length(pars)),
                     method = "martinez",
                     nboot = 0) {
  UseMethod("ODEsobol", odenet)
}

#' @method ODEsobol ODEnetwork
#' @export

ODEsobol.ODEnetwork <- function(odenet,
                                pars,
                                times,
                                ode_method = "lsoda",
                                seed = 2015,
                                n = 1000,
                                rfuncs = rep("runif", length(pars)),
                                rargs = rep("min = 0, max = 1", length(pars)),
                                method = "martinez",
                                nboot = 0){
  
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
    rfuncs <- rfuncs[!duplicated(pars)]
    rargs <- rargs[!duplicated(pars)]
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
  assertIntegerish(n)
  assertCharacter(rfuncs, len = length(pars))
  assertCharacter(rargs, len = length(pars))
  rfuncs_exist <- sapply(rfuncs, exists)
  if(!all(rfuncs_exist)) stop(paste("At least one of the supplied functions",
                                    "in \"rfuncs\" was not found"))
  stopifnot(method %in% c("jansen", "martinez"))
  assertIntegerish(nboot)
  
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
  
  rfunc_calls <- paste0(rfuncs, "(n, ", rargs, ")", collapse = ", ")
  X1 <- matrix(eval(parse(text = paste0("c(", rfunc_calls, ")"))), ncol = k)
  X2 <- matrix(eval(parse(text = paste0("c(", rfunc_calls, ")"))), ncol = k)
  colnames(X1) <- colnames(X2) <- pars
  
  ##### Sensitivity analysis #########################################
  
  # Sensitivity analysis with the functions from package "sensitivity":
  if(method == "jansen"){
    x <- soboljansen_list(model = model_fit, X1, X2, nboot = nboot)
  } else if(method == "martinez"){
    x <- sobolmartinez_list(model = model_fit, X1, X2, nboot = nboot)
  }
  
  # Process the results:
  ST_original_by_y <- lapply(x$ST_by_y, function(L){
    ST_original <- sapply(L, function(ST_col){
      c(ST_col$S[, 1], ST_col$T[, 1])
    })
    return(ST_original)
  })
  # Split ST_original again in 2 matrices S and T:
  ST_by_y <- lapply(ST_original_by_y, function(ST_original){
    S <- rbind(times, ST_original[1:k, ])
    T <- rbind(times, ST_original[(k+1):(2*k), ])
    rownames(S) <- rownames(T) <- c("time", pars)
    return(list(S = S, T = T))
  })
  
  # Return:
  res <- list(ST_by_y = ST_by_y, method = method)
  class(res) <- "sobolRes"
  return(res)
}
