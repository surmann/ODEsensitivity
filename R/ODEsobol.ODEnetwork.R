#' @title Sobol' SA for Objects of Class \code{ODEnetwork}
#'
#' @description
#' \code{ODEsobol.ODEnetwork} performs a variance-based sensitivity analysis 
#' for objects of class \code{ODEnetwork} according to either the Sobol'-Jansen-
#' or the Sobol'-Martinez-method. Package \code{ODEnetwork} is required for this 
#' function to work.
#'
#' @param mod [\code{ODEnetwork}]\cr
#'   list of class \code{ODEnetwork}.
#' @param pars [\code{character(k)}]\cr
#'   vector of \code{k} input variable names. All parameters must be 
#'   contained in \code{names(ODEnetwork::createParamVec(mod))} and must not
#'   be derivable from other parameters supplied (e.g., \code{"k.2.1"} can be 
#'   derived from \code{"k.1.2"}, so supplying \code{"k.1.2"} suffices).
#' @param times [\code{numeric}]\cr
#'   points of time at which the SA should be executed (vector of arbitrary 
#'   length). The first point of time must be greater than zero.
#' @param seed [\code{numeric(1)}]\cr
#'   seed.
#' @param n [\code{integer(1)}]\cr
#'   number of random parameter values (\code{n} per input factor) used to 
#'   estimate the variance-based sensitivity indices by Monte Carlo method.
#'   (Variance-based methods for sensitivity analysis rely on 
#'   Monte Carlo simulation to estimate the integrals needed for the calculation
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
#' @param sobol_method [\code{character(1)}]\cr
#'   either \code{"jansen"} or \code{"martinez"}, specifying which modification
#'   of the variance-based Sobol' method shall be used. Defaults to 
#'   \code{"martinez"}.
#' @param ode_method [\code{character(1)}]\cr
#'   method to be used for solving the differential equations, see 
#'   \code{\link[deSolve]{ode}}. Defaults to \code{"lsoda"}.
#' @param ode_parallel [\code{logical(1)}]\cr
#'   logical indicating if a parallelization shall be done for computing the
#'   \code{\link[deSolve]{ode}}-results for the different parameter combinations
#'   generated for Monte Carlo estimation of the sensitivity indices.
#' @param ode_parallel_ncores [\code{integer(1)}]\cr
#'   number of processor cores to be used for parallelization. Only applies if
#'   \code{ode_parallel = TRUE}. If set to \code{NA} (as per default) and 
#'   \code{ode_parallel = TRUE}, 1 processor core is used.
#' @param ... further arguments passed to or from other methods.
#'
#' @return 
#'   List of length \code{2 * nrow(mod$state)} and of class 
#'   \code{sobolRes} containing in each element a list of the Sobol' SA results 
#'   for the corresponding state variable (i.e. first order sensitivity indices
#'   \code{S} and total sensitivity indices \code{T}) for every point of time of 
#'   the \code{times} vector.
#'
#' @details 
#'   The sensitivity analysis is done for all state variables since the 
#'   simulation of the network implies all state variables anyway.
#'
#' @note 
#'   It might be helpful to try different types of ODE-solvers (argument 
#'   \code{ode_method}). Problems are known for the \code{ode_method}s 
#'   \code{"euler"}, \code{"rk4"} and \code{"ode45"}. 
#'   In contrast, the \code{ode_method}s \code{"vode"}, \code{"bdf"}, 
#'   \code{"bdf_d"}, \code{"adams"}, \code{"impAdams"} and \code{"impAdams_d"} 
#'   might be even faster than the standard \code{ode_method} \code{"lsoda"}.
#'   
#'   If \code{n} is too low, the Monte Carlo estimation of the sensitivity 
#'   indices might be very bad and even produce first order indices < 0 or
#'   total indices > 1. First order indices in the interval [-0.05, 0) and total 
#'   indices in (1, 1.05] are considered as minor deviations and set to 0 
#'   resp. 1 without a warning. First order indices < -0.05 or total indices 
#'   > 1.05 are considered as major deviations. They remain unchanged and a 
#'   warning is thrown. Up to now, first order indices > 1 or total indices < 0
#'   haven't occured yet. If this should be the case, please contact the package
#'   author.
#'
#' @author Frank Weber
#' @seealso \code{\link[sensitivity]{soboljansen}, 
#'   \link[sensitivity]{sobolmartinez},
#'   \link{plot.sobolRes}}
#' 
#' @examples
#' ##### A network of ordinary differential equations #####
#' # Definition of the network using the package "ODEnetwork":
#' library(ODEnetwork)
#' masses <- c(1, 1)
#' dampers <- diag(c(1, 1))
#' springs <- diag(c(1, 1))
#' springs[1, 2] <- 1
#' distances <- diag(c(0, 2))
#' distances[1, 2] <- 1
#' lfonet <- ODEnetwork(masses, dampers, springs, 
#'                      cartesian = TRUE, distances = distances)
#' lfonet <- setState(lfonet, c(0.5, 1), c(0, 0))
#' LFOpars <- c("m.1", "d.1", "k.1", "k.1.2", "m.2", "d.2", "k.2")
#' LFOtimes <- seq(0.01, 20, by = 0.1)
#' LFObinf <- rep(0.001, length(LFOpars))
#' LFObsup <- c(2, 1.5, 6, 6, 2, 1.5, 6)
#' 
#' LFOres <- ODEsobol(lfonet, 
#'                    LFOpars, 
#'                    LFOtimes, 
#'                    seed = 2015, 
#'                    n = 1000,
#'                    rfuncs = c(rep("runif", length(LFObinf)), "rnorm"),
#'                    rargs = c(paste0("min = ", LFObinf, ", max = ", LFObsup),
#'                              "mean = 3, sd = 0.8"),
#'                    sobol_method = "martinez",
#'                    ode_method = "adams",
#'                    ode_parallel = TRUE,
#'                    ode_parallel_ncores = 2)
#' 
#' @import checkmate
#' @importFrom deSolve ode
#' @importFrom sensitivity soboljansen
#' @importFrom sensitivity sobolmartinez
#' @method ODEsobol ODEnetwork
#' @export
#' 

ODEsobol.ODEnetwork <- function(mod,
                                pars,
                                times,
                                seed = 2015,
                                n = 1000,
                                rfuncs = rep("runif", length(pars)),
                                rargs = rep("min = 0, max = 1", length(pars)),
                                sobol_method = "martinez",
                                ode_method = "lsoda",
                                ode_parallel = FALSE,
                                ode_parallel_ncores = NA, ...){
  
  ##### Package checks #################################################
  
  if(!requireNamespace("ODEnetwork", quietly = TRUE)){
    stop(paste("Package \"ODEnetwork\" needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }
  
  ##### Input checks ###################################################
  
  assertClass(mod, "ODEnetwork")
  assertCharacter(pars)
  stopifnot(all(pars %in% names(ODEnetwork::createParamVec(mod))))
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
  assertNumeric(seed)
  assertIntegerish(n)
  assertCharacter(rfuncs, len = length(pars))
  assertCharacter(rargs, len = length(pars))
  rfuncs_exist <- sapply(rfuncs, exists)
  if(!all(rfuncs_exist)) stop(paste("At least one of the supplied functions",
                                    "in \"rfuncs\" was not found"))
  stopifnot(sobol_method %in% c("jansen", "martinez"))
  stopifnot(ode_method %in% c("lsoda", "lsode", "lsodes","lsodar","vode", 
                              "daspk", "euler", "rk4", "ode23", "ode45", 
                              "radau", "bdf", "bdf_d", "adams", "impAdams", 
                              "impAdams_d" ,"iteration"))
  assertLogical(ode_parallel, len = 1)
  assertIntegerish(ode_parallel_ncores, len = 1, lower = 1)
  if(ode_parallel && is.na(ode_parallel_ncores)){
    ode_parallel_ncores <- 1
  }
  
  ##### Preparation ####################################################
  
  set.seed(seed)
  # Initial state:
  state_init <- ODEnetwork::createState(mod)
  # Number of parameters:
  k <- length(pars)
  # Number of state variables:
  z <- length(state_init)
  # Number of timepoints:
  timesNum <- length(times)
  
  # Adapt the ODE-model for argument "model" of soboljansen() resp.
  # sobolmartinez():
  model_fit <- function(X){
    # Input: Matrix X with k columns, containing the random parameter 
    # combinations.
    colnames(X) <- pars
    one_par <- function(i){
      # Get the parameter values:
      pars_upd <- X[i, ]
      names(pars_upd) <- pars
      # Update the model function in the ODEnetwork-object "mod":
      mod_parmod <- ODEnetwork::updateOscillators(mod,
                                                  ParamVec = pars_upd)
      # Simulate the network (the results correspond to those by ode() in the 
      # package "deSolve"):
      simnet_res <- ODEnetwork::simuNetwork(mod_parmod, 
                                            c(0, times), 
                                            method = ode_method)
      return(simnet_res$simulation$results[2:(timesNum + 1), 2:(z + 1)])
    }
    if(ode_parallel){
      # Run one_par() on parallel nodes:
      ode_cl <- parallel::makeCluster(rep("localhost", ode_parallel_ncores), 
                                      type = "PSOCK")
      parallel::clusterExport(ode_cl, 
                              varlist = c("X", "pars", "mod", "z",
                                          "times", "ode_method", "timesNum"),
                              envir = environment())
      res_per_par <- parallel::parLapply(ode_cl, 1:nrow(X), one_par)
      parallel::stopCluster(ode_cl)
    } else{
      # Just use lapply():
      res_per_par <- lapply(1:nrow(X), one_par)
    }
    if(timesNum == 1){
      # Correction needed if timesNum == 1:
      res_vec <- unlist(res_per_par)
      res_matrix <- matrix(res_vec, ncol = 1)
    } else{
      # Transpose the matrix of the results, so that each column represents
      # one timepoint:
      res_matrix <- t(do.call(cbind, res_per_par))
    }
    rownames(res_matrix) <- NULL
    # Convert the results matrix to a list (one element for each state
    # variable):
    nrow_res_matrix <- nrow(res_matrix)
    res_per_state <- lapply(1:z, function(i){
      res_matrix[seq(i, nrow_res_matrix, z), , drop = FALSE]
    })
    names(res_per_state) <- names(state_init)
    return(res_per_state)
  }
  
  # Create the two matrices containing the parameter samples for Monte Carlo
  # estimation:
  rfunc_calls <- paste0(rfuncs, "(n, ", rargs, ")", collapse = ", ")
  X1 <- matrix(eval(parse(text = paste0("c(", rfunc_calls, ")"))), ncol = k)
  X2 <- matrix(eval(parse(text = paste0("c(", rfunc_calls, ")"))), ncol = k)
  colnames(X1) <- colnames(X2) <- pars
  
  ##### Sensitivity analysis #########################################
  
  # Sensitivity analysis with either soboljansen_list() or sobolmartinez_list()
  # from package "sensitivity":
  if(sobol_method == "jansen"){
    x <- soboljansen_list(model = model_fit, X1, X2, nboot = 0)
  } else if(sobol_method == "martinez"){
    x <- sobolmartinez_list(model = model_fit, X1, X2, nboot = 0)
  }
  
  # Process the results:
  ST_by_state <- sobol_process(x, pars, times)
  
  # Return:
  res <- list(ST_by_state = ST_by_state, sobol_method = sobol_method)
  class(res) <- "sobolRes"
  return(res)
}
