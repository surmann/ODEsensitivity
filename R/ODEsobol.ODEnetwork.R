#' @title Sobol' Sensitivity Analysis for Objects of Class \code{ODEnetwork}
#'
#' @description
#' \code{ODEsobol.ODEnetwork} performs the variance-based Sobol' sensitivity 
#' analysis for objects of class \code{ODEnetwork}. Package 
#' \code{ODEnetwork} is required for this function to work.
#'
#' @param mod [\code{ODEnetwork}]\cr
#'   list of class \code{ODEnetwork}.
#' @param pars [\code{character(k)}]\cr
#'   names of the parameters to be included as input variables in the Sobol'
#'   sensitivity analysis. All parameters must be contained in 
#'   \code{names(ODEnetwork::createParamVec(mod))} and must not be derivable 
#'   from other parameters supplied (e.g. \code{"k.2.1"} can be derived from 
#'   \code{"k.1.2"}, so supplying \code{"k.1.2"} suffices).
#' @param times [\code{numeric}]\cr
#'   points of time at which the sensitivity analysis should be executed (vector
#'   of arbitrary length). The first point of time must be greater than zero.
#' @param n [\code{integer(1)}]\cr
#'   number of random parameter values used to estimate the Sobol' sensitivity 
#'   indices by Monte Carlo simulation. Defaults to 1000.
#' @param rfuncs [\code{character(1} or \code{k)}]\cr
#'   names of the functions used to generate the \code{n} random values
#'   for the \code{k} parameters. Can be of length 1 or \code{k}. If of length 
#'   1, the same function is used for all parameters. Defaults to 
#'   \code{"runif"}, so a uniform distribution is assumed for all parameters.
#' @param rargs [\code{character(1} or \code{k)}]\cr
#'   arguments to be passed to the functions in \code{rfuncs}. Can be of length 
#'   1 or \code{k}. If of length 1, the same arguments are used for all 
#'   parameters. Each element of \code{rargs} has to be a string of the form 
#'   \code{"tag1 = value1, tag2 = value2, ..."}, see example below. Default is 
#'   \code{"min = 0, max = 1"}, so (together with the default value of 
#'   \code{rfuncs}) a uniform distribution on [0, 1] is assumed for all 
#'   parameters.
#' @param sobol_method [\code{character(1)}]\cr
#'   either \code{"Jansen"} or \code{"Martinez"}, specifying which modification
#'   of the variance-based Sobol' method shall be used. Defaults to 
#'   \code{"Martinez"}.
#' @param ode_method [\code{character(1)}]\cr
#'   method to be used for solving the ODEs in situations where the solution has
#'   to be determined numerically, see \code{\link[deSolve]{ode}} for details.
#'   Defaults to \code{"lsoda"}.
#' @param parallel_eval [\code{logical(1)}]\cr
#'   logical indicating if the evaluation of the ODE model shall be performed
#'   parallelized.
#' @param parallel_eval_ncores [\code{integer(1)}]\cr
#'   number of processor cores to be used for parallelization. Only applies if
#'   \code{parallel_eval = TRUE}. If set to \code{NA} (as per default) and 
#'   \code{parallel_eval = TRUE}, 1 processor core is used.
#' @param ... further arguments passed to or from other methods.
#'
#' @return 
#'   List of length \code{2 * nrow(mod$state)} and of class 
#'   \code{ODEsobol} containing in each element a list of the Sobol' sensitivity
#'   analysis results for the corresponding state variable (i.e. first order 
#'   sensitivity indices \code{S} and total sensitivity indices \code{T}) for 
#'   every point of time in the \code{times} vector. This list has an extra 
#'   attribute \code{"sobol_method"} where the value of argument 
#'   \code{sobol_method} is stored (either \code{"Jansen"} or 
#'   \code{"Martinez"}).
#'
#' @details
#'   If the object of class \code{ODEnetwork} supplied for \code{mod} doesn't
#'   include any events, the solution of the ODE network is determined 
#'   analytically using \code{\link[ODEnetwork]{simuNetwork}}. In the presence
#'   of events, \code{\link[ODEnetwork]{simuNetwork}} uses 
#'   \code{\link[deSolve]{ode}} to solve the ODE network numerically.
#' 
#'   The sensitivity analysis is done for all state variables and all
#'   timepoints simultaneously. If \code{sobol_method = "Jansen"},
#'   \code{\link[sensitivity]{soboljansen}} from the package \code{sensitivity}
#'   is used to estimate the Sobol' sensitivity indices and if 
#'   \code{sobol_method = "Martinez"}, \code{\link[sensitivity]{sobolmartinez}}
#'   is used (also from the package \code{sensitivity}).
#'
#' @note 
#'   In situations where the solution of the ODE model has to be determined 
#'   numerically, it might be helpful to try a different type of ODE-solver 
#'   (argument \code{ode_method}) if the simulation of the model takes too long.
#'   The \code{ode_method}s \code{"vode"}, \code{"bdf"}, \code{"bdf_d"}, 
#'   \code{"adams"}, \code{"impAdams"} and \code{"impAdams_d"} 
#'   might be faster than the standard \code{ode_method} \code{"lsoda"}.
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
#'   \link{plot.ODEsobol}}
#' 
#' @examples
#' ##### A network of 4 mechanical oscillators connected in a circle #####
#' # Definition of the network using the package "ODEnetwork":
#' M_mat <- rep(2, 4)
#' K_mat <- diag(rep(2 * (2*pi*0.17)^2, 4))
#' K_mat[1, 2] <- K_mat[2, 3] <- 
#'   K_mat[3, 4] <- K_mat[1, 4] <- 2 * (2*pi*0.17)^2 / 10
#' D_mat <- diag(rep(0.05, 4))
#' library("ODEnetwork")
#' lfonet <- ODEnetwork(masses = M_mat, dampers = D_mat, springs = K_mat)
#' # The parameters to be included in the sensitivity analysis and their lower
#' # and upper boundaries:
#' LFOpars <- c("k.1", "k.2", "k.3", "k.4",
#'              "d.1", "d.2", "d.3", "d.4")
#' LFObinf <- c(rep(0.2, 4), rep(0.01, 4))
#' LFObsup <- c(rep(20, 4), rep(0.1, 4))
#' # Setting of the initial values of the state variables:
#' lfonet <- setState(lfonet, state1 = rep(2, 4), state2 = rep(0, 4))
#' # The timepoints of interest:
#' LFOtimes <- seq(25, 150, by = 2.5)
#' # Sobol' sensitivity analysis (here only with n = 500 since n = 1000 takes
#' # a lot more time; warnings are suppressed):
#' set.seed(1739)
#' suppressWarnings(
#'   LFOres_sobol <- ODEsobol(mod = lfonet,
#'                            pars = LFOpars,
#'                            times = LFOtimes,
#'                            n = 500,
#'                            rfuncs = "runif",
#'                            rargs = paste0("min = ", LFObinf,
#'                                           ", max = ", LFObsup),
#'                            sobol_method = "Martinez",
#'                            parallel_eval = TRUE,
#'                            parallel_eval_ncores = 2)
#' )
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
                                n = 1000,
                                rfuncs = "runif",
                                rargs = "min = 0, max = 1",
                                sobol_method = "Martinez",
                                ode_method = "lsoda",
                                parallel_eval = FALSE,
                                parallel_eval_ncores = NA, ...){
  
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
    if(length(rfuncs) == length(pars)){
      rfuncs <- rfuncs[!duplicated(pars)]
    }
    if(length(rfuncs) == length(pars)){
      rargs <- rargs[!duplicated(pars)]
    }
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
  assertIntegerish(n, lower = 2)
  assertCharacter(rfuncs)
  if(! length(rfuncs) %in% c(1, length(pars))){
    stop("Argument \"rfuncs\" must be of length 1 or of the same length as ",
         "\"pars\"")
  }
  assertCharacter(rargs)
  if(! length(rargs) %in% c(1, length(pars))){
    stop("Argument \"rargs\" must be of length 1 or of the same length as ",
         "\"pars\"")
  }
  if(! all(sapply(rfuncs, exists))){
    stop("At least one of the supplied functions in \"rfuncs\" was not found")
  }
  stopifnot(sobol_method %in% c("Jansen", "Martinez"))
  stopifnot(ode_method %in% c("lsoda", "lsode", "lsodes","lsodar","vode", 
                              "daspk", "euler", "rk4", "ode23", "ode45", 
                              "radau", "bdf", "bdf_d", "adams", "impAdams", 
                              "impAdams_d" ,"iteration"))
  assertLogical(parallel_eval, len = 1)
  assertIntegerish(parallel_eval_ncores, len = 1, lower = 1)
  if(parallel_eval && is.na(parallel_eval_ncores)){
    parallel_eval_ncores <- 1
  }
  
  ##### Preparation ####################################################
  
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
      return(simnet_res$simulation$results[2:(timesNum + 1), 
                                           2:(z + 1), drop = FALSE])
    }
    if(parallel_eval){
      # Run one_par() on parallel nodes:
      local_cluster <- parallel::makePSOCKcluster(names = parallel_eval_ncores)
      parallel::clusterExport(local_cluster, 
                              varlist = c("pars", "mod", "z", "X",
                                          "times", "timesNum", "ode_method"),
                              envir = environment())
      res_per_par <- parallel::parSapply(local_cluster, 1:nrow(X), one_par, 
                                         simplify = "array")
      parallel::stopCluster(local_cluster)
    } else{
      # Just use sapply() with "simplify = "array"":
      res_per_par <- sapply(1:nrow(X), one_par, simplify = "array")
    }
    res_per_state <- aperm(res_per_par, perm = c(3, 1, 2))
    dimnames(res_per_state) <- list(NULL, paste0("time", 1:timesNum), 
                                    names(state_init))
    return(res_per_state)
  }
  
  # Create the two matrices containing the parameter samples for Monte Carlo
  # estimation:
  if(length(rfuncs) == 1){
    rfuncs <- rep(rfuncs, length(pars))
  }
  if(length(rargs) == 1){
    rargs <- rep(rargs, length(pars))
  }
  rfunc_calls <- paste0(rfuncs, "(n, ", rargs, ")", collapse = ", ")
  X1 <- matrix(eval(parse(text = paste0("c(", rfunc_calls, ")"))), ncol = k)
  X2 <- matrix(eval(parse(text = paste0("c(", rfunc_calls, ")"))), ncol = k)
  colnames(X1) <- colnames(X2) <- pars
  
  ##### Sensitivity analysis #########################################
  
  # Sensitivity analysis with either soboljansen() or sobolmartinez()
  # from package "sensitivity":
  if(sobol_method == "Jansen"){
    x <- soboljansen(model = model_fit, X1, X2, nboot = 0)
  } else if(sobol_method == "Martinez"){
    x <- sobolmartinez(model = model_fit, X1, X2, nboot = 0)
  }
  
  # Process the results:
  ST_by_state <- sobol_process(x, pars, times)
  
  # Return:
  class(ST_by_state) <- "ODEsobol"
  attr(ST_by_state, "sobol_method") <- sobol_method
  return(ST_by_state)
}
