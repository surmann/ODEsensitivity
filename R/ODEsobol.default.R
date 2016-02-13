#' @title Sobol' SA for general ODE models
#'
#' @description
#' \code{ODEsobol.default} performs a variance-based sensitivity analysis for
#' ordinary differential equations according to either the Sobol'-Jansen- or the
#' Sobol'-Martinez-method. The analysis is done for all output variables at all 
#' timepoints simultaneously using \code{\link[sensitivity]{soboljansen}} 
#' or \code{\link[sensitivity]{sobolmartinez}} from the package 
#' \code{sensitivity}.
#'
#' @param mod [\code{function(Time, State, Pars)}]\cr
#'   model to examine, cf. example below.
#' @param pars [\code{character(k)}]\cr
#'   vector of \code{k} input variable names.
#' @param state_init [\code{numeric(z)}]\cr
#'   vector of \code{z} initial values. Must be named (with unique names).
#' @param times [\code{numeric}]\cr
#'   points of time at which the SA should be executed
#'   (vector of arbitrary length). Also the
#'   first point of time must be positive.
#' @param seed [\code{numeric(1)}]\cr
#'   seed.
#' @param n [\code{integer(1)}]\cr
#'   number of random parameter values (\code{n} per input factor) used to 
#'   estimate the variance-based sensitivity indices by Monte Carlo method.
#'   (Variance-based methods for sensitivity analysis rely on 
#'   Monte-Carlo-simulation to estimate the integrals needed for the calculation
#'   of the sensitivity indices.) Defaults to 1000.
#' @param rfuncs [\code{character(k)}]\cr
#'   names of the \code{k} functions used to generate the \code{n} random values
#'   for the \code{k} parameters. This way, different distributions can be 
#'   used for the \code{k} parameters. Defaults to \code{"runif"} for each of
#'   the \code{k} parameters.
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
#' @param nboot [\code{integer(1)}]\cr
#'   parameter \code{nboot} used in \code{\link{soboljansen}} resp.
#'   \code{\link{sobolmartinez}}, i.e. the number of bootstrap 
#'   replicates. Defaults to 0, so no bootstrapping is done.
#' @param ode_method [\code{character(1)}]\cr
#'   method to be used for solving the differential equations, see 
#'   \code{\link[deSolve]{ode}}. Defaults to \code{"lsoda"}.
#' @param ode_parallel [\code{logical(1)}]\cr
#'   logical indicating if a parallelization shall be done for computing the
#'   \code{\link[deSolve]{ode}}-results for the different parameter combinations
#'   generated for Monte Carlo estimation of the SA indices.
#' @param ode_parallel_ncores [\code{integer(1)}]\cr
#'   number of processor cores to be used for parallelization. Only applies if
#'   \code{ode_parallel = TRUE}. If set to \code{NA} (as per default) and 
#'   \code{ode_parallel = TRUE}, 1 processor core is used.
#' @param ... further arguments passed to or from other methods.
#'
#' @return 
#'   List of length \code{length(state_init)} and of class \code{sobolRes} 
#'   containing in each element a list of the Sobol' SA results for the 
#'   corresponding \code{state_init}-variable (i.e. first order sensitivity 
#'   indices \code{S} and total sensitivity indices \code{T}) for every point of
#'   time of the \code{times} vector.
#'
#' @details 
#'   \code{ODEsobol} uses \code{\link[sensitivity]{soboljansen}} resp.
#'   \code{\link[sensitivity]{sobolmartinez}} which can handle lists 
#'   as output for their model functions. Thus, each element of the list can be 
#'   used to contain the results for one output variable. This saves time since 
#'   \code{\link[deSolve]{ode}} from the package \code{deSolve} does its 
#'   calculations for all output variables anyway, so \code{\link[deSolve]{ode}} 
#'   only needs to be executed once.
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
#' @references J. O. Ramsay, G. Hooker, D. Campbell and J. Cao, 2007,
#'   \emph{Parameter estimation for differential equations: a generalized 
#'   smoothing approach}, Journal of the Royal Statistical Society, Series B, 
#'   69, Part 5, 741--796.
#' @seealso \code{\link[sensitivity]{soboljansen},
#'   \link[sensitivity]{sobolmartinez},
#'   \link{plot.sobolRes}}
#' 
#' @examples
#' ##### FitzHugh-Nagumo equations (Ramsay et al., 2007) #####
#' # Definition of the model itself, parameters, initial state values
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
#' FHNstate  <- c(Voltage = -1, Current = 1)
#' FHNtimes <- seq(0.1, 50, by = 5)
#' 
#' ### The following code might take a long time!
#' FHNres <- ODEsobol(mod = FHNmod,
#'                    pars = c("a", "b", "s"),
#'                    state_init = FHNstate,
#'                    times = FHNtimes,
#'                    seed = 2015,
#'                    n = 1000,
#'                    rfuncs = c("runif", "runif", "rnorm"),
#'                    rargs = c(rep("min = 0.18, max = 0.22", 2),
#'                              "mean = 3, sd = 0.2 / 3"),
#'                    sobol_method = "martinez",
#'                    nboot = 0,
#'                    ode_method = "adams",
#'                    ode_parallel = TRUE,
#'                    ode_parallel_ncores = 2)
#'
#' @import checkmate
#' @importFrom deSolve ode
#' @importFrom sensitivity soboljansen
#' @importFrom sensitivity sobolmartinez
#' @export
#'

ODEsobol.default <- function(mod,
                             pars,
                             state_init,
                             times,
                             seed = 2015,
                             n = 1000,
                             rfuncs = rep("runif", length(pars)),
                             rargs = rep("min = 0, max = 1", length(pars)),
                             sobol_method = "martinez",
                             nboot = 0,
                             ode_method = "lsoda",
                             ode_parallel = FALSE,
                             ode_parallel_ncores = NA, ...){

  ##### Input checks ###################################################
  
  assertFunction(mod)
  assertCharacter(pars)
  assertNumeric(state_init)
  assertNamed(state_init, type = "unique")
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
  assertIntegerish(nboot)
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
  # Number of parameters:
  k <- length(pars)
  # Number of state variables:
  z <- length(state_init)
  # Number of timepoints:
  timesNum <- length(times)
  
  # Adapt the ODE model for argument "model" of soboljansen() resp.
  # sobolmartinez():
  model_fit <- function(X){
    # Input: Matrix X with k columns, containing the random parameter 
    # combinations.
    colnames(X) <- pars
    one_par <- function(i){
      # Resolve the ODE system by using ode() from the package "deSolve":
      ode(state_init, times = c(0, times), mod, parms = X[i, ], 
          method = ode_method)[2:(timesNum + 1), 2:(z + 1)]
    }
    if(ode_parallel){
      # Run one_par() on parallel nodes:
      ode_cl <- parallel::makeCluster(rep("localhost", ode_parallel_ncores), 
                                      type = "SOCK")
      parallel::clusterExport(ode_cl, 
                              varlist = c("ode", "mod", "state_init", "z", "X",
                                          "times", "timesNum", "ode_method"),
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
  
  ##### Sensitivity analysis ###########################################
  
  # Sensitivity analysis with either soboljansen() or sobolmartinez()
  # from package "sensitivity":
  if(sobol_method == "jansen"){
    x <- soboljansen(model = model_fit, X1, X2, nboot = nboot)
  } else if(sobol_method == "martinez"){
    x <- sobolmartinez(model = model_fit, X1, X2, nboot = nboot)
  }
  
  # Process the results:
  ST_by_state <- sobol_process(x, pars, times)
  
  # Return:
  res <- list(ST_by_state = ST_by_state, sobol_method = sobol_method)
  class(res) <- "sobolRes"
  return(res)
}
