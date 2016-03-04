#' @title Sobol' SA for General ODE Models
#'
#' @description
#' \code{ODEsobol.default} is the default method of \code{\link{ODEmorris}}. It
#' performs a variance-based sensitivity analysis for general ODE models.
#'
#' @param mod [\code{function(Time, State, Pars)}]\cr
#'   model to examine, cf. example below.
#' @param pars [\code{character(k)}]\cr
#'   vector of \code{k} input variable names.
#' @param state_init [\code{numeric(z)}]\cr
#'   vector of \code{z} initial values. Must be named (with unique names).
#' @param times [\code{numeric}]\cr
#'   points of time at which the SA should be executed (vector of arbitrary 
#'   length). The first point of time must be greater than zero.
#' @param seed [\code{numeric(1)}]\cr
#'   seed.
#' @param n [\code{integer(1)}]\cr
#'   number of random parameter values (\code{n} per input factor) used to 
#'   estimate the variance-based sensitivity indices by Monte Carlo method.
#'   (Variance-based methods for sensitivity analysis rely on 
#'   Monte-Carlo-simulation to estimate the integrals needed for the calculation
#'   of the sensitivity indices.) Defaults to 1000.
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
#'   The sensitivity analysis is done for all state variables and all
#'   timepoints simultaneously using either
#'   \code{\link[sensitivity]{soboljansen}} or
#'   \code{\link[sensitivity]{sobolmartinez}} from the package 
#'   \code{sensitivity} (depending on \code{sobol_method}). 
#'   \code{\link[sensitivity]{soboljansen}} and
#'   \code{\link[sensitivity]{sobolmartinez}} can handle three-dimensional
#'   arrays as output for their model functions. Each element of the third 
#'   dimension of the output array is used to contain the results for one 
#'   state variable of the ODE system. Each element of the second dimension of
#'   the output array is used for one timepoint. The use of an array as output
#'   saves time (compared to looping over all state variables and all 
#'   timepoints) since \code{\link[deSolve]{ode}} from the package 
#'   \code{deSolve} does its calculations for all state variables and all
#'   timepoints anyway, so \code{\link[deSolve]{ode}} only needs to be executed 
#'   once.
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
#' # Warning: The following code might take a long time!
#' 
#' # Arguments "rfuncs" and "rargs" being of length 1:
#' FHNres_1 <- ODEsobol(mod = FHNmod,
#'                      pars = c("a", "b", "s"),
#'                      state_init = FHNstate,
#'                      times = FHNtimes,
#'                      seed = 2015,
#'                      n = 1000,
#'                      rfuncs = "runif",
#'                      rargs = "min = 0.1, max = 1.5",
#'                      sobol_method = "Martinez",
#'                      ode_method = "adams",
#'                      ode_parallel = TRUE,
#'                      ode_parallel_ncores = 2)
#' 
#' # Arguments "rfuncs" and "rargs" being of the same length as "pars" (the
#' # distributions and their corresponding arguments are chosen more or less 
#' # arbitrarily and shall only demonstrate the use of different distributions):
#' FHNres <- ODEsobol(mod = FHNmod,
#'                    pars = c("a", "b", "s"),
#'                    state_init = FHNstate,
#'                    times = FHNtimes,
#'                    seed = 2015,
#'                    n = 1000,
#'                    rfuncs = c("runif", "rnorm", "rexp"),
#'                    rargs = c("min = 0.18, max = 0.22", 
#'                              "mean = 0.2, sd = 0.2 / 3",
#'                              "rate = 1 / 3"),
#'                    sobol_method = "Martinez",
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
                             rfuncs = "runif",
                             rargs = "min = 0, max = 1",
                             sobol_method = "Martinez",
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
          method = ode_method)[2:(timesNum + 1), 2:(z + 1), drop = FALSE]
    }
    if(ode_parallel){
      # Run one_par() on parallel nodes:
      local_cluster <- parallel::makePSOCKcluster(names = ode_parallel_ncores)
      parallel::clusterExport(local_cluster, 
                              varlist = c("ode", "mod", "state_init", "z", "X",
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
  
  ##### Sensitivity analysis ###########################################
  
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
  class(ST_by_state) <- "sobolRes"
  attr(ST_by_state, "sobol_method") <- sobol_method
  return(ST_by_state)
}
