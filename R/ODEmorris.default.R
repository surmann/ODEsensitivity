#' @title Morris Screening for General ODE Models
#'
#' @description
#' \code{ODEmorris.default} is the default method of \code{\link{ODEmorris}}. It
#' performs a sensitivity analysis for general ODE models using the Morris 
#' screening method.
#'
#' @param mod [\code{function(Time, State, Pars)}]\cr
#'   model to examine, supplied in the manner as needed for 
#'   \code{\link[deSolve]{ode}} (see example below).
#' @param pars [\code{character(k)}]\cr
#'   names of the parameters to be included as input variables in Morris 
#'   screening.
#' @param state_init [\code{numeric(z)}]\cr
#'   vector of \code{z} initial values. Must be named (with unique names).
#' @param times [\code{numeric}]\cr
#'   points of time at which the sensitivity analysis should be executed (vector
#'   of arbitrary length). The first point of time must be greater than zero.
#' @param binf [\code{character(1} or \code{k)}]\cr
#'   vector of lower borders of possible input parameter values.
#'   If they are all equal, a single value can be set.
#' @param bsup [\code{character(1} or \code{k)}]\cr
#'   vector of upper borders of possible input parameter values.
#'   If they are all equal, a single value can be set.
#' @param r [\code{integer(1 or 2)}]\cr
#'   if of length 1, the number of repetitions of the \code{design}. If of 
#'   length 2, a space-filling optimization of the sampling design is used, see 
#'   \code{\link[sensitivity]{morris}}. However, this space-filling optimization
#'   might lead to long runtimes, so length 1 is recommended for \code{r}. 
#'   Defaults to 500.
#' @param design [\code{list}]\cr
#'   a list specifying the design type and its parameters,
#'   cf. \code{\link[sensitivity]{morris}}.
#' @param scale [\code{logical(1)}]\cr
#'   if \code{TRUE}, scaling is done for the input design of experiments after 
#'   building the design and before calculating the elementary effects,
#'   cf. \code{\link[sensitivity]{morris}}. Defaults to \code{TRUE}, which is
#'   highly recommended if the factors have different orders of magnitude, see
#'   \code{\link[sensitivity]{morris}}.
#' @param ode_method [\code{character(1)}]\cr
#'   method to be used for solving the differential equations, see 
#'   \code{\link[deSolve]{ode}}. Defaults to \code{"lsoda"}.
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
#'   List of class \code{ODEmorris} of length \code{length(state_init)} 
#'   containing in each element a matrix for one state variable. The
#'   matrices itself contain the Morris screening results for all timepoints 
#'   (rows: \code{mu, mu.star} and \code{sigma} for every parameter; columns: 
#'   timepoints).
#'
#' @details
#'   Function \code{\link[deSolve]{ode}} from \code{\link[deSolve]{deSolve}} is 
#'   used to solve the ODE system.
#'   
#'   The sensitivity analysis is done for all state variables and all
#'   timepoints simultaneously using \code{\link[sensitivity]{morris}} from the 
#'   package \code{sensitivity}.
#'   
#'   For non-ODE models, values for \code{r} are typically between 10 and 50.
#'   However, much higher values are recommended for ODE models (the default is
#'   \code{r = 500}).
#' 
#' @note 
#'   If the evaluation of the model function takes too long, it might be helpful 
#'   to try another ODE-solver (argument \code{ode_method}). The 
#'   \code{ode_method}s \code{"vode"}, \code{"bdf"}, \code{"bdf_d"}, 
#'   \code{"adams"}, \code{"impAdams"} and \code{"impAdams_d"} might be faster 
#'   than the default \code{"lsoda"}.
#'   
#'   If \code{\link[sensitivity]{morris}} throws a warning message stating
#'   "In ... keeping ... repetitions out of ...", try using a bigger number of 
#'   \code{levels} in the \code{design} argument (only possible for OAT design).
#'
#' @author Stefan Theers, Frank Weber
#' @references J. O. Ramsay, G. Hooker, D. Campbell and J. Cao, 2007,
#'   \emph{Parameter estimation for differential equations: a generalized 
#'   smoothing approach}, Journal of the Royal Statistical Society, Series B, 
#'   69, Part 5, 741--796.
#' @seealso \code{\link[sensitivity]{morris}, \link{plot.ODEmorris}}
#' 
#' @examples
#' ##### Lotka-Volterra equations #####
#' # The model function:
#' LVmod <- function(Time, State, Pars) {
#'   with(as.list(c(State, Pars)), {
#'     Ingestion    <- rIng  * Prey * Predator
#'     GrowthPrey   <- rGrow * Prey * (1 - Prey/K)
#'     MortPredator <- rMort * Predator
#'     
#'     dPrey        <- GrowthPrey - Ingestion
#'     dPredator    <- Ingestion * assEff - MortPredator
#'     
#'     return(list(c(dPrey, dPredator)))
#'   })
#' }
#' # The parameters to be included in the sensitivity analysis and their lower 
#' # and upper boundaries:
#' LVpars  <- c("rIng", "rGrow", "rMort", "assEff", "K")
#' LVbinf <- c(0.05, 0.05, 0.05, 0.05, 1)
#' LVbsup <- c(1.00, 3.00, 0.95, 0.95, 20)
#' # The initial values of the state variables:
#' LVinit  <- c(Prey = 1, Predator = 2)
#' # The timepoints of interest:
#' LVtimes <- c(0.01, seq(1, 50, by = 1))
#' # Morris screening:
#' set.seed(7292)
#' \dontrun{
#' LVres_morris <- ODEmorris(mod = LVmod,
#'                           pars = LVpars,
#'                           state_init = LVinit,
#'                           times = LVtimes,
#'                           binf = LVbinf,
#'                           bsup = LVbsup,
#'                           r = 500,
#'                           design = list(type = "oat", 
#'                                         levels = 10, grid.jump = 1),
#'                           scale = TRUE,
#'                           ode_method = "lsoda",
#'                           parallel_eval = TRUE,
#'                           parallel_eval_ncores = 2)
#' }
#' 
#' ##### FitzHugh-Nagumo equations (Ramsay et al., 2007) #####
#' FHNmod <- function(Time, State, Pars) {
#'   with(as.list(c(State, Pars)), {
#'     
#'     dVoltage <- s * (Voltage - Voltage^3 / 3 + Current)
#'     dCurrent <- - 1 / s *(Voltage - a + b * Current)
#'     
#'     return(list(c(dVoltage, dCurrent)))
#'   })
#' }
#' \dontrun{
#' FHNres_morris <- ODEmorris(mod = FHNmod,
#'                            pars = c("a", "b", "s"),
#'                            state_init = c(Voltage = -1, Current = 1),
#'                            times = seq(0.1, 50, by = 5),
#'                            binf = c(0.18, 0.18, 2.8),
#'                            bsup = c(0.22, 0.22, 3.2),
#'                            r = 500,
#'                            design = list(type = "oat", 
#'                                          levels = 50, grid.jump = 1),
#'                            scale = TRUE,
#'                            ode_method = "adams",
#'                            parallel_eval = TRUE,
#'                            parallel_eval_ncores = 2)
#' }
#'
#' @import checkmate
#' @importFrom deSolve ode
#' @importFrom sensitivity morris
#' @export
#'

ODEmorris.default <- function(mod,
                              pars,
                              state_init,
                              times,
                              binf = 0,
                              bsup = 1,
                              r = 500,
                              design =
                                list(type = "oat", levels = 10, grid.jump = 1),
                              scale = TRUE,
                              ode_method = "lsoda",
                              parallel_eval = FALSE,
                              parallel_eval_ncores = NA, ...){
  
  ##### Input checks ###################################################
  
  assertFunction(mod)
  assertCharacter(pars)
  assertNumeric(state_init)
  assertNamed(state_init, type = "unique")
  assertNumeric(times, lower = 0, finite = TRUE, unique = TRUE)
  times <- sort(times)
  stopifnot(!any(times == 0))
  assertNumeric(binf)
  if(length(binf) != length(pars) && length(binf) != 1)
    stop("\"binf\" must be of length 1 or of the same length as \"pars\".")
  assertNumeric(bsup)
  if(length(bsup) != length(pars) & length(bsup) != 1)
    stop("\"bsup\" must be of length 1 or of the same length as \"pars\".")
  if(any(bsup < binf)){
    idx_swap <- bsup < binf
    bsup_tmp <- bsup[idx_swap]
    bsup[idx_swap] <- binf[idx_swap]
    binf[idx_swap] <- bsup_tmp
    warning("At least one element of \"bsup\" was lower than the ", 
            "corresponding element of \"binf\". Elements were swapped.")
  }
  assertIntegerish(r, lower = 1, min.len = 1, max.len = 2)
  assertList(design)
  assertLogical(scale, len = 1)
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
  
  # Number of parameters:
  k <- length(pars)
  # Number of state variables:
  z <- length(state_init)
  # Number of timepoints:
  timesNum <- length(times)
  
  # Adapt the ODE model for argument "model" of morris():
  model_fit <- function(X){
    # Input: Matrix X with k columns, containing the random parameter 
    # combinations.
    colnames(X) <- pars
    one_par <- function(i){
      # Resolve the ODE system by using ode() from the package "deSolve":
      ode(state_init, times = c(0, times), mod, parms = X[i, ], 
          method = ode_method)[2:(timesNum + 1), 2:(z + 1), drop = FALSE]
    }
    if(parallel_eval){
      # Run one_par() on parallel nodes:
      local_cluster <- parallel::makePSOCKcluster(names = parallel_eval_ncores)
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
  
  ##### Sensitivity analysis ###########################################
  
  # Sensitivity analysis with function morris() from package "sensitivity":
  x <- morris(model = model_fit, factors = pars, r = r, design = design, 
              binf = binf, bsup = bsup, scale = scale)
  
  # Process the results:
  mu <- lapply(1:dim(x$ee)[4], function(i){
    apply(x$ee[, , , i, drop = FALSE], 3, function(M){
      apply(M, 2, mean)
    })
  })
  mu.star <- lapply(1:dim(x$ee)[4], function(i){
    apply(abs(x$ee)[, , , i, drop = FALSE], 3, function(M){
      apply(M, 2, mean)
    })
  })
  sigma <- lapply(1:dim(x$ee)[4], function(i){
    apply(x$ee[, , , i, drop = FALSE], 3, function(M){
      apply(M, 2, sd)
    })
  })
  names(mu) <- names(mu.star) <- names(sigma) <- dimnames(x$ee)[[4]]
  
  out_all_states <- lapply(1:length(mu), function(i){
    one_state <- rbind(times, mu[[i]], mu.star[[i]], sigma[[i]])
    rownames(one_state) <- c("time", paste0("mu_", pars), 
                             paste0("mu.star_", pars),
                             paste0("sigma_", pars))
    return(one_state)
  })
  names(out_all_states) <- names(mu)
  
  # Throw a warning if NAs occur (probably there are parameter combinations
  # which are not suitable, so the ODE system can't be solved):
  NA_check_mu <- function(M){
    any(is.na(M[1:(1 + k*2), ]))
  }
  NA_check_sigma <- function(M){
    all(is.na(M[(2 + k*2):(1 + k*3), ]))
  }
  if(any(unlist(lapply(out_all_states, NA_check_mu)))){
    warning(paste("The ODE system can't be solved. This might be due to", 
      "arising unrealistic parameters by means of Morris's Screening. Use",
      "ODEsobol() instead or set binf and bsup differently together with",
      "scale = TRUE. It might also be helpful to try another ODE-solver by",
      "using the \"ode_method\"-argument."))
  } else if(all(unlist(lapply(out_all_states, NA_check_sigma))) && r[1] == 1){
    warning("Calculation of sigma requires r >= 2.")
  } else{
    NA_check_sigma_any <- function(M){
      any(is.na(M[(2 + k*2):(1 + k*3), ]))
    }
    if(any(unlist(lapply(out_all_states, NA_check_sigma_any)))){
      warning("NAs for sigma. This might be due to r being too small.")
    }
  }
  
  # Return:
  class(out_all_states) <- "ODEmorris"
  return(out_all_states)
}
