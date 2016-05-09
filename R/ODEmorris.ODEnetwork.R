#' @title Morris Screening for Objects of Class \code{ODEnetwork}
#'
#' @description
#' \code{ODEmorris.ODEnetwork} performs a sensitivity analysis for objects of 
#' class \code{ODEnetwork} using the Morris screening method.
#' Package \code{ODEnetwork} is required for this function to work.
#'
#' @param mod [\code{ODEnetwork}]\cr
#'   list of class \code{ODEnetwork}.
#' @param pars [\code{character(k)}]\cr
#'   names of the parameters to be included as input variables in Morris 
#'   screening. All parameter names must be contained in 
#'   \code{names(ODEnetwork::createParamVec(mod))} and must not be derivable 
#'   from other parameters supplied (e.g. \code{"k.2.1"} can be derived from 
#'   \code{"k.1.2"}, so supplying \code{"k.1.2"} suffices).
#' @param times [\code{numeric}]\cr
#'   points of time at which the sensitivity analysis should be executed (vector
#'   of arbitrary length). The first point of time must be greater than zero.
#' @param binf [\code{character(1} or \code{k)}]\cr
#'   vector of lower borders of possible values for the \code{k} input 
#'   parameters. If they are all equal, a single value can be set.
#' @param bsup [\code{character(1} or \code{k)}]\cr
#'   vector of upper borders of possible values for the \code{k} input 
#'   parameters. If they are all equal, a single value can be set.
#' @param r [\code{integer(1)}]\cr
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
#'   List of class \code{ODEmorris} of length \code{2 * nrow(mod$state)} 
#'   containing in each element a matrix for one state variable (all components 
#'   of the 2 state variables are analyzed independently). The matrices 
#'   themselves contain the Morris screening results for all timepoints (rows: 
#'   \code{mu, mu.star} and \code{sigma} for every parameter; columns: 
#'   timepoints).
#'
#' @details
#'   If the object of class \code{ODEnetwork} supplied for \code{mod} doesn't
#'   include any events, the solution of the ODE network is determined 
#'   analytically using \code{\link[ODEnetwork]{simuNetwork}}. In the presence
#'   of events, \code{\link[ODEnetwork]{simuNetwork}} uses 
#'   \code{\link[deSolve]{ode}} to solve the ODE network numerically.
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
#'   In situations where the solution of the ODE model has to be determined 
#'   numerically, it might be helpful to try another ODE-solver if the 
#'   evaluation of the model function takes too long, (argument 
#'   \code{ode_method}). The \code{ode_method}s \code{"vode"}, \code{"bdf"}, 
#'   \code{"bdf_d"}, \code{"adams"}, \code{"impAdams"} and \code{"impAdams_d"} 
#'   might be faster than the default \code{"lsoda"}.
#'   
#'   If \code{\link[sensitivity]{morris}} throws a warning message stating
#'   "In ... keeping ... repetitions out of ...", try using a bigger number of 
#'   \code{levels} in the \code{design} argument (only possible for OAT design).
#'
#' @author Frank Weber
#' @seealso \code{\link[sensitivity]{morris}, \link{plot.ODEmorris}}
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
#' # Morris screening:
#' set.seed(283)
#' \dontrun{
#' LFOres_morris <- ODEmorris(mod = lfonet,
#'                            pars = LFOpars,
#'                            times = LFOtimes,
#'                            binf = LFObinf,
#'                            bsup = LFObsup,
#'                            r = 500,
#'                            design = list(type = "oat", 
#'                                          levels = 10, grid.jump = 1),
#'                            scale = TRUE,
#'                            parallel_eval = TRUE,
#'                            parallel_eval_ncores = 2)
#' }
#'
#' @import checkmate
#' @importFrom deSolve ode
#' @importFrom sensitivity morris
#' @method ODEmorris ODEnetwork
#' @export
#' 

ODEmorris.ODEnetwork <- function(mod,
                                 pars,
                                 times,
                                 binf = 0,
                                 bsup = 1,
                                 r = 500,
                                 design = list(type = "oat", levels = 10, 
                                               grid.jump = 1),
                                 scale = TRUE, 
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
  
  # Initial state:
  state_init <- ODEnetwork::createState(mod)
  # Number of parameters:
  k <- length(pars)
  # Number of output variables (state variables):
  z <- length(state_init)
  # Number of timepoints:
  timesNum <- length(times)
  
  # Adapt the ODE-model for argument "model" of morris():
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
  
  ##### Sensitivity analysis #########################################
  
  # Sensitivity analysis with function morris() from package "sensitivity":
  x <- morris(model = model_fit, factors = pars, r = r, 
              design = design, binf = binf, bsup = bsup, scale = scale)
  
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