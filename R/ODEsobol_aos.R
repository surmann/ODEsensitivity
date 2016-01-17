#' @title Sobol' SA for ODEs for All Output Variables and All Timepoints
#' Simultaneously
#'
#' @description
#' \code{ODEsobol_aos} performs a variance-based sensitivity analysis for
#' ordinary differential equations according to either the Sobol'-Jansen- or the
#' Sobol'-Martinez-method. The analysis is done for all output variables at all 
#' timepoints simultaneously using \code{\link[sensitivity]{soboljansen_list}} 
#' or \code{\link[sensitivity]{sobolmartinez_list}} from the package 
#' \code{sensitivity}.
#'
#' @param mod [\code{function(Time, State, Pars)}]\cr
#'   model to examine, cf. example below.
#' @param pars [\code{character(k)}]\cr
#'   vector of \code{k} input variable names.
#' @param yini [\code{numeric(z)}]\cr
#'   vector of \code{z} initial values. Must be named (with unique names).
#' @param times [\code{numeric}]\cr
#'   points of time at which the SA should be executed
#'   (vector of arbitrary length). Also the
#'   first point of time must be positive.
#' @param ode_method [\code{character(1)}]\cr
#'   method to be used for solving the differential equations, see 
#'   \code{\link[deSolve]{ode}}. Defaults to \code{"lsoda"}.
#' @param ode_parallel [\code{logical(1)}]\cr
#'   logical indicating if a parallelization shall be done for computing the
#'   \code{\link[deSolve]{ode}}-results for the different parameter combinations
#'   generated for Monte Carlo estimation of the SA indices.
#' @param ode_parallel_ncores [\code{integer(1)}]\cr
#'   number of processor cores to be used for parallelization. Only applies if
#'   \code{ode_parallel = TRUE}. Default is 1.
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
#' @param method [\code{character(1)}]\cr
#'   either \code{"jansen"} or \code{"martinez"}, specifying which modification
#'   of the variance-based Sobol' method shall be used. Defaults to 
#'   \code{"martinez"}, which is slightly faster than \code{"jansen"}.
#' @param nboot [\code{integer(1)}]\cr
#'   parameter \code{nboot} used in \code{\link{soboljansen_list}} resp.
#'   \code{\link{sobolmartinez_list}}, i.e. the number of bootstrap 
#'   replicates. Defaults to 0, so no bootstrapping is done.
#'
#' @return List of length \code{length(yini)} and of class \code{sobolRes_aos} 
#'   containing in each element a list of the Sobol' SA results for the 
#'   corresponding \code{yini}-variable (i.e. 1st order sensitivity indices
#'   \code{S} and total sensitivity indices \code{T}) for every point of
#'   time of the \code{times} vector.
#'
#' @details \code{ODEsobol_aos} uses 
#' \code{\link[sensitivity]{soboljansen_list}} resp.
#' \code{\link[sensitivity]{sobolmartinez_list}} which can handle lists 
#' as output for their model functions. Thus, each element of the list can be 
#' used to contain the results for one output variable. This saves time since 
#' \code{\link[deSolve]{ode}} from the package \code{deSolve} does its 
#' calculations for all output variables anyway, so \code{\link[deSolve]{ode}} 
#' only needs to be executed once.
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
#'   indices might be very bad and even produce negative sensitivity indices (up
#'   to now, this problem only occured for first order indices). Sensitivity 
#'   indices in the interval [-0.05, 0) are considered as minor deviations and 
#'   set to 0 without a warning. Sensitivity indices lower than -0.05 are 
#'   considered as major deviations. They remain unchanged and a warning is 
#'   thrown.
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
#' ##### FitzHugh-Nagumo equations (Ramsay et al., 2007)
#' # definition of the model itself, parameters, initial values
#' # and the times vector:
#'
#' @import checkmate
#' @importFrom deSolve ode
#' @importFrom sensitivity soboljansen_list
#' @importFrom sensitivity sobolmartinez_list
#' @export
#'

ODEsobol_aos <- function(mod,
                         pars,
                         yini,
                         times,
                         ode_method = "lsoda",
                         ode_parallel = FALSE,
                         ode_parallel_ncores = 1,
                         seed = 2015,
                         n = 1000,
                         rfuncs = rep("runif", length(pars)),
                         rargs = rep("min = 0, max = 1", length(pars)),
                         method = "martinez",
                         nboot = 0) {

  ##### Input-Checks #################################################
  assertFunction(mod)
  assertCharacter(pars)
  assertNumeric(yini)
  assertNamed(yini, type = "unique")
  assertNumeric(times, lower = 0, finite = TRUE, unique = TRUE)
  times <- sort(times)
  stopifnot(!any(times == 0))
  stopifnot(ode_method %in% c("lsoda", "lsode", "lsodes","lsodar","vode", 
                              "daspk", "euler", "rk4", "ode23", "ode45", 
                              "radau", "bdf", "bdf_d", "adams", "impAdams", 
                              "impAdams_d" ,"iteration"))
  assertLogical(ode_parallel, len = 1)
  assertIntegerish(ode_parallel_ncores, len = 1, lower = 1)
  assertNumeric(seed)
  assertIntegerish(n)
  assertCharacter(rfuncs, len = length(pars))
  assertCharacter(rargs, len = length(pars))
  rfuncs_exist <- sapply(rfuncs, exists)
  if(!all(rfuncs_exist)) stop(paste("At least one of the supplied functions",
                                    "in \"rfuncs\" was not found"))
  stopifnot(method %in% c("jansen", "martinez"))
  assertIntegerish(nboot)

  ##### Vorarbeiten ####################################################
  set.seed(seed)
  # Anzahl Parameter:
  k <- length(pars)
  # Anzahl Outputgroessen:
  z <- length(yini)
  # Anzahl Zeitpunkte von Interesse:
  timesNum <- length(times)
  # Forme DGL-Modell um, sodass fuer soboljansen_list()- bzw.
  # sobolmartinez_list()-Argument "model" passend:
  model_fit <- function(X){
    # X   - (nxk)-Matrix mit den n einzugebenden Parameter-Konstellationen
    #       als Zeilen
    colnames(X) <- pars
    one_par <- function(i){
      ode(yini, times = c(0, times), mod, parms = X[i, ], 
          method = ode_method)[2:(timesNum + 1), 2:(z + 1)]
    }
    if(ode_parallel){
      ode_cl <- parallel::makeCluster(rep("localhost", ode_parallel_ncores), 
                                      type = "SOCK")
      parallel::clusterExport(ode_cl, varlist = c("ode", "mod", "yini", "z", 
                                                  "X", "times", "timesNum"),
                              envir = environment())
      res_per_par <- parallel::parLapply(ode_cl, 1:nrow(X), one_par)
      parallel::stopCluster(ode_cl)
    } else{
      res_per_par <- lapply(1:nrow(X), one_par)
    }
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
  
  ##### Sensitivitaet ##################################################
  rfunc_calls <- paste0(rfuncs, "(n, ", rargs, ")", collapse = ", ")
  X1 <- matrix(eval(parse(text = paste0("c(", rfunc_calls, ")"))), ncol = k)
  X2 <- matrix(eval(parse(text = paste0("c(", rfunc_calls, ")"))), ncol = k)
  colnames(X1) <- colnames(X2) <- pars
  
  # Listen der Sensitivitaetsindizes (Haupteffekt, total) zu den
  # interessierenden Zeitpunkten:
  # S <- T <- matrix(nrow = 1 + k, ncol = timesNum)
  
  # Durchfuehrung der Sensitivitaetsanalyse mit den Funktionen aus dem Paket
  # "sensitivity":
  if(method == "jansen"){
    x <- soboljansen_list(model = model_fit, X1, X2, nboot = nboot)
  } else if(method == "martinez"){
    x <- sobolmartinez_list(model = model_fit, X1, X2, nboot = nboot)
  }
  
  # Verarbeitung der Ergebnisse:
  ST_original_by_y <- lapply(x$ST_by_y, function(L){
    ST_original <- sapply(L, function(ST_col){
      c(ST_col$S[, 1], ST_col$T[, 1])
    })
    return(ST_original)
  })
  # ST_original wieder in 2 Matrizen S und T aufspalten:
  ST_by_y <- lapply(ST_original_by_y, function(ST_original){
    S <- rbind(times, ST_original[1:k, , drop = FALSE])
    T <- rbind(times, ST_original[(k+1):(2*k), , drop = FALSE])
    rownames(S) <- rownames(T) <- c("time", pars)
    return(list(S = S, T = T))
  })
  
  # Handling of negative first order SA indices ("minor": >= -0.05 and < 0, 
  # "major": < -0.05):
  check_negative <- sapply(ST_by_y, function(ST){
    c(minor = any(-0.05 <= ST$S & ST$S < 0), major = any(ST$S < -0.05))
  })
  if(any(check_negative["minor", ])){
    # "Repair" minor negative SA indices by setting them to zero:
    for(i in seq_along(ST_by_y)[check_negative[1, ]]){
      check_minor <- -0.05 <= ST_by_y[[i]]$S & ST_by_y[[i]]$S < 0
      ST_by_y[[i]]$S[check_minor] <- 0
    }
  }
  if(any(check_negative["major", ])){
    warning("Negative sensitivity indices (< -0.05) detected. Argument \"n\" ",
            "might be too low. If using a higher value for \"n\" does not ",
            "help, please check the parameter distributions (\"rfuncs\") and ",
            "their arguments (\"rargs\").")
  }
  
  # Rueckgabe:
  res <- list(ST_by_y = ST_by_y, method = method)
  class(res) <- "sobolRes_aos"
  return(res)
}
