#' @title Sobol SA for ODEs at All Timepoints Simultaneously
#'
#' @description
#' \code{ODEsobol_ats} performs a sensitivity analysis for
#' ordinary differential equations using the variance-based Sobol-Jansen method.
#' The analysis is done for one output variable at all timepoints simultaneously
#' using by choice either \code{\link[sensitivity]{soboljansen_matrix}} or 
#' \code{\link[sensitivity]{sobolmartinez_matrix}} from the package 
#' \code{sensitivity}.
#'
#' @param mod [\code{function(Time, State, Pars)}]\cr
#'   model to examine, cf. example below.
#' @param pars [\code{character(k)}]\cr
#'   vector of \code{k} input variable names.
#' @param yini [\code{numeric(z)}]\cr
#'   vector of \code{z} initial values.
#' @param times [\code{numeric}]\cr
#'   points of time at which the SA should be executed
#'   (vector of arbitrary length). Also the
#'   first point of time must be positive.
#' @param ode_method [\code{character(1)}]\cr
#'   method to be used for solving the differential equations, see 
#'   \code{\link[deSolve]{ode}}. Defaults to \code{"lsoda"}.
#' @param y_idx [\code{integer(1)}]\cr
#'   index of the output variable to be analyzed. Defaults to 1, so the first
#'   output variable is used.
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
#'   of the variance-based Sobol method shall be used. Defaults to 
#'   \code{"martinez"}, which is slightly faster than \code{"jansen"}.
#' @param nboot [\code{integer(1)}]\cr
#'   parameter \code{nboot} used in \code{\link{soboljansen_matrix}} resp.
#'   \code{\link{sobolmartinez_matrix}}, i.e. the number of bootstrap 
#'   replicates. Defaults to 0, so no bootstrapping is done.
#'
#' @return List of Sobol SA results (i.e. 1st order sensitivity indices
#'   \code{S} and total sensitivity indices \code{T}) for every point of
#'   time of the \code{times} vector, of class \code{sobolRes_ats}.
#'
#' @details \code{ODEsobol_ats} is faster than \code{ODEsobol} since 
#' \code{\link[sensitivity]{soboljansen_matrix}} and
#' \code{\link[sensitivity]{sobolmartinez_matrix}} can handle matrix output for 
#' their model functions. In \code{ODEsobol_ats}, a model function is
#' created returning a matrix with the \code{\link[deSolve]{ode}}-results for 
#' all timepoints (one per column). Thus, \code{\link[deSolve]{ode}} only needs
#' to be executed once.
#'
#' @note \code{ODEsobol_ats} is only purposed for analysing one output 
#' variable.
#' 
#'   Sometimes, it is also helpful to try another ODE-solver (argument 
#'   \code{ode_method}). Problems are known for the
#'   \code{ode_method}s \code{"euler"}, \code{"rk4"} and \code{"ode45"}. 
#'   In contrast, the \code{ode_method}s \code{"vode"}, \code{"bdf"}, 
#'   \code{"bdf_d"}, \code{"adams"}, \code{"impAdams"} and \code{"impAdams_d"} 
#'   might be even faster than the standard \code{ode_method} \code{"lsoda"}.
#'
#' @author Frank Weber
#' @examples
#' ##### FitzHugh-Nagumo equations (Ramsay et al., 2007)
#' # definition of the model itself, parameters, initial values
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
#'
#' FHNyini  <- c(Voltage = -1, Current = 1)
#' FHNtimes <- seq(0.1, 20, by = 0.5)
#'
#' FHNres <- ODEsobol_ats(mod = FHNmod,
#'                        pars = c("a", "b", "s"),
#'                        yini = FHNyini,
#'                        times = FHNtimes,
#'                        ode_method = "adams",
#'                        y_idx = 1,            # only voltage
#'                        seed = 2015,
#'                        n = 10,               # use n >> 10!
#'                        rfuncs = c("runif", "runif", "rnorm"),
#'                        rargs = c(rep("min = 0.18, max = 0.22", 2),
#'                                  "mean = 3, sd = 0.2 / 3"),
#'                        method = "martinez",
#'                        nboot = 0)
#'
#' @seealso \code{\link[sensitivity]{sobol}},
#'   \code{\link[sensitivity]{soboljansen_matrix}},
#'   \code{\link[sensitivity]{sobolmartinez_matrix}},
#'   \code{\link{plot.sobolRes_ats}}
#'
#' @export
#' @import
#'   checkmate
#' @importFrom deSolve ode
#' @importFrom sensitivity soboljansen_matrix
#' @importFrom sensitivity sobolmartinez_matrix
#'

ODEsobol_ats <- function(mod,
                         pars,
                         yini,
                         times,
                         ode_method = "lsoda",
                         y_idx = 1,
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
  assertNumeric(times, lower = 0, finite = TRUE, unique = TRUE)
  times <- sort(times)
  stopifnot(!any(times == 0))
  stopifnot(ode_method %in% c("lsoda", "lsode", "lsodes","lsodar","vode", 
                              "daspk", "euler", "rk4", "ode23", "ode45", 
                              "radau", "bdf", "bdf_d", "adams", "impAdams", 
                              "impAdams_d" ,"iteration"))
  assertIntegerish(y_idx)
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
  # Forme DGL-Modell um, sodass fuer soboljansen_matrix()- bzw.
  # sobolmartinez_matrix()-Argument "model" passend:
  model_fit <- function(X){
    # X   - (nxk)-Matrix mit den n einzugebenden Parameter-Konstellationen
    #       als Zeilen
    colnames(X) <- pars
    res <- apply(X, 1, function(x){
      ode(yini, times = c(0, times), mod, parms = x, 
          method = ode_method)[2:(timesNum + 1), y_idx + 1]
    })
    # (Jede Zeile in "res" steht fuer einen Zeitpunkt.)
    
    if(timesNum == 1){
      # Korrektur noetig, falls timesNum == 1:
      res <- matrix(res)
    } else{
      # Transponiere die Ergebnis-Matrix, sodass jede Spalte fuer einen
      # Zeitpunkt steht:
      res <- t(res)
    }
    return(res)
  }
  
  ##### Sensitivitaet ##################################################
  rfunc_calls <- paste0(rfuncs, "(n, ", rargs, ")", collapse = ", ")
  X1 <- matrix(eval(parse(text = paste0("c(", rfunc_calls, ")"))), ncol = k)
  X2 <- matrix(eval(parse(text = paste0("c(", rfunc_calls, ")"))), ncol = k)
  colnames(X1) <- colnames(X2) <- pars
  
  # Listen der Sensitivitaetsindizes (Haupteffekt, total) zu den
  # interessierenden Zeitpunkten:
  S <- T <- matrix(nrow = 1 + k, ncol = timesNum)
  
  # Durchfuehrung der Sensitivitaetsanalyse mit den Funktionen aus dem Paket
  # "sensitivity":
  if(method == "jansen"){
    x <- soboljansen_matrix(model = model_fit, X1, X2, nboot = nboot)
  } else if(method == "martinez"){
    x <- sobolmartinez_matrix(model = model_fit, X1, X2, nboot = nboot)
  }
  
  # Verarbeitung der Ergebnisse:
  ST_original <- sapply(x$ST_by_col, function(ST_col){
    c(ST_col$S[, 1], ST_col$T[, 1])
  })
  # ST_original wieder in 2 Matrizen S und T aufspalten:
  S <- rbind(times, ST_original[1:k, ])
  T <- rbind(times, ST_original[(k+1):(2*k), ])
  rownames(S) <- rownames(T) <- c("time", pars)
  
  # Rueckgabe:
  res <- list(S = S, T = T, method = method, y_idx = y_idx)
  class(res) <- "sobolRes_ats"
  return(res)
}
