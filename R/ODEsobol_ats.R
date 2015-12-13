#' @title Sobol SA for ODEs
#'
#' @description
#' \code{ODEsobol_ats} performs a sensitivity analysis for
#' ordinary differential equations at all timepoints simultaneously using
#' the variance-based Sobol method and 
#' \code{\link[sensitivity]{soboljansen_matrix}}.
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
#' @param y_idx [\code{integer(1)}]\cr
#'   index of the output variable to be analyzed. Defaults to 1, so the first
#'   output variable is used.
#' @param seed [\code{numeric(1)}]\cr
#'   seed.
#' @param n [\code{integer(1)}]\cr
#'   parameter \code{nboot} used in \code{\link[sensitivity]{soboljansen}},
#'   i.e. the number of bootstrap replicates. Can also be zero for no 
#'   bootstrapping.
#'
#' @return list of Sobol SA results (i.e. 1st order sensitivity indices
#'   \code{S} and total sensitivity indices \code{T}) for every point of
#'   time of the \code{times} vector, of class \code{sobolRes_ats}.
#'
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
#' FHNpars  <- c(a = 0.2,     # parameter a
#'               b = 0.3,     # parameter b
#'               s = 3)       # parameter s (= c in the original notation)
#'
#' FHNyini  <- c(Voltage = -1, Current = 1)
#' FHNtimes <- seq(0.1, 20, by = 0.5)
#'
#' FHNres <- ODEsobol_ats(mod = FHNmod,
#'                        pars = names(FHNpars),
#'                        yini = FHNyini,
#'                        times = FHNtimes,
#'                        y_idx = 1,            # only Voltage
#'                        seed = 2015,
#'                        n = 10)               # use n >> 10!
#'
#' @seealso \code{\link[sensitivity]{sobol}},
#'   \code{\link[sensitivity]{soboljansen_matrix}},
#'   \code{\link{plot.sobolRes_ats}}
#'
#' @export
#' @import
#'   checkmate
#'   deSolve
#'   sensitivity
#'   boot
#'

ODEsobol_ats <- function(mod,
                         pars,
                         yini,
                         times,
                         y_idx = 1,
                         seed = 2015,
                         n = 1000) {

  ##### Check input #################################################
  ## stopifnot(!missing(...))
  assertFunction(mod)
  assertCharacter(pars)
  assertNumeric(yini)
  assertNumeric(times, lower = 0, finite = TRUE, unique = TRUE)
  times <- sort(times)
  stopifnot(!any(times == 0))
  assertIntegerish(y_idx)
  assertNumeric(seed)
  assertIntegerish(n)

  ##### Vorarbeiten ####################################################
  set.seed(seed)
  # Anzahl Parameter:
  k <- length(pars)
  # Anzahl Outputgroessen:
  z <- length(yini)
  # Anzahl Zeitpunkte von Interesse:
  timesNum <- length(times)
  # Forme DGL-Modell um, sodass fuer soboljansen_matrix()-Argument 
  # "model_matrix" passend:
  model_fit <- function(X){
    # X   - (nxk)-Matrix mit den n einzugebenden Parameter-Konstellationen
    #       als Zeilen
    colnames(X) <- pars
    res <- apply(X, 1, function(x){
      ode(yini, times = c(0, times), 
          mod, parms = x)[2:(timesNum + 1), y_idx + 1]
    })
    # (Each row in "res" represents one timepoint.)
    
    if(timesNum == 1){
      # Correction needed if timesNum == 1:
      res <- matrix(res)
    } else{
      # Transpose the results matrix, so columns represent timepoints:
      res <- t(res)
    }
    return(res)
  }
  
  ##### Sensitivitaet ##################################################
  X1 <- data.frame(matrix(runif(k * n), nrow = n))
  X2 <- data.frame(matrix(runif(k * n), nrow = n))
  colnames(X1) <- colnames(X2) <- pars
  # Listen der Sensitivitaetsindizes (Haupteffekt, total) zu den
  # interessierenden Zeitpunkten:
  S <- T <- matrix(nrow = 1 + k, ncol = timesNum)
  
  x <- soboljansen_matrix(model_matrix = model_fit, X1, X2)
  ST_original <- sapply(x$ST_by_col, function(ST_col){
    c(ST_col$S[, 1], ST_col$T[, 1])
  })
  
  # ST_original wieder in 2 Matrizen S und T aufspalten:
  S <- rbind(times, ST_original[1:k, ])
  T <- rbind(times, ST_original[(k+1):(2*k), ])
  rownames(S) <- rownames(T) <- c("time", pars)
  
  # Rueckgabe:
  res <- list(S = S, T = T)
  class(res) <- "sobolRes_ats"
  return(res)
}
