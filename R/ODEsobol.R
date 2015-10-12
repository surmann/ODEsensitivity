#' @title Sobol SA for ODEs
#'
#' @description
#' \code{ODEsobol} performs a sensitivity analysis for
#' ordinary differential equations using
#' the variance-based Sobol method.
#'
#' @param mod [\code{function(Time, State, Pars)}]\cr
#'   model to examine, cf. example below.
#' @param pars [\code{numeric(k)}]\cr
#'   vector of \code{k} input variable names.
#' @param yini [\code{numeric(z)}]\cr
#'   vector of \code{z} initial values.
#' @param times [\code{numeric}]\cr
#'   points of time at which the SA should be executed
#'   (vector of arbitrary length). Also the
#'   first point of time must be positive.
#' @param seed [\code{numeric(1)}]\cr
#'   seed.
#' @param n [\code{integer(1)}]\cr
#'   parameter \code{nboot} used in \code{\link{soboljansen}},
#'   i.e. the number of bootstrap replicates.
#' @param trafo [\code{function}]\cr
#'   function to transform \code{z > 1} output variables to
#'   IR [only needed, if \code{z > 1}]. Must be able to deal with a
#'   matrix.
#' @param ncores [\code{integer(1)}]\cr
#'   number of processor cores to be used for calculating the sensitivity
#'   indices. Must be between 1 and 4.
#'
#' @return list of Sobol SA results (i.e. 1st order sensitivity indices
#'   \code{S} and total sensitivity indices \code{T}) for every point of
#'   time of the \code{times} vector, of class \code{sobolRes}.
#'
#' @examples
#' ##### FitzHugh-Nagumo equations (Ramsay et al, 2007)
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
#' FHNpars  <- c(a = 0.2,     # paramter a
#'               b = 0.3,     # paramter b
#'               s = 3)       # paramter s (= c in the original notation)
#'
#' FHNyini  <- c(Voltage = -1, Current = 1)
#' FHNtimes <- seq(0.1, 20, by = 0.5)
#'
#' FHNres <- ODEsobol(mod = FHNmod,
#'                    pars = names(FHNpars),
#'                    yini = FHNyini,
#'                    times = FHNtimes,
#'                    seed = 2015,
#'                    n = 10,                        # use n >> 10!
#'                    trafo = function(Y) Y[, 1],    # voltage only
#'                    ncores = 4)
#'
#' @seealso \code{\link[sensitivity]{sobol}},
#'   \code{\link[sensitivity]{soboljansen}},
#'   \code{\link{plot.sobolRes}}
#'
#' @export
#' @import
#'   checkmate
#'   deSolve
#'   sensitivity
#'   boot
#'   parallel
#'   BBmisc
#'

ODEsobol <- function(mod,
                     pars,
                     yini,
                     times,
                     seed = 2015,
                     n = 1000,    # default eigentlich 1000
                     trafo = function(Y) rowSums(Y^2),
                     ncores = 1) {

  ##### Plausibilitaet #################################################
  ## stopifnot(!missing(...))
  assertFunction(mod)
  assertCharacter(pars)
  assertNumeric(yini)
  assertNumeric(times, lower = 0, finite = TRUE, unique = TRUE)
  times <- sort(times)
  stopifnot(!any(times == 0))
  assertNumeric(seed)
  assertIntegerish(n)
  assertFunction(trafo)
  notOk <- !testVector(trafo(matrix(1:30, nrow = 6)), len = 6)
  if(notOk)
    stop("Make sure that trafo() transforms matrices to suitable vectors!")
  assertIntegerish(ncores, lower = 1L, upper = 4L)

  ##### Vorarbeiten ####################################################
  set.seed(seed)
  # Anzahl Parameter:
  k <- length(pars)
  # Anzahl Outputgroessen:
  z <- length(yini)
  # Anzahl Zeitpunkte von Interesse:
  timesNum <- length(times)
  # Umformen DGL-Modell, sodass fuer soboljansen()-Argument model passend
  modFun <- function(X, pot) {
    # X   - (nxk)-Matrix
    # pot - point of time
    colnames(X) <- pars
    res <-
      t(apply(X, 1, function(x)
              ode(yini, times = c(0, pot), mod, parms = x)[2, 2:(z+1)]))
    # Transformation der Output-Variablen nach IR:
    trafo(res)
    ## res <- trafo(res)
    ## # Das "BE CAREFUL!" aus der R-Doku zu sobol2007() beachten, d.h.
    ## # Zentrieren:
    ## res - mean(res)
  }

  ##### Sensitivitaet ##################################################
  X1 <- data.frame(matrix(runif(k * n), nrow = n))
  X2 <- data.frame(matrix(runif(k * n), nrow = n))
  colnames(X1) <- colnames(X2) <- pars
  # Listen der Sensitivitaetsindizes (Haupteffekt, total) zu den
  # interessierenden Zeitpunkten:
  S <- T <- matrix(nrow = 1 + k, ncol = timesNum)

  # Calculates SA indices S and T for the i-th point of time:
  STForPot <- function(i) {
    pot <- times[i]
    res <- soboljansen(model = modFun, X1, X2, pot = times[i])
    return(c(res$S[, 1], res$T[, 1]))
  }

  # Durchlaufe alle Zeitpunkte und bestimme die Sensitivitaet:
  ## # ohne Parallelisierung:
  ## for(i in 1:timesNum) {
  ##   res <- soboljansen(model = modFun, X1, X2, pot = times[i])
  ##   S[, i] <- c(times[i], res$S[, 1])
  ##   T[, i] <- c(times[i], res$T[, 1])
  ## }
  cl <- makeCluster(rep("localhost", ncores), type = "SOCK")
  clusterSetRNGStream(cl)
  clusterExport(cl, list("mod", "modFun", "X1", "X2", "times",
                         "timesNum", "pars", "yini", "z", "STForPot",
                         "S", "T", "soboljansen", "ode", "trafo"),
  ## clusterExport(cl, list(ls(), "soboljansen", "ode"),
    envir = environment())
  res <- parSapply(cl, 1:timesNum, STForPot)
  stopCluster(cl)

  # res wieder in 2 Matrizen S und T aufspalten:
  S <- rbind(times, res[1:k, ])
  T <- rbind(times, res[(k+1):(2*k), ])
  rownames(S) <- rownames(T) <- c("time", pars)

  # Rueckgabe:
  res <- setClasses(list(S = S, T = T), "sobolRes")
  return(res)
}
