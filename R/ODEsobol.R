#' @title Sobol SA for ODEs
#'
#' @description
#' \code{ODEsobol} performs a sensitivity analysis for
#' ordinary differential equations using
#' the variance-based Sobol method.
#'
#' @param mod model to examine.
#' @param pars vector of input variable names.
#' @param yini vector of initial values.
#' @param times points of time at which the SA should be executed.
#' @param seed seed.
#' @param n parameter \code{nboot} used in \code{\link{sobol2007}},
#'   i.e. the number of bootstrap replicates.
#' @param trafo function to transform \code{z > 1} output variables to
#'   IR [only needed, if \code{z > 1}]. Must be able to deal with a
#'   matrix.
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
#'                    pars = c("a", "b", "s"),
#'                    yini = FHNyini,
#'                    times = FHNtimes,
#'                    seed = 2015,
#'                    n = 10,
#'                    trafo = function(Y) rowSums(Y^2))
#'
#' @seealso \code{\link[sensitivity]{sobol}},
#'   \code{\link[sensitivity]{sobol2007}}
#'
#' @export
#' @import
#'   checkmate
#'   deSolve
#'   sensitivity
#'   boot
#'   parallel
#'

ODEsobol <- function(mod,
                     pars,
                     yini,
                     times,
                     seed = 2015,
                     n = 1000,    # default eigentlich 1000
                     trafo = function(Y) rowSums(Y^2)) {

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

  ##### Vorarbeiten ####################################################
  # Anzahl Parameter:
  k <- length(pars)
  # Anzahl Outputgroessen:
  z <- length(yini)
  # Anzahl Zeitpunkte von Interesse:
  timesNum <- length(times)
  # Umformen DGL-Modell, sodass fuer sobol2007()-Argument model passend
  modFun <- function(X, pot) {
    # X   - (nxk)-Matrix
    # pot - point of time
    colnames(X) <- pars
    res <-
      t(apply(X, 1, function(x)
              ode(yini, times = c(0, pot), mod, parms = x)[2, 2:(z+1)]))
    # Transformation der Output-Variablen nach IR:
    res <- trafo(res)
    # Das "BE CAREFUL!" aus der R-Doku zu sobol2007() beachten, d.h.
    # Zentrieren:
    res - mean(res)
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
    res <- sobol2007(model = modFun, X1, X2, pot = times[i])
    return(c(res$S[, 1], res$T[, 1]))
  }

  # Durchlaufe alle Zeitpunkte und bestimme die Sensitivitaet:
  ## # ohne Parallelisierung:
  ## for(i in 1:timesNum) {
  ##   res <- sobol2007(model = modFun, X1, X2, pot = times[i])
  ##   S[, i] <- c(times[i], res$S[, 1])
  ##   T[, i] <- c(times[i], res$T[, 1])
  ## }
  cl <- makeCluster(rep("localhost", 4), type = "SOCK")
  clusterSetRNGStream(cl)
  clusterExport(cl, list("mod", "modFun", "X1", "X2", "times",
                         "timesNum", "pars", "yini", "z", "STForPot",
                         "S", "T", "sobol2007", "ode", "trafo"),
  ## clusterExport(cl, list(ls(), "sobol2007", "ode"),
    envir = environment())
  res <- parSapply(cl, 1:timesNum, STForPot)
  stopCluster(cl)

  # res wieder in 2 Matrizen S und T aufspalten:
  S <- rbind(times, res[1:k, ])
  T <- rbind(times, res[(k+1):(2*k), ])
  rownames(S) <- rownames(T) <- c("time", pars)

  # Rueckgabe:
  sobolRes <- function(S, T) {
    res <- list(S = S, T = T)
    # Name der Klasse:
    class(res) <- append(class(res), "sobolRes")
    return(res)
  }
  res <- sobolRes(S = S, T = T)
  return(res)
}
