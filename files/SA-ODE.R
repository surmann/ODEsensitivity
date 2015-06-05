##--------------------------------------------------------------------##
## FOR1511: R-Code                                                    ##
##                                                                    ##
##   SA for ODE                                                       ##
##--------------------------------------------------------------------##

options(help_type = "html")
source("Package-deSolve.R")
source("Package-sensitivity.R")
source("Package-ODEnetwork.R")
dev.off()
ls()


##----------------------------------------------------------------------
## Einfachste vorstellbare DGL:   y'(t) = alpha * y(t)
##----------------------------------------------------------------------

# Modell aufstellen fuer die DGL :
dglFun <- function(X, my.t) {
  exp(X[, 1] * my.t) + X[, 2] * my.t + (X[, 1] * X[, 3]) / 1e4
  ## exp(X[, 1] * my.t)            # natuerlich die wahre Lsg.
}

# SA durchfuehren (somit kann my.t uebergeben werden):
sobol2007(model = dglFun, X1[, 1:3], X2[, 1:3], my.t = 1)
(res <- morris(model = dglFun, factors = 2, r = 25, my.t = 0,
               design = list(type = "oat", levels = 100, grid.jump = 1)))
(res <- morris(model = dglFun, factors = 2, r = 25, my.t = 17,
               design = list(type = "oat", levels = 100, grid.jump = 1)))
res$X
res$y
res$ee
plot(res)

# Fuehre fuer diesen einfachen Fall eine SA zu verschiedenen Zeitpunkten
# durch:
t.vec <- 0:10
x.fun <- function(my.t) {
  morris(model = dglFun, factors = 2, r = 100, my.t = my.t,
         design = list(type = "oat", levels = 100, grid.jump = 1))
}
one.run <- function(x.fun, my.t) {
  x <- x.fun(my.t)
  # analog zur Hilfeseite von morris()/ hoestgradig primitiv:
  k <- ncol(x$ee)
  mu <- mu.star <- sigma <- numeric(k)
  for(i in 1:k) {
    mu[i]      <- mean(x$ee[, i])
    mu.star[i] <- mean(abs(x$ee[, i]))
    sigma[i]   <- sd(x$ee[, i])
  }
  # Ergebnisse:
  res <- c(my.t, mu, mu.star, sigma)
  names(res) <- c("time",
                  paste("mu", 1:k, sep = ""),
                  paste("mu.star", 1:k, sep = ""),
                  paste("sigma", 1:k, sep = ""))
  return(res)
}
set.seed(2015)
res <- sapply(t.vec, one.run, x.fun = x.fun)

## # Visualisieren von mu.star und sigma ueber die Zeit:
## plot.inone <- function(res, ...) {
##   par(mar=c(5,4,4,5)+.1)
##   plot(t.vec, y = res[2, ], type = "b", col = "red", lwd = 2,
##        main = "mu.star und sigma in Abhaengigkeit von t",
##        xlab = "Zeit t", ylab = "mu.star")
##   par(new = TRUE)
##   plot(t.vec, y = res[3, ], type = "b", col = "darkgreen", lwd = 2,
##        xlab = "", ylab = "", yaxt = "n")
##   axis(4)
##   mtext("sigma", side = 4, line = 3)
##   legend("topleft", legend = c("mu.star", "sigma"),
##          col = c("red", "darkgreen"), lty = 1, bg = "white")
## }
## plot.inone(res)
## # Einen exponentiellen Verlauf wuerde man erwarten, da die Bedeutung
## # von alpha mit der Zeit exponentiell anwaechst.

plot.intwo <- function(res, ...) {
  t.vec <- res[1, ]
  k     <- (nrow(res) - 1) / 3
  my.cols <- rainbow(k)
  # Teile den Plot in mu.star und sigma:
  par(mfrow = c(1, 2))
  # mu.star:
  plot(t.vec, y = res[2, ], type = "b", col = my.cols[1], lwd = 1,
       main = "mu.star in Abhaengigkeit von t",
       xlab = "Zeit t", ylab = "mu.star")
  if(k > 1) {
    j <- 1
    for(i in (2*k):(3*k-1)) {
      lines(t.vec, y = res[i, ], col = my.cols[j], lwd = 1, type = "b")
      j <- j + 1
    }
  }
  legend("topleft", legend = paste("X", 1:k, sep = ""), lty = 1,
         col = my.cols)
  # sigma:
  plot(t.vec, y = res[2+2*k, ], type = "b", col = my.cols[1], lwd = 1,
       main = "sigma in Abhaengigkeit von t",
       xlab = "Zeit t", ylab = "sigma")
  if(k > 1) {
    j <- 1
    for(i in (3*k):(4*k-1)) {
      lines(t.vec, y = res[i, ], col = my.cols[j], lwd = 1, type = "b")
      j <- j + 1
    }
  }
  legend("topleft", legend = paste("X", 1:k, sep = ""), lty = 1,
         col = my.cols)
  par(mfrow = c(1, 1))
}
plot.intwo(res)

# Visualisieren mithilfe von Trajektoren:
plot.trajectories <- function(res, ...) {
  t.vec <- res[1, ]
  k     <- (nrow(res) - 1) / 3
  my.cols <- rainbow(k)
  # Zeichne Trajektoren:
  plot(x = res[k+2, ], y = res[2+2*k, ], type = "b", col = my.cols[1], lwd = 1,
       main = "Trajektoren: mu.star gegen sigma",
       xlab = "mu.star", ylab = "sigma")
  if(k > 1) {
    j <- 1
    for(i in (2*k):(3*k-1)) {
      lines(x = res[i, ], y = res[i+k, ], col = my.cols[j], lwd = 1, type = "b")
      j <- j + 1
    }
  }
  legend("topleft", legend = paste("X", 1:k, sep = ""), lty = 1,
         col = my.cols)
}
plot.trajectories(res)


##----------------------------------------------------------------------
## Sobol SA for ODEs
##----------------------------------------------------------------------

# Fuehre eine SA fuer die FHN-DGL und andere durch. Nutze Sobol.

#' @title Sobol SA for ODEs
#'
#' @description
#' \code{ODEsobol} performs a Sensitivity Analysis for
#' oridnary differential equations using
#' the variance-based Sobol method.
#'
#' @param mod model to examine
#' @param pars vector of input variable names
#' @param yini vector of initial values
#' @param times points of time at which the SA should be executed
#' @param seed Seed
#' @param n Parameter \code{n} used in \link{\code{sobol2007}}
#' @param trafo function to transform \code{z > 1} output variables to
#'   IR [only needed, if \code{z > 1}].
#'
#' @return list of Sobol SA results (i.e. 1st order \code{S} and total
#'   sensitivity indices \code{T}) for every point of time out of the
#'   \code{times} vector.
#'
#' @export
#' @import checkmate, deSolve, sensitivity, boot
#'

ODEsobol <- function(mod = LVmod,
                     pars = c("rIng", "rGrow", "rMort", "assEff", "K"),
                     yini = LVyini,
                     times = seq(1, 100, 5),
                     seed = 2015,
                     n = 1000,
                     trafo = function(Y) rowSums(Y^2)) {

  ##### Plausibilitaet #################################################
  assertNumeric(yini)

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
    trafo(res)
  }

  ##### Sensitivitaet ##################################################
  X1 <- data.frame(matrix(runif(k * n), nrow = n))
  X2 <- data.frame(matrix(runif(k * n), nrow = n))
  colnames(X1) <- colnames(X2) <- pars
  # Listen der Sensitivitaetsindizes (Haupteffekt, total) zu den
  # interessierenden Zeitpunkten:
  S <- T <- matrix(nrow = 1 + k, ncol = timesNum)
  # Durchlaufe alle Zeitpunkte und bestimme die Sensitivitaet:
  for(i in 1:timesNum) {
    res <- sobol2007(model = modFun, X1, X2, pot = times[i])
    S[, i] <- c(times[i], res$S[, 1])
    T[, i] <- c(times[i], res$T[, 1])
  }
  rownames(S) <- rownames(T) <- c("time", pars)

  # Rueckgabe:
  res <- list(S = S, T = T)
  setClasses(res, "sobolRes")
  return(res)
}

LVres <- ODEsobol(n = 50)



##----------------------------------------------------------------------
## Plot der Sobol SA-Ergebnisse
##----------------------------------------------------------------------

#' @title
#' Plotting the results of Sobol SA
#'
#' @description
#' \code{plot.sobolRes} plots the results of Sobol SA.
#'
#' @details
#' 1st order and total Sobol SA indices are plotted for each input
#' parameter against time.
#'
#' @param res resulting output of \link{\code{ODEsobol}}.
#' @param ... additional arguments.
#'
#' @return NULL
#'
#' @method plot sobolRes
#'
#' @examples
#' # none.
#'
#' @author Stefan Theers
#' @seealso \link[sensitivity]{sobol}, \link[sensitivity]{sobol2007}
#'
#' @export
#' @import checkmate
#'
plot.sobolRes <- function(res) {
  k <- nrow(res$S) - 1
  pars <- rownames(res$S)[- 1]
  parsCols <- rainbow(k)
  # Extrema SA Indizes:
  minMaxS <- c(0.95 * min(res$S[-1, ]), 1.05 * max(res$S[-1, ]))
  minMaxT <- c(0.95 * min(res$T[-1, ]), 1.05 * max(res$T[-1, ]))

  par(mfrow = c(1, 2))

  ##### 1st order SA indices ###########################################
  # Plot fuer ersten Parameter:
  plot(x = res$S[1, ], y = res$S[2, ],
       main = "1st order Sobol SA indices", xlab = "time",
       ylab = "1st order Sobol SA indices",
       type = "b", col = parsCols[1], ylim = minMaxS)
  # Plots fuer alle weiteren Parameter:
  for(i in 2:k) {
    lines(x = res$S[1, ], y = res$S[i + 1, ],
          type = "b", col = parsCols[i])
  }
  # Legende:
  legend("topleft", legend = pars, col = parsCols, bg = "white",
         lty = 1, pch = 1)

  ##### total SA indices ###############################################
  # Plot fuer ersten Parameter:
  plot(x = res$T[1, ], y = res$T[2, ],
       main = "Total Sobol SA indices", xlab = "time",
       ylab = "Total Sobol SA indices",
       type = "b", col = parsCols[1], ylim = minMaxT)
  # Plots fuer alle weiteren Parameter:
  for(i in 2:k) {
    lines(x = res$T[1, ], y = res$T[i + 1, ],
          type = "b", col = parsCols[i])
  }
  # Legende:
  legend("topleft", legend = pars, col = parsCols, bg = "white",
         lty = 1, pch = 1)

  par(mfrow = c(1, 1))
  return(invisible(TRUE))
}

plot.sobolRes(LVres)



##----------------------------------------------------------------------
## weitere Beispiele
##----------------------------------------------------------------------

# FitzHugh-Nagumo:
FHNres <- ODEsobol(mod = FHNmod,
                   pars = c("a", "b", "s"),
                   yini = FHNyini,
                   times = seq(1, 100, 5),
                   seed = 2015,
                   n = 50,
                   trafo = function(Y) rowSums(Y^2))
plot.sobolRes(FHNres)

## # ODEnetwork/ Oscillator:
## OSres <- ODEsobol(mod = OSmod,
##                   pars = c("cTime", "cState", "cParameters"),
##                   yini = OSyini,
##                   times = 1:3,
##                   seed = 2015,
##                   n = 50,
##                   trafo = function(Y) Y)
## plot.sobolRes(OSres)


