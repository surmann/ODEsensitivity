##### auxiliary function: plotting mu and sigma seperately ###########
plotSep <- function(res, ...) {
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

##### auxiliary function: plotting trajectories ######################
plotTrajectories <- function(res, ...) {
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

#' @title
#' Plotting the results of Morris SA
#'
#' @description
#' \code{plot} plots the results of Morris SA.
#'
#' @details
#' \code{plot} with \code{type = "sep"} plots mu.star and
#'   sigma seperately versus time.
#'
#' \code{plot} with \code{type = "trajec"} plots mu.star versus
#'   sigma for every point of time.
#'
#' @param res resulting ranking of class \code{morrisRes}.
#' @param type plot type, choose between \code{"sep"} and \code{"trajec"}.
#' @param ... additional arguments.
#'
#' @return NULL
#'
#' @method plot morrisRes
#'
#' @examples
#' # Modell aufstellen fuer die DGL y'(t) = alpha * y(t):
#' dglFun <- function(X, my.t) {
#'   exp(X[, 1] * my.t) + X[, 2] * my.t
#'   ## exp(X[, 1] * my.t)            # natuerlich die wahre Lsg.
#' }
#'
#' # Fuehre fuer diesen einfachen Fall eine SA zu verschiedenen Zeitpunkten
#' # durch:
#' t.vec <- 0:10
#' xFun <- function(my.t) {
#'   morris(model = dglFun, factors = 2, r = 100, my.t = my.t,
#'          design = list(type = "oat", levels = 100, grid.jump = 1))
#' }
#' oneRun <- function(xFun, my.t) {
#'   x <- xFun(my.t)
#'   # analog zur Hilfeseite von morris()/ hoestgradig primitiv:
#'   k <- ncol(x$ee)
#'   mu <- mu.star <- sigma <- numeric(k)
#'   for(i in 1:k) {
#'     mu[i]      <- mean(x$ee[, i])
#'     mu.star[i] <- mean(abs(x$ee[, i]))
#'     sigma[i]   <- sd(x$ee[, i])
#'   }
#'   # Ergebnisse:
#'   res <- c(my.t, mu, mu.star, sigma)
#'   names(res) <- c("time",
#'                   paste("mu", 1:k, sep = ""),
#'                   paste("mu.star", 1:k, sep = ""),
#'                   paste("sigma", 1:k, sep = ""))
#'   return(res)
#' }
#' set.seed(2015)
#' res <- setClasses(sapply(t.vec, oneRun, xFun = xFun), "morrisRes")
#'
#' # Plots:
#' plot(res, type = "sep")
#' plot(res, type = "trajec")
#'
#' @author Stefan Theers
#' @seealso \link[sensitivity]{morris}
#'
#' @export
#' @import checkmate
#'
plot.morrisRes <- function(res, type, ...) {
  if(type == "sep")    plotSep(res, ...)
  if(type == "trajec") plotTrajectories(res, ...)
  # for testing purposes:
  return(invisible(TRUE))
}
