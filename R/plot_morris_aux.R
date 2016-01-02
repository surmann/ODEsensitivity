# Helper functions for Morris-SA

##### auxiliary function: plotting mu.star and sigma separately ###########
plotSep <- function(res, pars, legendPos, ...) {
  t.vec <- res[1, ]
  k     <- (nrow(res) - 1) / 3
  my.cols <- rainbow(k)
  # mu.star:
  plot(t.vec, y = res[k + 2, ], type = "l", col = my.cols[1], lwd = 1,
       ylim = c(min(res[(k+2):(2*k+1), ], na.rm = TRUE),
                max(res[(k+2):(2*k+1), ], na.rm = TRUE)),
       ## min(max(res[(k+2):(2*k+1), ], na.rm = TRUE),  # willkuerlich!
       ##     5 * median(res[(k+2):(2*k+1), ], na.rm = TRUE)) ),
       xlab = "Time", ylab = "mu.star", ...)
  if(k > 1) {
    j <- 2
    for(i in (k+3):(2*k+1)) {
      lines(t.vec, y = res[i, ], col = my.cols[j], lwd = 1, type = "l")
      j <- j + 1
    }
  }
  legend(legendPos, legend = pars, lty = 1, col = my.cols)
  # sigma:
  plot(t.vec, y = res[2+2*k, ], type = "l", col = my.cols[1], lwd = 1,
       ylim = c(min(res[(2*k+2):(3*k+1), ], na.rm = TRUE),
                max(res[(2*k+2):(3*k+1), ], na.rm = TRUE)),
       xlab = "Time", ylab = "sigma", ...)
  if(k > 1) {
    j <- 2
    for(i in (2*k+3):(3*k+1)) {
      lines(t.vec, y = res[i, ], col = my.cols[j], lwd = 1, type = "l")
      j <- j + 1
    }
  }
  legend(legendPos, legend = pars, lty = 1, col = my.cols)
  par(mfrow = c(1, 1))
}

##### auxiliary function: plotting trajectories ######################
plotTrajectories <- function(res, pars, legendPos, main_title, ...) {
  t.vec <- res[1, ]
  k     <- length(pars)
  my.cols <- rainbow(k)
  # Zeichne Trajektoren:
  plot(x = res[k+2, ], y = res[2+2*k, ], type = "b", col = my.cols[1], lwd = 1,
       main = main_title,
       xlim = c(min(res[(k+2):(2*k+1), ], na.rm = TRUE),
                max(res[(k+2):(2*k+1), ], na.rm = TRUE)),
       ylim = c(min(res[(2*k+2):(3*k+1), ], na.rm = TRUE),
                max(res[(2*k+2):(3*k+1), ], na.rm = TRUE)),
       xlab = "mu.star", ylab = "sigma", ...)
  if(k > 1) {
    j <- 2
    for(i in (k+3):(2*k+1)) {
      lines(x = res[i, ], y = res[i+k, ], col = my.cols[j], lwd = 1, type = "b")
      j <- j + 1
    }
  }
  legend(legendPos, legend = pars, lty = 1, col = my.cols)
}