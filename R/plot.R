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
