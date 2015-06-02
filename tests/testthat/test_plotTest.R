context("Test of plot")

plotTest <- function() {
  # irgendeine Funktion:
  dglFun <- function(X, my.t) {
    exp(X[, 1] * my.t) + X[, 2] * my.t
    ## exp(X[, 1] * my.t)            # natuerlich die wahre Lsg.
  }

  # Fuehre fuer diesen einfachen Fall eine SA zu verschiedenen Zeitpunkten
  # durch:
  t.vec <- 0:10
  xFun <- function(my.t) {
    morris(model = dglFun, factors = 2, r = 100, my.t = my.t,
           design = list(type = "oat", levels = 100, grid.jump = 1))
  }
  oneRun <- function(xFun, my.t) {
    x <- xFun(my.t)
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

  res <- setClasses(sapply(t.vec, oneRun, xFun = xFun), "morrisRes")

  expect_true(plot(res))
  expect_true(plot(res), type = "trajec")
  expect_true(plot(res), type = "sep")
  expect_true(plot(res), main = "Something.")
}
