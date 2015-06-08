##--------------------------------------------------------------------##
## FOR1511: R-Code                                                    ##
##                                                                    ##
##   Saltelli test functions                                          ##
##   (Saltelli et al, 2010, Kap. 2)                                   ##
##--------------------------------------------------------------------##

options(help_type = "html")
library("sensitivity")                         # fuer SA
library("boot")                                # fuer SA, Sobol
source("Package-sensitivity.R")
dev.off()


##----------------------------------------------------------------------
## SA fuer lineare Testprobleme: model1
##----------------------------------------------------------------------

# Modell aufstellen mit (levels x k)-Matrix X:
model1.fun <- function(X) rowSums(X)
model1.fun <- function(X) {
  rowSums(X + rnorm(length(X), 0, 0.01))
}
model1.fun <- function(X) {
  k <- ncol(X)
  l <- nrow(X)   # fuer levels
  x.mean <- sapply(1:k, function(j) 3^(j - 1))
  x.var  <- 0.5 * x.mean
  X.add  <- sapply(1:k,
    function(j) {runif(l, x.mean[j] - x.var[j], x.mean[j] + x.var[j])})
  rowSums(X + X.add)
}

# SA durchfuehren:
set.seed(2015)
res <- morris(model = model1.fun, factors = 3, r = 50,
              design = list(type = "oat", levels = 100, grid.jump = 3))
# Ergebnis plotten:
plot(res, main = "Ergebnis Morris-SA fuer model1")
##plot3d.morris(res)
# hochgradig langweilig!


##----------------------------------------------------------------------
## SA fuer monotone Testprobleme: model4
##----------------------------------------------------------------------

# Modell aufstellen mit (levels x k)-Matrix X:
model2.fun <- function(X) X[, 1] + (X[, 2])^4

# SA durchfuehren:
res <- morris(model = model2.fun, factors = 3, r = 4,
              design = list(type = "oat", levels = 100, grid.jump = 3))
# Ergebnis plotten:
plot(res, main = "Ergebnis Morris-SA fuer model4")
##plot3d.morris(res)
# X3 kann voellig vernachlaessigt werden, X2 streut viel, aber mit
# geringerem Einfluss als X1 (klar: auf (0, 1)).


##----------------------------------------------------------------------
## SA fuer nicht-monotone Testprobleme: model5
##----------------------------------------------------------------------

# Modell aufstellen mit (levels x k)-Matrix X:
model5.fun <- function(X) {
  ## z.B:  X <- matrix(rep(1:6, 5), ncol = 6, byrow = TRUE)
  k <- 6
  b <- c(1.5, rep(0.9, 5))
  I.k <- prod((exp(b) - 1) / b)
  apply(X, 1, function(x) sum(b * x)) - I.k
}

# SA durchfuehren:
res <- morris(model = model5.fun, factors = 6, r = 4,
              design = list(type = "oat", levels = 100, grid.jump = 3))
# Ergebnis plotten:
plot(res, main = "Ergebnis Morris-SA fuer model5")
plot3d.morris(res)
# X1 mit groesstem Einfluss. Alle uebrigen weniger.


##----------------------------------------------------------------------
## Untersuchung vorimplementierte Saltelli-Testfunktionen
##----------------------------------------------------------------------

# Hilfeseiten zu Standard-Testfunktionen:
?sobol.fun
sobol.fun

y <- sobol2002(model = sobol.fun, X1, X2, nboot = 1e6)
# auch hier: negative ST_i !?
