##--------------------------------------------------------------------##
## FOR1511: R-Code                                                    ##
##                                                                    ##
##   sensitivity Package                                              ##
##--------------------------------------------------------------------##

options(help_type = "html")
library("sensitivity")                         # fuer SA
library("boot")                                # fuer SA, Sobol


##----------------------------------------------------------------------
## morris(.)
##----------------------------------------------------------------------

## # Testmodelle aus dem Paket:
## ??test.models

# Test case : the non-monotonic function of Morris
## ?morris
x <- morris(model = morris.fun, factors = 20, r = 4,
            design = list(type = "oat", levels = 5, grid.jump = 3))
print(x)
plot(x)
## ## Not run:
## library(rgl)
## plot3d.morris(x) # (requires the package 'rgl')
## ## End(Not run)


##----------------------------------------------------------------------
## sobol(.)
##----------------------------------------------------------------------

# Test case : the non-monotonic Sobol g-function
# The method of sobol requires 2 samples
# There are 8 factors, all following the uniform distribution
# on [0,1]
n <- 1000
X1 <- data.frame(matrix(runif(8 * n), nrow = n))
X2 <- data.frame(matrix(runif(8 * n), nrow = n))

# sensitivity analysis
x <- sobol2002(model = sobol.fun, X1, X2, nboot = 100)
print(x)
plot(x)

# 1st order indices:
x$S

# Nutze besser sobol2002, sobol2007!


##----------------------------------------------------------------------
## sobol2007(.)
##----------------------------------------------------------------------

# Test case : the non-monotonic Sobol g-function
# The method of sobol requires 2 samples
# There are 8 factors, all following the uniform distribution
# on [0,1]
## n <- 1000
## X1 <- data.frame(matrix(runif(8 * n), nrow = n))
## X2 <- data.frame(matrix(runif(8 * n), nrow = n))

# sensitivity analysis
y <- sobol2002(model = sobol.fun, X1, X2, nboot = 1000)
print(y)
plot(y)

# Haupteffekt-Sensitivitaetsindizes S_i fuer x_i:
y$S
# totaler Sensitivitaetsindex ST_i fuer x_i:
y$T

