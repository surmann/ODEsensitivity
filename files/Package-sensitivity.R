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
set.seed(2015)
x <- morris(model = morris.fun, factors = 20, r = 10,
            design = list(type = "oat", levels = 5, grid.jump = 1))
## Eingabe:
##   factors      - k, Anzahl Inputs
##   r            - r, Anzahl Trajektorien
##   levels       - p, Anzahl Levels
##   grid.jump    - Hooehe eines Sprungs bei Levels

## ## Not run:
## library(rgl)
## plot3d.morris(x) # (requires the package 'rgl')
## ## End(Not run)


##### Plotting the results for the report ##############################
set.seed(2015)
x <- morris(model = sobol.fun, factors = 8, r = 50,
            design = list(type = "oat", levels = 10, grid.jump = 1))
print(x)
plot(x)

mu <- apply(x$ee, 2, mean)
mu.star <- apply(x$ee, 2, function(x) mean(abs(x)))
sigma <- apply(x$ee, 2, sd)

pdf("C:/Users/theers/ALLE/07_Studium/2012-2015_TU_Dortmund_Statistik_BSc/HiWi-Jobs/LehrstuhlCStat/Bericht/img/sobolfun-morris-mu.pdf",
    width = 5, height = 5)
# Plot mu vs. sigma:
plot(mu, sigma, xlim = c(min(mu), max(mu + 0.1)),
     xlab = expression(paste(mu)), ylab = expression(sigma))
text(mu, sigma, labels = paste("x", seq_along(sigma), sep = ""),
     pos = 4)
dev.off()

pdf("C:/Users/theers/ALLE/07_Studium/2012-2015_TU_Dortmund_Statistik_BSc/HiWi-Jobs/LehrstuhlCStat/Bericht/img/sobolfun-morris-mustar.pdf",
    width = 5, height = 5)
# Plot mu.star vs. sigma:
plot(mu.star, sigma, xlim = c(0, max(mu.star + 0.5)),
     xlab = expression(paste(mu, "*")), ylab = expression(sigma))
text(mu.star, sigma, labels = paste("x", seq_along(sigma), sep = ""),
     pos = 4)
dev.off()


##### Scatterplots for Sobol g-function ################################
pdf("C:/Users/theers/ALLE/07_Studium/2012-2015_TU_Dortmund_Statistik_BSc/HiWi-Jobs/LehrstuhlCStat/Bericht/img/sobolfun-x1.pdf",
    width = 10, height = 5)
set.seed(2015)
n <- 500
X1 <- data.frame(matrix(runif(8 * n), nrow = n))
y <- sobol.fun(X1)
# Scatterplot x1 vs. y:
plot(X1[, 1], y, xlab = expression(x[1]), ylab = "y")
dev.off()


##----------------------------------------------------------------------
## sobol(.)
##----------------------------------------------------------------------

# Test case : the non-monotonic Sobol g-function
# The method of sobol requires 2 samples
# There are 8 factors, all following the uniform distribution
# on [0,1]
set.seed(2015)
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
y <- sobol2007(model = sobol.fun, X1, X2, nboot = 1000)
print(y)
plot(y)

# Haupteffekt-Sensitivitaetsindizes S_i fuer x_i:
y$S
# totaler Sensitivitaetsindex ST_i fuer x_i:
y$T


##### Plotting the results for the report ##############################
pdf("C:/Users/theers/ALLE/07_Studium/2012-2015_TU_Dortmund_Statistik_BSc/HiWi-Jobs/LehrstuhlCStat/Bericht/img/sobolfun-sobol2007.pdf",
    width = 10, height = 5)
y <- soboljansen(model = sobol.fun, X1, X2, nboot = 1000, conf = 0.95)
print(y)
plot(y)
abline(h=0, lty = 2, col = "gray")
dev.off()


##----------------------------------------------------------------------
## soboljansen(.) und sobolmartinez(.)
##----------------------------------------------------------------------

# Ursprünglich wurde wegen
#  http://stats.stackexchange.com/questions/43504/interpreting-results-from-sobol-sensitivity-analysis-in-r
# die soboljansen()-Implementierung (mit gleichen Argumenten wie sobol() und
# sobol2007()) in ODEsobol() verwendet. Seit Version 1.11 
# des Pakets "sensitivity" gibt es jedoch auch die Funktion 
# sobolmartinez(), die etwas schneller arbeitet als soboljansen(). Dem 
# Anwender wird die Wahl zwischen den beiden Methoden gelassen, jedoch
# ist sobolmartinez() der Standard.
