# Testing ODEmorris() for objects of class "ODEnetwork":

library(ODEnetwork)

masses <- c(1, 1)
dampers <- diag(c(1, 1))
springs <- diag(c(1, 1))
springs[1, 2] <- 1
distances <- diag(c(0, 2))
distances[1, 2] <- 1
odenet <- ODEnetwork(masses, dampers, springs, 
                     cartesian = TRUE, distances = distances)
odenet <- setState(odenet, c(0.5, 1), c(0, 0))

ODEpars <- c("m.1", "d.1", "k.1", "k.1.2", "m.2", "d.2", "k.2")
ODEtimes <- seq(0.01, 20, by = 0.1)
ODEbinf <- rep(0.001, length(ODEpars) - 1)
ODEbsup <- c(2, 1.5, 6, 6, 2, 1.5)

devtools::load_all(".")

system.time(
  ODEres <- ODEsobol(odenet, ODEpars, ODEtimes, ode_method = "adams", 
                     seed = 2015, n = 10,
                     rfuncs = c(rep("runif", length(ODEbinf)), "rnorm"),
                     rargs = c(paste0("min = ", ODEbinf, ", max = ", ODEbsup),
                               "mean = 3, sd = 0.8"),
                     method = "martinez",
                     nboot = 0)
)
#  User      System verstrichen 
# 21.18        0.01       30.36

# save(ODEres, file = "test_ODEsobol.Rdata")

# pdf("test_ODEsobol.pdf", width = 12, height = 10)
# Checking defaults:
plot(ODEres)
# Custom arguments:
plot(ODEres, state_plot = "x.2", type = "l", legendPos = "topleft")
plot(ODEres, state_plot = "x.2", type = "l", legendPos = "outside")
# Palette "Dark2" from the package "RColorBrewer" with some 
# additional colors:
my_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
             "#E6AB02", "#A6761D", "#666666", "black", "firebrick",
             "darkblue", "darkgreen")
plot(ODEres, state_plot = "x.2", type = "l", colors_pars = my_cols, 
     legendPos = "outside")
plot(ODEres, state_plot = "x.2", type = "p", colors_pars = my_cols, 
     main_title = "Big Title", legendPos = "outside", cex.axis = 2, cex = 4, 
     main = "Small Title", cex.main = 0.5)
# dev.off()
