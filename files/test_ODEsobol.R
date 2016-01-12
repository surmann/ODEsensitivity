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

ODEtimes <- seq(0.01, 20, by = 0.1)
ODEbinf <- c(rep(0.001, 9), -4, 0.001, -10)
ODEbsup <- c(2, 1.5, 6, 4, 6, 10, 2, 1.5, 6, -0.001, 6, -0.001)

devtools::load_all(".")

system.time(
  ODEres <- ODEsobol(odenet, ODEtimes, ode_method = "adams", seed = 2015, 
                     n = 10,
#                      rfuncs = c(rep("runif", length(ODEbinf) - 1), "rnorm"),
#                      rargs = c(rep("min = 0.18, max = 0.22", 2),
#                                paste0("mean = ", sd = 0.2 / 3"),
                     method = "martinez",
                     nboot = 0)
)
#  User      System verstrichen 
# 31.22        0.05       32.04

# save(ODEres, file = "test_ODEsobol.Rdata")

# pdf("test_ODEsobol.pdf", width = 12, height = 10)
# Checking defaults:
plot(ODEres)
# Standard (separate plots):
plot(ODEres, state_plot = "x.2", type = "sep", legendPos = "topleft")
plot(ODEres, state_plot = "x.2", type = "sep", legendPos = "outside")
# Palette "Dark2" from the package "RColorBrewer" with some 
# additional colors:
my_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
             "#E6AB02", "#A6761D", "#666666", "black", "firebrick",
             "darkblue", "darkgreen")
plot(ODEres, state_plot = "x.2", type = "sep", colors_pars = my_cols, 
     legendPos = "outside")

# Trajectories:
plot(ODEres, state_plot = "x.2", type = "trajec", legendPos = "topleft")
plot(ODEres, state_plot = "x.2", type = "trajec", legendPos = "outside")
plot(ODEres, state_plot = "x.2", type = "trajec", colors_pars = my_cols, 
     legendPos = "outside")
# dev.off()
