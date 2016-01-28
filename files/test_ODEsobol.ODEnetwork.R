# Testing ODEmorris() for objects of class "ODEnetwork":

library(ODEnetwork)
masses <- c(1, 1)
dampers <- diag(c(1, 1))
springs <- diag(c(1, 1))
springs[1, 2] <- 1
distances <- diag(c(0, 2))
distances[1, 2] <- 1
lfonet <- ODEnetwork(masses, dampers, springs, 
                     cartesian = TRUE, distances = distances)
lfonet <- setState(lfonet, c(0.5, 1), c(0, 0))
LFOpars <- c("m.1", "d.1", "k.1", "k.1.2", "m.2", "d.2", "k.2")
LFOtimes <- seq(0.01, 20, by = 0.1)
LFObinf <- rep(0.001, length(LFOpars) - 1)
LFObsup <- c(2, 1.5, 6, 6, 2, 1.5)

devtools::load_all(".")

# Martinez ohne Parallelisierung:
system.time(
  LFOres_martinez <- ODEsobol(lfonet, 
                              LFOpars, 
                              LFOtimes, 
                              seed = 2015, 
                              n = 10,
                              rfuncs = c(rep("runif", length(LFObinf)), 
                                         "rnorm"),
                              rargs = c(paste0("min = ", LFObinf, 
                                               ", max = ", LFObsup),
                                        "mean = 3, sd = 0.8"),
                              sobol_method = "martinez",
                              nboot = 0,
                              ode_method = "adams",
                              ode_parallel = FALSE,
                              ode_parallel_ncores = NA)
)
#  User      System verstrichen 
# 20.77        0.04       20.92

# Jansen ohne Parallelisierung:
system.time(
  LFOres_jansen <- ODEsobol(lfonet, 
                            LFOpars, 
                            LFOtimes, 
                            seed = 2015, 
                            n = 10,
                            rfuncs = c(rep("runif", length(LFObinf)), "rnorm"),
                            rargs = c(paste0("min = ", LFObinf, 
                                             ", max = ", LFObsup),
                                      "mean = 3, sd = 0.8"),
                            sobol_method = "jansen",
                            nboot = 0,
                            ode_method = "adams",
                            ode_parallel = FALSE,
                            ode_parallel_ncores = NA)
)
# User      System verstrichen 
# 5.20        0.00        5.36

# Martinez mit Parallelisierung und n = 1000:
system.time(
  LFOres <- ODEsobol(lfonet, 
                     LFOpars, 
                     LFOtimes, 
                     seed = 2015, 
                     n = 1000,
                     rfuncs = c(rep("runif", length(LFObinf)), 
                                "rnorm"),
                     rargs = c(paste0("min = ", LFObinf,
                                      ", max = ", LFObsup),
                               "mean = 3, sd = 0.8"),
                     sobol_method = "martinez",
                     nboot = 0,
                     ode_method = "adams",
                     ode_parallel = TRUE,
                     ode_parallel_ncores = 2)
)
#  User      System verstrichen 
# 31.25        1.53      149.10

# save(LFOres, file = "test_ODEsobol.ODEnetwork.Rdata")

# pdf("test_ODEsobol.ODEnetwork.pdf", width = 12, height = 10)
# Checking defaults:
plot(LFOres_martinez)
plot(LFOres_jansen)
plot(LFOres)
# Custom arguments:
plot(LFOres, state_plot = "x.2", legendPos = "topleft")
# Palette "Dark2" from the package "RColorBrewer" with some 
# additional colors:
my_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
             "#E6AB02", "#A6761D", "#666666", "black", "firebrick",
             "darkblue", "darkgreen")
plot(LFOres, state_plot = "x.2", colors_pars = my_cols)
# Checking the passing of arguments:
plot(LFOres, state_plot = "x.2", type = "p", colors_pars = my_cols, 
     main_title = "Big Title", legendPos = "outside", cex.axis = 2, cex = 4, 
     main = "Small Title", cex.main = 0.5)
# dev.off()
