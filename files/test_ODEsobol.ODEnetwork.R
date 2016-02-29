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
LFOpars <- c("k.1", "k.2", "k.1.2")
LFOtimes <- seq(0.01, 20, by = 0.1)

devtools::load_all(".")

# Martinez without parallelization:
system.time(
  LFOres_martinez <- ODEsobol(lfonet, 
                              LFOpars, 
                              LFOtimes, 
                              seed = 2015, 
                              n = 10,
                              rfuncs = "runif",
                              rargs = "min = 0.001, max = 6",
                              sobol_method = "martinez",
                              ode_method = "adams",
                              ode_parallel = FALSE,
                              ode_parallel_ncores = NA)
)
# User      System verstrichen 
# 2.29        0.02        2.32

# Jansen without parallelization:
system.time(
  LFOres_jansen <- ODEsobol(lfonet, 
                            LFOpars, 
                            LFOtimes, 
                            seed = 2015, 
                            n = 10,
                            rfuncs = "runif",
                            rargs = "min = 0.001, max = 6",
                            sobol_method = "jansen",
                            ode_method = "adams",
                            ode_parallel = FALSE,
                            ode_parallel_ncores = NA)
)
# User      System verstrichen 
# 1.11        0.00        1.21

# Martinez with parallelization, n = 1000 and "rfuncs" and "rargs" being of
# length 1:
system.time(
  LFOres_1 <- ODEsobol(lfonet, 
                       LFOpars, 
                       LFOtimes, 
                       seed = 2015, 
                       n = 1000,
                       rfuncs = "runif",
                       rargs = "min = 0.001, max = 6",
                       sobol_method = "martinez",
                       ode_method = "adams",
                       ode_parallel = TRUE,
                       ode_parallel_ncores = 2)
)
# User      System verstrichen 
# 6.55        0.71       61.95
# (A warning is thrown, concerning the state variables "v.1" and "v.2". Assuming
# that we are only interested in a sensitivity analysis of "x.1" and "x.2", this
# warning is ignored.)

# Martinez with parallelization, n = 1000 and "rfuncs" and "rargs" being of the
# same length as "pars":
system.time(
  LFOres <- ODEsobol(lfonet, 
                     LFOpars, 
                     LFOtimes, 
                     seed = 2015, 
                     n = 1000,
                     rfuncs = c("runif", "rnorm", "rexp"),
                     rargs = c("min = 0.001, max = 6",
                               "mean = 3, sd = 0.5",
                               "rate = 1 / 3"),
                     sobol_method = "martinez",
                     ode_method = "adams",
                     ode_parallel = TRUE,
                     ode_parallel_ncores = 2)
)
# User      System verstrichen 
# 5.32        0.46       58.78
# (A warning is thrown, concerning the state variables "v.1" and "v.2". Assuming
# that we are only interested in a sensitivity analysis of "x.1" and "x.2", this
# warning is ignored.)

# save(LFOres, file = "test_ODEsobol.ODEnetwork.Rdata")

# pdf("test_ODEsobol.ODEnetwork.pdf", width = 12, height = 10)
# Checking defaults:
plot(LFOres_martinez)
plot(LFOres_jansen)
plot(LFOres)
# Custom arguments:
plot(LFOres, state_plot = "x.2", legendPos = "topleft")
# Custom colors:
my_cols <- c("firebrick", "chartreuse3", "dodgerblue")
plot(LFOres, state_plot = "x.2", colors_pars = my_cols)
# Checking the passing of arguments:
plot(LFOres, state_plot = "x.2", type = "p", colors_pars = my_cols, 
     main_title = "Big Title", legendPos = "outside", cex.axis = 2, cex = 4, 
     main = "Small Title", cex.main = 0.5)
# dev.off()
