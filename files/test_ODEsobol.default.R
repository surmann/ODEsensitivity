# Testing ODEsobol.default():

devtools::load_all()

FHNmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dVoltage <- s * (Voltage - Voltage^3 / 3 + Current)
    dCurrent <- - 1 / s *(Voltage - a + b * Current)
    
    return(list(c(dVoltage, dCurrent)))
  })
}
FHNstate  <- c(Voltage = -1, Current = 1)
FHNtimes <- seq(0.1, 50, by = 5)

# Martinez without parallelization:
system.time(
  FHNres_martinez <- ODEsobol(mod = FHNmod,
                              pars = c("a", "b", "s"),
                              state_init = FHNstate,
                              times = FHNtimes,
                              seed = 2015,
                              n = 10,
                              rfuncs = c("runif", "runif", "rnorm"),
                              rargs = c(rep("min = 0.18, max = 0.22", 2),
                                        "mean = 3, sd = 0.2 / 3"),
                              sobol_method = "martinez",
                              ode_method = "adams",
                              ode_parallel = FALSE,
                              ode_parallel_ncores = NA)
)
# User      System verstrichen 
# 7.72        0.03        7.91 

# Jansen without parallelization:
system.time(
  FHNres_jansen <- ODEsobol(mod = FHNmod,
                            pars = c("a", "b", "s"),
                            state_init = FHNstate,
                            times = FHNtimes,
                            seed = 2015,
                            n = 10,
                            rfuncs = c("runif", "runif", "rnorm"),
                            rargs = c(rep("min = 0.18, max = 0.22", 2),
                                      "mean = 3, sd = 0.2 / 3"),
                            sobol_method = "jansen",
                            ode_method = "adams",
                            ode_parallel = FALSE,
                            ode_parallel_ncores = NA)
)
# User      System verstrichen 
# 7.53        0.00        7.57 

# Martinez with parallelization, n = 1000 and default values for "rfuncs" and
# "rargs":
system.time(
  FHNres_default <- ODEsobol(mod = FHNmod,
                             pars = c("a", "b", "s"),
                             state_init = FHNstate,
                             times = FHNtimes,
                             seed = 2015,
                             n = 1000,
                             rfuncs = "runif",
                             rargs = "min = 0.1, max = 1.5",
                             sobol_method = "martinez",
                             ode_method = "adams",
                             ode_parallel = TRUE,
                             ode_parallel_ncores = 2)
)
# User      System verstrichen 
# 0.47        0.30      148.73

# Martinez with parallelization, n = 1000 and changed arguments "rfuncs" and
# "rargs":
system.time(
  FHNres <- ODEsobol(mod = FHNmod,
                     pars = c("a", "b", "s"),
                     state_init = FHNstate,
                     times = FHNtimes,
                     seed = 2015,
                     n = 1000,
                     rfuncs = c("runif", "rnorm", "rexp"),
                     rargs = c("min = 0.18, max = 0.22", 
                               "mean = 0.2, sd = 0.2 / 3",
                               "rate = 1 / 3"),
                     sobol_method = "martinez",
                     ode_method = "adams",
                     ode_parallel = TRUE,
                     ode_parallel_ncores = 2)
)
# User      System verstrichen 
# 0.61        0.29      400.08

# save(FHNres, file = "test_ODEsobol.default.RData")

# pdf("test_ODEsobol.default.pdf", width = 10, height = 7)
# Checking defaults:
plot(FHNres_martinez)
plot(FHNres_jansen)
plot(FHNres)
# Custom arguments:
plot(FHNres, state_plot = "Current", legendPos = "topleft")
# Palette "Dark2" from the package "RColorBrewer" with some 
# additional colors:
my_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
             "#E6AB02", "#A6761D", "#666666", "black", "firebrick",
             "darkblue", "darkgreen")
plot(FHNres, state_plot = "Current", colors_pars = my_cols)
# dev.off()
