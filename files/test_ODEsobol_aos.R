# Testing ODEsobol_aos():

devtools::load_all()

FHNmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dVoltage <- s * (Voltage - Voltage^3 / 3 + Current)
    dCurrent <- - 1 / s *(Voltage - a + b * Current)
    
    return(list(c(dVoltage, dCurrent)))
  })
}

FHNyini  <- c(Voltage = -1, Current = 1)
FHNtimes <- seq(0.1, 50, by = 5)

# Martinez ohne Parallelisierung:

system.time(
  FHNres_aos <- ODEsobol_aos(mod = FHNmod,
                             pars = c("a", "b", "s"),
                             yini = FHNyini,
                             times = FHNtimes,
                             ode_method = "adams",
                             ode_parallel = FALSE,
                             seed = 2015,
                             n = 10,
                             rfuncs = c("runif", "runif", "rnorm"),
                             rargs = c(rep("min = 0.18, max = 0.22", 2),
                                       "mean = 3, sd = 0.2 / 3"),
                             method = "martinez",
                             nboot = 0)
)
# Mit ode_method = "adams":
# User      System verstrichen
# 4.27        0.02        4.36
# Mit "FHNtimes <- seq(0.1, 50, by = 5)":
# User      System verstrichen 
# 7.38        0.00        7.40

# Jansen ohne Parallelisierung:

system.time(
  FHNres_aos_jansen <- ODEsobol_aos(mod = FHNmod,
                                    pars = c("a", "b", "s"),
                                    yini = FHNyini,
                                    times = FHNtimes,
                                    ode_method = "adams",
                                    seed = 2015,
                                    n = 10,
                                    rfuncs = c("runif", "runif", "rnorm"),
                                    rargs = c(rep("min = 0.18, max = 0.22", 2),
                                              "mean = 3, sd = 0.2 / 3"),
                                    method = "jansen",
                                    nboot = 0)
)
# User      System verstrichen 
# 7.37        0.00        7.44

# save(FHNres_aos, file = "test_ODEsobol_aos.RData")

# pdf("test_ODEsobol_aos.pdf", width = 10, height = 7)
plot(FHNres_aos, y_plot = "Voltage", type = "l", legendPos = "topright")
plot(FHNres_aos_jansen, y_plot = "Voltage", type = "l", legendPos = "topright")
# dev.off()

# Martinez mit Parallelisierung und n = 1000:

system.time(
  FHNres_aos_par <- ODEsobol_aos(mod = FHNmod,
                                 pars = c("a", "b", "s"),
                                 yini = FHNyini,
                                 times = FHNtimes,
                                 ode_method = "adams",
                                 ode_parallel = TRUE,
                                 ode_parallel_ncores = 4,
                                 seed = 2015,
                                 n = 1000,
                                 rfuncs = c("runif", "runif", "rnorm"),
                                 rargs = c(rep("min = 0.18, max = 0.22", 2),
                                           "mean = 3, sd = 0.2 / 3"),
                                 method = "martinez",
                                 nboot = 0)
)
# User      System verstrichen 
# 8.38        1.51      425.36

plot(FHNres_aos_par, y_plot = "Voltage", type = "l", legendPos = "topright")
