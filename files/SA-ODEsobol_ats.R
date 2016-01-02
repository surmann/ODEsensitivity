# Testing ODEsobol_ats():

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

system.time(
FHNres_ats <- ODEsobol_ats(mod = FHNmod,
                           pars = c("a", "b", "s"),
                           yini = FHNyini,
                           times = FHNtimes,
                           ode_method = "adams",
                           y_idx = 1,
                           seed = 2015,
                           n = 10,
                           rfuncs = c("runif", "runif", "rnorm"),
                           rargs = c(rep("min = 0.18, max = 0.22", 2),
                                     "mean = 3, sd = 0.2 / 3"),
                           method = "martinez",
                           nboot = 0)
)
# User      System verstrichen 
# 1.77        0.00        1.92
# Mit neuen Argumenten "rfuncs" und "rargs":
# User      System verstrichen 
# 5.56        0.00        5.75
# Mit neuen Argumenten "rfuncs" und "rargs" sowie method = "jansen":
# User      System verstrichen 
# 5.55        0.00        6.26
# Mit neuen Argumenten "rfuncs" und "rargs" sowie method = "martinez":
# User      System verstrichen 
# 5.66        0.00        6.00
# Mit ode_method = "adams":
# User      System verstrichen 
# 3.88        0.00        4.11
# Mit "FHNtimes <- seq(0.1, 50, by = 5)":
# User      System verstrichen 
# 7.44        0.02        7.62

# save(FHNres_ats, file = "SA-ODEsobol_ats.RData")

# pdf("SA-ODEsobol_ats.pdf", width = 10, height = 7)
plot(FHNres_ats, type = "l", legendPos = "topright")
# dev.off()
