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

system.time(
FHNres_aos <- ODEsobol_aos(mod = FHNmod,
                           pars = c("a", "b", "s"),
                           yini = FHNyini,
                           times = FHNtimes,
                           ode_method = "adams",
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
# 7.91        0.02        8.69

# save(FHNres_aos, file = "SA-ODEsobol_aos.RData")

# pdf("SA-ODEsobol_aos.pdf", width = 10, height = 7)
plot(FHNres_aos, y_idx = 1, type = "l", legendPos = "topright")
# dev.off()
