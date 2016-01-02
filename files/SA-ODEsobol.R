# Make sure that ODEsobol() still works as before:

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
FHNres <- ODEsobol(mod = FHNmod,
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
                   nboot = 0,
                   trafo = function(Y) Y[, 1],    # voltage only
                   ncores = 2)
)
# User      System verstrichen 
# 0.18        0.09       28.64
# Mit neuen Argumenten "rfuncs" und "rargs" sowie method = "jansen":
# User      System verstrichen 
# 0.17        0.03      102.44
# Mit neuen Argumenten "rfuncs" und "rargs" sowie method = "martinez":
# User      System verstrichen 
# 0.25        0.11       83.64
# Mit ode_method = "adams":
# User      System verstrichen 
# 0.11        0.08       53.84
# Mit "FHNtimes <- seq(0.1, 50, by = 5)":
# User      System verstrichen 
# 0.17        0.09       31.50

# save(FHNres, file = "SA-ODEsobol.RData")

# pdf("SA-ODEsobol.pdf", width = 10, height = 7)
plot(FHNres, type = "l", legendPos = "topright")
# dev.off()
