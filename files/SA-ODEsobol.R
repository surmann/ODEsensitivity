# Make sure that ODEsobol() still works as before:

devtools::load_all()

FHNmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dVoltage <- s * (Voltage - Voltage^3 / 3 + Current)
    dCurrent <- - 1 / s *(Voltage - a + b * Current)
    
    return(list(c(dVoltage, dCurrent)))
  })
}

FHNpars  <- c(a = 0.2,     # parameter a
              b = 0.3,     # parameter b
              s = 3)       # parameter s (= c in the original notation)

FHNyini  <- c(Voltage = -1, Current = 1)
FHNtimes <- seq(0.1, 20, by = 0.5)

system.time(
FHNres <- ODEsobol(mod = FHNmod,
                   pars = names(FHNpars),
                   yini = FHNyini,
                   times = FHNtimes,
                   seed = 2015,
                   n = 10,                        # use n >> 10!
                   trafo = function(Y) Y[, 1],    # voltage only
                   ncores = 2)
)
# User      System verstrichen 
# 0.18        0.09       28.64

# save(FHNres, file = "SA-ODEsobol.RData")

# pdf("SA-ODEsobol.pdf", width = 10, height = 7)
plot(FHNres, type = "l", legendPos = "topright")
# dev.off()