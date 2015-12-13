# Testing ODEsobol_ats():

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
FHNres_ats <- ODEsobol_ats(mod = FHNmod,
                           pars = names(FHNpars),
                           yini = FHNyini,
                           times = FHNtimes,
                           y_idx = 1,
                           seed = 2015,
                           n = 10)
)
# User      System verstrichen 
# 1.77        0.00        1.92

# save(FHNres_ats, file = "SA-ODEsobol_ats.RData")

# pdf("SA-ODEsobol_ats.pdf", width = 10, height = 7)
plot(FHNres_ats, type = "l", legendPos = "topright")
# dev.off()
