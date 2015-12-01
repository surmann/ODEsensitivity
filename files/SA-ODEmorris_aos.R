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
FHNtimes <- seq(0.1, 100, by = 10)

system.time(
FHNres_aos <- ODEmorris_aos(mod = FHNmod,
                    pars = names(FHNpars),
                    yini = FHNyini,
                    times = FHNtimes,
                    seed = 2015,
                    binf = c(0.18, 0.18, 2.8),
                    bsup = c(0.22, 0.22, 3.2),
                    r = 10,
                    design =
                      list(type = "oat", levels = 30, grid.jump = 1),
                    scale = FALSE)
)
#  User      System verstrichen 
# 19.01        0.02       19.77

# save(FHNres_aos, file = "SA-ODEmorris_aos.RData")

# pdf("SA-ODEmorris_aos.pdf", width = 10, height = 7)
# Voltage:
plot(FHNres_aos, y_idx = 1, type = "sep")
plot(FHNres_aos, y_idx = 1, type = "trajec")
# Current:
plot(FHNres_aos, y_idx = 2, type = "sep")
plot(FHNres_aos, y_idx = 2, type = "trajec")
# dev.off()
