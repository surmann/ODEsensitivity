# Testing ODEmorris_aos():

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
FHNres_aos <- ODEmorris_aos(mod = FHNmod,
                            pars = c("a", "b", "s"),
                            yini = FHNyini,
                            times = FHNtimes,
                            ode_method = "lsoda",
                            seed = 2015,
                            binf = c(0.18, 0.18, 2.8),
                            bsup = c(0.22, 0.22, 3.2),
                            r = 10,
                            design =
                              list(type = "oat", levels = 30, grid.jump = 1),
                            scale = TRUE)
)
# User      System verstrichen 
# 9.18        0.02        9.51

# save(FHNres_aos, file = "SA-ODEmorris_aos.RData")

# pdf("SA-ODEmorris_aos.pdf", width = 10, height = 7)
# Voltage:
plot(FHNres_aos, y_idx = 1, type = "sep")
plot(FHNres_aos, y_idx = 1, type = "trajec")
# Current:
plot(FHNres_aos, y_idx = 2, type = "sep")
plot(FHNres_aos, y_idx = 2, type = "trajec")
# dev.off()
