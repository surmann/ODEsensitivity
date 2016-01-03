# Make sure that ODEmorris_trafo() still works as before:

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
FHNres_trafo <- ODEmorris_trafo(mod = FHNmod,
                    pars = c("a", "b", "s"),
                    yini = FHNyini,
                    times = FHNtimes,
                    ode_method = "adams",
                    seed = 2015,
                    binf = c(0.18, 0.18, 2.8),
                    bsup = c(0.22, 0.22, 3.2),
                    r = 10,
                    design =
                      list(type = "oat", levels = 30, grid.jump = 1),
                    scale = TRUE,
                    trafo = function(Y) Y[, 1],    # voltage only
                    ncores = 2)
)
# User      System verstrichen 
# 0.17        0.05       25.58

# save(FHNres_trafo, file = "SA-ODEmorris.RData")

# pdf("SA-ODEmorris_trafo.pdf", width = 10, height = 7)
plot(FHNres_trafo, type = "sep")
plot(FHNres_trafo, type = "trajec")
# dev.off()
