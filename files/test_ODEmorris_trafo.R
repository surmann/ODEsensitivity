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

# save(FHNres_trafo, file = "test_ODEmorris_trafo.RData")

# pdf("test_ODEmorris_trafo.pdf", width = 10, height = 7)
# Palette "Dark2" from the package "RColorBrewer" with some 
# additional colors:
my_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
             "#E6AB02", "#A6761D", "#666666", "black", "firebrick",
             "darkblue", "darkgreen")
plot(FHNres_trafo, type = "sep")
plot(FHNres_trafo, type = "sep", legendPos = "topleft")
plot(FHNres_trafo, type = "sep", colors_pars = my_cols)
plot(FHNres_trafo, type = "trajec")
plot(FHNres_trafo, type = "trajec", legendPos = "topleft")
plot(FHNres_trafo, type = "trajec", colors_pars = my_cols)
# dev.off()
