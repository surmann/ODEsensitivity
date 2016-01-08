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

# save(FHNres_aos, file = "test_ODEmorris_aos.RData")

# pdf("test_ODEmorris_aos.pdf", width = 10, height = 7)
# Palette "Dark2" from the package "RColorBrewer" with some 
# additional colors:
my_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
             "#E6AB02", "#A6761D", "#666666", "black", "firebrick",
             "darkblue", "darkgreen")
# Voltage:
plot(FHNres_aos, y_plot = "Voltage", type = "sep")
plot(FHNres_aos, y_plot = "Voltage", type = "sep", legendPos = "topleft")
plot(FHNres_aos, y_plot = "Voltage", type = "sep", colors_pars = my_cols)
plot(FHNres_aos, y_plot = "Voltage", type = "trajec")
plot(FHNres_aos, y_plot = "Voltage", type = "trajec", legendPos = "topleft")
plot(FHNres_aos, y_plot = "Voltage", type = "trajec", colors_pars = my_cols)
# Current:
plot(FHNres_aos, y_plot = "Current", type = "sep")
plot(FHNres_aos, y_plot = "Current", type = "sep", legendPos = "topleft")
plot(FHNres_aos, y_plot = "Current", type = "sep", colors_pars = my_cols)
plot(FHNres_aos, y_plot = "Current", type = "trajec")
plot(FHNres_aos, y_plot = "Current", type = "trajec", legendPos = "topleft")
plot(FHNres_aos, y_plot = "Current", type = "trajec", colors_pars = my_cols)
# dev.off()
