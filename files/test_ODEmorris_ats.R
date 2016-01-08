# Testing ODEmorris_ats():

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
FHNres_ats <- ODEmorris_ats(mod = FHNmod,
                            pars = c("a", "b", "s"),
                            yini = FHNyini,
                            times = FHNtimes,
                            y_analyzed = "Voltage",
                            ode_method = "adams",
                            seed = 2015,
                            binf = c(0.18, 0.18, 2.8),
                            bsup = c(0.22, 0.22, 3.2),
                            r = 10,
                            design =
                              list(type = "oat", levels = 30, grid.jump = 1),
                            scale = TRUE)
)
# User      System verstrichen 
# 5.71        0.00        5.74

# save(FHNres_ats, file = "test_ODEmorris_ats.RData")

# pdf("test_ODEmorris_ats.pdf", width = 10, height = 7)
# Palette "Dark2" from the package "RColorBrewer" with some 
# additional colors:
my_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
             "#E6AB02", "#A6761D", "#666666", "black", "firebrick",
             "darkblue", "darkgreen")
plot(FHNres_ats, type = "sep")
plot(FHNres_ats, type = "sep", legendPos = "topleft")
plot(FHNres_ats, type = "sep", colors_pars = my_cols)
plot(FHNres_ats, type = "trajec")
plot(FHNres_ats, type = "trajec", legendPos = "topleft")
plot(FHNres_ats, type = "trajec", colors_pars = my_cols)
# dev.off()
