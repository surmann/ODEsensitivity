# Testing ODEmorris.default():

devtools::load_all()

FHNmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dVoltage <- s * (Voltage - Voltage^3 / 3 + Current)
    dCurrent <- - 1 / s *(Voltage - a + b * Current)
    
    return(list(c(dVoltage, dCurrent)))
  })
}
FHNstate  <- c(Voltage = -1, Current = 1)
FHNtimes <- seq(0.1, 50, by = 5)

# r = 10, ode_parallel = FALSE:
system.time(
  FHNres10 <- ODEmorris(mod = FHNmod,
                        pars = c("a", "b", "s"),
                        state_init = FHNstate,
                        times = FHNtimes,
                        seed = 2015,
                        binf = c(0.18, 0.18, 2.8),
                        bsup = c(0.22, 0.22, 3.2),
                        r = 10,
                        design =
                          list(type = "oat", levels = 100, grid.jump = 1),
                        scale = TRUE,
                        ode_method = "adams",
                        ode_parallel = FALSE,
                        ode_parallel_ncores = NA)
)
# User      System verstrichen 
# 5.89        0.00        5.90

# r = 10, ode_parallel = FALSE, simplex design:
system.time(
  FHNres10_simplex <- ODEmorris(mod = FHNmod,
                                pars = c("a", "b", "s"),
                                state_init = FHNstate,
                                times = FHNtimes,
                                seed = 2015,
                                binf = c(0.18, 0.18, 2.8),
                                bsup = c(0.22, 0.22, 3.2),
                                r = 10,
                                design = 
                                  list(type = "simplex", scale.factor = 0.01),
                                scale = TRUE,
                                ode_method = "adams",
                                ode_parallel = FALSE,
                                ode_parallel_ncores = NA)
)
# User      System verstrichen 
# 5.84        0.00        5.91

# r = 50, ode_parallel = FALSE:
system.time(
  FHNres50_nopar <- ODEmorris(mod = FHNmod,
                              pars = c("a", "b", "s"),
                              state_init = FHNstate,
                              times = FHNtimes,
                              seed = 2015,
                              binf = c(0.18, 0.18, 2.8),
                              bsup = c(0.22, 0.22, 3.2),
                              r = 50,
                              design = list(type = "oat", levels = 100, 
                                            grid.jump = 1),
                              scale = TRUE,
                              ode_method = "adams",
                              ode_parallel = FALSE,
                              ode_parallel_ncores = NA)
)
#  User      System verstrichen 
# 29.32        0.01       29.58

# r = 50, ode_parallel = FALSE, simplex design:
system.time(
  FHNres50_nopar_simplex <- ODEmorris(mod = FHNmod,
                                      pars = c("a", "b", "s"),
                                      state_init = FHNstate,
                                      times = FHNtimes,
                                      seed = 2015,
                                      binf = c(0.18, 0.18, 2.8),
                                      bsup = c(0.22, 0.22, 3.2),
                                      r = 50,
                                      design = list(type = "simplex", 
                                                    scale.factor = 0.01),
                                      scale = TRUE,
                                      ode_method = "adams",
                                      ode_parallel = FALSE,
                                      ode_parallel_ncores = NA)
)
#  User      System verstrichen 
# 29.30        0.00       29.47

# r = 50, ode_parallel = TRUE:
system.time(
  FHNres <- ODEmorris(mod = FHNmod,
                      pars = c("a", "b", "s"),
                      state_init = FHNstate,
                      times = FHNtimes,
                      seed = 2015,
                      binf = c(0.18, 0.18, 2.8),
                      bsup = c(0.22, 0.22, 3.2),
                      r = 50,
                      design = list(type = "oat", levels = 100, 
                                    grid.jump = 1),
                      scale = TRUE,
                      ode_method = "adams",
                      ode_parallel = TRUE,
                      ode_parallel_ncores = 2)
)
# User      System verstrichen 
# 0.20        0.11       18.13

# r = 50, ode_parallel = TRUE, simplex design:
system.time(
  FHNres_simplex <- ODEmorris(mod = FHNmod,
                              pars = c("a", "b", "s"),
                              state_init = FHNstate,
                              times = FHNtimes,
                              seed = 2015,
                              binf = c(0.18, 0.18, 2.8),
                              bsup = c(0.22, 0.22, 3.2),
                              r = 50,
                              design = 
                                list(type = "simplex", scale.factor = 0.01),
                              scale = TRUE,
                              ode_method = "adams",
                              ode_parallel = TRUE,
                              ode_parallel_ncores = 2)
)
# User      System verstrichen 
# 0.22        0.07       18.65

# save(FHNres, file = "test_ODEmorris.default.RData")

# pdf("test_ODEmorris.default.pdf", width = 10, height = 7)
# Checking defaults:
plot(FHNres)
# Palette "Dark2" from the package "RColorBrewer" with some 
# additional colors:
my_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
             "#E6AB02", "#A6761D", "#666666", "black", "firebrick",
             "darkblue", "darkgreen")
# Voltage:
plot(FHNres, state_plot = "Voltage", kind = "sep", type = "b")
plot(FHNres, state_plot = "Voltage", kind = "sep", legendPos = "topleft")
plot(FHNres, state_plot = "Voltage", kind = "sep", colors_pars = my_cols)
plot(FHNres, state_plot = "Voltage", kind = "trajec")
plot(FHNres, state_plot = "Voltage", kind = "trajec", legendPos = "topleft")
plot(FHNres, state_plot = "Voltage", kind = "trajec", colors_pars = my_cols)
# Current:
plot(FHNres, state_plot = "Current", kind = "sep")
plot(FHNres, state_plot = "Current", kind = "sep", legendPos = "topleft")
plot(FHNres, state_plot = "Current", kind = "sep", colors_pars = my_cols)
plot(FHNres, state_plot = "Current", kind = "trajec")
plot(FHNres, state_plot = "Current", kind = "trajec", type = "b")
plot(FHNres, state_plot = "Current", kind = "trajec", legendPos = "topleft")
plot(FHNres, state_plot = "Current", kind = "trajec", colors_pars = my_cols)
# dev.off()

# Checking other variations:
plot(FHNres10)
plot(FHNres10_simplex)
plot(FHNres50_nopar)
plot(FHNres50_nopar_simplex)
plot(FHNres)
plot(FHNres_simplex)
