# Testing ODEmorris.ODEnetwork():

library(ODEnetwork)
masses <- c(1, 1)
dampers <- diag(c(1, 1))
springs <- diag(c(1, 1))
springs[1, 2] <- 1
distances <- diag(c(0, 2))
distances[1, 2] <- 1
lfonet <- ODEnetwork(masses, dampers, springs, 
                     cartesian = TRUE, distances = distances)
lfonet <- setState(lfonet, c(0.5, 1), c(0, 0))
LFOpars <- c("m.1", "d.1", "k.1", "k.1.2", "m.2", "d.2", "k.2")
LFOtimes <- seq(0.01, 20, by = 0.1)
LFObinf <- rep(0.001, length(LFOpars))
LFObsup <- c(2, 1.5, 6, 6, 2, 1.5, 6)

devtools::load_all(".")

# r = 50, ode_parallel = TRUE:
system.time(
  LFOres <- ODEmorris(lfonet, LFOpars, LFOtimes, 
                      seed = 2015, binf = LFObinf, bsup = LFObsup, r = 50, 
                      design = list(type = "oat", levels = 100, grid.jump = 1),
                      scale = TRUE, ode_method = "adams", 
                      ode_parallel = TRUE, ode_parallel_ncores = 2)
)
# User      System verstrichen 
# 1.46        0.13        9.30

# save(LFOres, file = "test_ODEmorris.ODEnetwork.Rdata")

# pdf("test_ODEmorris.ODEnetwork.pdf", width = 12, height = 10)
# Checking defaults:
plot(LFOres)
# Separate plots:
plot(LFOres, state_plot = "x.2", kind = "sep", type = "b")
plot(LFOres, state_plot = "x.2", kind = "sep", legendPos = "topleft")
plot(LFOres, state_plot = "x.2", kind = "sep", legendPos = "outside")
# Palette "Dark2" from the package "RColorBrewer" with some 
# additional colors:
my_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
             "#E6AB02", "#A6761D", "#666666", "black", "firebrick",
             "darkblue", "darkgreen")
plot(LFOres, state_plot = "x.2", kind = "sep", colors_pars = my_cols)

# Trajectories:
plot(LFOres, state_plot = "x.2", kind = "trajec", legendPos = "topleft")
plot(LFOres, state_plot = "x.2", kind = "trajec", legendPos = "outside")
plot(LFOres, state_plot = "x.2", kind = "trajec", colors_pars = my_cols)
# Testing the passing of arguments:
plot(LFOres, state_plot = "x.2", kind = "sep", colors_pars = my_cols, 
     main_title = "Big Title", legendPos = "outside", cex.axis = 2, 
     main = "Small Title", cex.main = 0.5)
# dev.off()
