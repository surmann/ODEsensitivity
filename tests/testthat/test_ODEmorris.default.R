context("Test of ODEmorris.default() and plot.morrisRes()")

# FHN-example:

FHNmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dVoltage <- s * (Voltage - Voltage^3 / 3 + Current)
    dCurrent <- - 1 / s *(Voltage - a + b * Current)
    
    return(list(c(dVoltage, dCurrent)))
  })
}
FHNstate  <- c(Voltage = -1, Current = 1)
FHNtimes1 <- seq(0.1, 20, by = 5)
FHNtimes2 <- 10

# Normal:
FHNres1 <- ODEmorris(mod = FHNmod,
                     pars = c("a", "b", "s"),
                     state_init = FHNstate,
                     times = FHNtimes1,
                     seed = 2015,
                     binf = c(0.18, 0.18, 2.8),
                     bsup = c(0.22, 0.22, 3.2),
                     r = 4,
                     design =
                       list(type = "oat", levels = 100, grid.jump = 1),
                     scale = TRUE,
                     ode_method = "adams",
                     ode_parallel = FALSE,
                     ode_parallel_ncores = NA)

# Only 1 point of time:
FHNres2 <- ODEmorris(mod = FHNmod,
                     pars = c("a", "b", "s"),
                     state_init = FHNstate,
                     times = FHNtimes2,
                     seed = 2015,
                     binf = c(0.18, 0.18, 2.8),
                     bsup = c(0.22, 0.22, 3.2),
                     r = 4,
                     design =
                       list(type = "oat", levels = 100, grid.jump = 1),
                     scale = TRUE,
                     ode_method = "adams",
                     ode_parallel = FALSE,
                     ode_parallel_ncores = NA)

# Only 1 parameter and 1 point of time:
FHNmod3 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dVoltage <- 3 * (Voltage - Voltage^3 / 3 + Current)
    dCurrent <- - 1 / 3 *(Voltage - a + 0.3 * Current)
    
    return(list(c(dVoltage, dCurrent)))
  })
}
FHNres3 <- ODEmorris(mod = FHNmod3,
                     pars = "a",
                     state_init = FHNstate,
                     times = FHNtimes2,
                     seed = 2015,
                     binf = 0.18,
                     bsup = 0.22,
                     r = 4,
                     design =
                       list(type = "oat", levels = 100, grid.jump = 1),
                     scale = TRUE,
                     ode_method = "adams",
                     ode_parallel = FALSE,
                     ode_parallel_ncores = NA)

# With parallelization:
FHNres_parallel <- ODEmorris(mod = FHNmod,
                             pars = c("a", "b", "s"),
                             state_init = FHNstate,
                             times = FHNtimes1,
                             seed = 2015,
                             binf = c(0.18, 0.18, 2.8),
                             bsup = c(0.22, 0.22, 3.2),
                             r = 4,
                             design = list(type = "oat", levels = 100, 
                                           grid.jump = 1),
                             scale = TRUE,
                             ode_method = "adams",
                             ode_parallel = TRUE,
                             ode_parallel_ncores = 2)

# Simplex design:
FHNres_simplex <- ODEmorris(mod = FHNmod,
                            pars = c("a", "b", "s"),
                            state_init = FHNstate,
                            times = FHNtimes1,
                            seed = 2015,
                            binf = c(0.18, 0.18, 2.8),
                            bsup = c(0.22, 0.22, 3.2),
                            r = 4,
                            design =
                              list(type = "simplex", scale.factor = 0.01),
                            scale = TRUE,
                            ode_method = "adams",
                            ode_parallel = FALSE,
                            ode_parallel_ncores = NA)

test_that("Type of ODEmorris.default()-result is correct", {
  # Normal case:
  expect_true(is.list(FHNres1))
  expect_equal(class(FHNres1), "morrisRes")
  expect_equal(length(FHNres1), length(FHNstate))
  expect_equal(names(FHNres1), names(FHNstate))
  expect_true(is.matrix(FHNres1$Voltage))
  expect_true(is.matrix(FHNres1$Current))
  expect_equal(dim(FHNres1$Voltage), 
               c(1 + 3*length(c("a", "b", "s")), length(FHNtimes1)))
  expect_equal(dim(FHNres1$Current), 
               c(1 + 3*length(c("a", "b", "s")), length(FHNtimes1)))
  
  # Only 1 point of time:
  expect_true(is.list(FHNres2))
  expect_equal(class(FHNres2), "morrisRes")
  expect_equal(length(FHNres2), length(FHNstate))
  expect_equal(names(FHNres2), names(FHNstate))
  expect_true(is.matrix(FHNres2$Voltage))
  expect_true(is.matrix(FHNres2$Current))
  expect_equal(dim(FHNres2$Voltage), 
               c(1 + 3*length(c("a", "b", "s")), length(FHNtimes2)))
  expect_equal(dim(FHNres2$Current), 
               c(1 + 3*length(c("a", "b", "s")), length(FHNtimes2)))
  
  # Only 1 parameter and 1 point of time:
  expect_true(is.list(FHNres3))
  expect_equal(class(FHNres3), "morrisRes")
  expect_equal(length(FHNres3), length(FHNstate))
  expect_equal(names(FHNres3), names(FHNstate))
  expect_true(is.matrix(FHNres3$Voltage))
  expect_true(is.matrix(FHNres3$Current))
  expect_equal(dim(FHNres3$Voltage), 
               c(1 + 3*length(c("a")), length(FHNtimes2)))
  expect_equal(dim(FHNres3$Current), 
               c(1 + 3*length(c("a")), length(FHNtimes2)))
  
  # With parallelization:
  expect_true(is.list(FHNres_parallel))
  expect_equal(class(FHNres_parallel), "morrisRes")
  expect_equal(length(FHNres_parallel), length(FHNstate))
  expect_equal(names(FHNres_parallel), names(FHNstate))
  expect_true(is.matrix(FHNres_parallel$Voltage))
  expect_true(is.matrix(FHNres_parallel$Current))
  expect_equal(dim(FHNres_parallel$Voltage), 
               c(1 + 3*length(c("a", "b", "s")), length(FHNtimes1)))
  expect_equal(dim(FHNres_parallel$Current), 
               c(1 + 3*length(c("a", "b", "s")), length(FHNtimes1)))
  
  # Simplex design:
  expect_true(is.list(FHNres_simplex))
  expect_equal(class(FHNres_simplex), "morrisRes")
  expect_equal(length(FHNres_simplex), length(FHNstate))
  expect_equal(names(FHNres_simplex), names(FHNstate))
  expect_true(is.matrix(FHNres_simplex$Voltage))
  expect_true(is.matrix(FHNres_simplex$Current))
  expect_equal(dim(FHNres_simplex$Voltage), 
               c(1 + 3*length(c("a", "b", "s")), length(FHNtimes1)))
  expect_equal(dim(FHNres_simplex$Current), 
               c(1 + 3*length(c("a", "b", "s")), length(FHNtimes1)))
})

test_that("ODEmorris.default() throws errors and warnings", {
  # bsup < binf:
  expect_warning(FHNres_binf_bsup <- 
                   ODEmorris(mod = FHNmod,
                             pars = c("a", "b", "s"),
                             state_init = FHNstate,
                             times = FHNtimes1,
                             seed = 2015,
                             binf = c(0.22, 0.18, 2.8),
                             bsup = c(0.18, 0.22, 3.2),
                             r = 4,
                             design =
                               list(type = "oat", levels = 100, grid.jump = 1),
                             scale = TRUE,
                             ode_method = "adams",
                             ode_parallel = FALSE,
                             ode_parallel_ncores = NA),
                 paste("At least one element of \"bsup\" was lower than the",
                       "corresponding element of \"binf\".", 
                       "Elements were swapped."))
  expect_equal(FHNres1, FHNres_binf_bsup)
  
  # r = 1:
  expect_warning(ODEmorris(mod = FHNmod,
                           pars = c("a", "b", "s"),
                           state_init = FHNstate,
                           times = FHNtimes2,
                           seed = 2015,
                           binf = c(0.18, 0.18, 2.8),
                           bsup = c(0.22, 0.22, 3.2),
                           r = 1,
                           design =
                             list(type = "oat", levels = 100, grid.jump = 1),
                           scale = TRUE,
                           ode_method = "adams",
                           ode_parallel = FALSE,
                           ode_parallel_ncores = NA),
                 "Calculation of sigma requires r >= 2.")
  
  # r = 0:
  expect_error(ODEmorris(mod = FHNmod,
                         pars = c("a", "b", "s"),
                         state_init = FHNstate,
                         times = FHNtimes2,
                         seed = 2015,
                         binf = c(0.18, 0.18, 2.8),
                         bsup = c(0.22, 0.22, 3.2),
                         r = 0,
                         design =
                           list(type = "oat", levels = 100, grid.jump = 1),
                         scale = TRUE,
                         ode_method = "adams",
                         ode_parallel = FALSE,
                         ode_parallel_ncores = NA),
               "\"r\" must be greater than or equal to 1.")
})

test_that("Plots are generated", {
  expect_true(plot(FHNres1))
  expect_true(plot(FHNres2))
  expect_true(plot(FHNres3))
  # Trajectories:
  expect_true(plot(FHNres1, kind = "trajec"))
  expect_true(plot(FHNres2, kind = "trajec"))
  expect_true(plot(FHNres3, kind = "trajec"))
  # Non-default arguments:
  expect_true(plot(FHNres1, state_plot = "Current", main_title = "Hi!", 
                   legendPos = "topleft", type = "b"))
  # Custom colors:
  my_cols <- c("firebrick", "chartreuse3", "dodgerblue")
  expect_true(plot(FHNres1, state_plot = "Current", colors_pars = my_cols))
  # Checking the passing of arguments:
  expect_true(plot(FHNres1, state_plot = "Current", cex.axis = 2, cex = 4, 
                   main = "Small Title", cex.main = 0.5))
})
