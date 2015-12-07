context("Test of ODEmorris_aos")

test_that("Type of result is correct", {
  # FHN-example:
  FHNmod <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      
      dVoltage <- s * (Voltage - Voltage^3 / 3 + Current)
      dCurrent <- - 1 / s *(Voltage - a + b * Current)
      
      return(list(c(dVoltage, dCurrent)))
    })
  }
  FHNpars  <- c(a = 0.2, b = 0.3, s = 3)
  FHNyini  <- c(Voltage = -1, Current = 1)
  FHNtimes1 <- seq(0.1, 101, by = 50)
  FHNres1 <- ODEmorris_aos(mod = FHNmod,
                       pars = names(FHNpars),
                       yini = FHNyini,
                       times = FHNtimes1,
                       seed = 2015,
                       binf = c(0.18, 0.18, 2.8),
                       bsup = c(0.22, 0.22, 3.2),
                       r = 4,
                       design =
                         list(type = "oat", levels = 100, grid.jump = 1),
                       scale = TRUE)
  
  expect_true(is.list(FHNres1))
  expect_equal(class(FHNres1), "morrisRes_aos")
  expect_equal(length(FHNres1), 2L)
  expect_equal(names(FHNres1), c("res", "pars"))
  expect_true(is.list(FHNres1$res))
  expect_equal(length(FHNres1$res), length(FHNyini))
  expect_equal(names(FHNres1$res), names(FHNyini))
  expect_true(is.matrix(FHNres1$res$Voltage))
  expect_true(is.matrix(FHNres1$res$Current))
  expect_equal(dim(FHNres1$res$Voltage), 
               c(1 + 3*length(FHNpars), length(FHNtimes1)))
  expect_equal(dim(FHNres1$res$Current), 
               c(1 + 3*length(FHNpars), length(FHNtimes1)))
  expect_true(is.vector(FHNres1$pars))
  expect_true(is.character(FHNres1$pars))
  expect_equal(length(FHNres1$pars), 3L)
  expect_equal(FHNres1$pars, names(FHNpars))
  
  # Only 1 point of time:
  FHNtimes2 <- 10
  FHNres2 <- ODEmorris_aos(mod = FHNmod,
                       pars = names(FHNpars),
                       yini = FHNyini,
                       times = FHNtimes2,
                       seed = 2015,
                       binf = c(0.18, 0.18, 2.8),
                       bsup = c(0.22, 0.22, 3.2),
                       r = 4,
                       design =
                         list(type = "oat", levels = 5, grid.jump = 1),
                       scale = TRUE)
  
  expect_true(is.list(FHNres2))
  expect_equal(class(FHNres2), "morrisRes_aos")
  expect_equal(length(FHNres2), 2L)
  expect_equal(names(FHNres2), c("res", "pars"))
  expect_true(is.list(FHNres2$res))
  expect_equal(length(FHNres2$res), length(FHNyini))
  expect_equal(names(FHNres2$res), names(FHNyini))
  expect_true(is.matrix(FHNres2$res$Voltage))
  expect_true(is.matrix(FHNres2$res$Current))
  expect_equal(dim(FHNres2$res$Voltage), 
               c(1 + 3*length(FHNpars), length(FHNtimes2)))
  expect_equal(dim(FHNres2$res$Current), 
               c(1 + 3*length(FHNpars), length(FHNtimes2)))
  expect_true(is.vector(FHNres2$pars))
  expect_true(is.character(FHNres2$pars))
  expect_equal(length(FHNres2$pars), 3L)
  expect_equal(FHNres2$pars, names(FHNpars))
})

test_that("Errors and warnings are correctly thrown", {
  # FHN-example:
  FHNmod <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      
      dVoltage <- s * (Voltage - Voltage^3 / 3 + Current)
      dCurrent <- - 1 / s *(Voltage - a + b * Current)
      
      return(list(c(dVoltage, dCurrent)))
    })
  }
  FHNpars  <- c(a = 0.2, b = 0.3, s = 3)
  FHNyini  <- c(Voltage = -1, Current = 1)
  FHNtimes2 <- 10
  
  expect_warning(ODEmorris_aos(mod = FHNmod,
                               pars = names(FHNpars),
                               yini = FHNyini,
                               times = FHNtimes2,
                               seed = 2015,
                               binf = c(0.18, 0.18, 2.8),
                               bsup = c(0.22, 0.22, 3.2),
                               r = 1,
                               design =
                                 list(type = "oat", levels = 5, grid.jump = 1),
                               scale = TRUE),
                 "Calculation of sigma requires r >= 2.")
  expect_error(ODEmorris_aos(mod = FHNmod,
                             pars = names(FHNpars),
                             yini = FHNyini,
                             times = FHNtimes2,
                             seed = 2015,
                             binf = c(0.18, 0.18, 2.8),
                             bsup = c(0.22, 0.22, 3.2),
                             r = 0,
                             design =
                               list(type = "oat", levels = 5, grid.jump = 1),
                             scale = TRUE),
               "r must be greater or equal to 1.")
})

