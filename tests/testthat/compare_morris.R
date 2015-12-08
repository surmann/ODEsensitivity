context("Comparison of ODEmorris, ODEmorris_ats and ODEmorris_aos")

test_that("Results are equal", {
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
  FHNtimes <- seq(0.1, 101, by = 50)
  
  FHNres <- ODEmorris(mod = FHNmod,
                      pars = names(FHNpars),
                      yini = FHNyini,
                      times = FHNtimes,
                      seed = 2015,
                      binf = c(0.18, 0.18, 2.8),
                      bsup = c(0.22, 0.22, 3.2),
                      r = 4,
                      design =
                        list(type = "oat", levels = 100, grid.jump = 1),
                      scale = TRUE,
                      trafo = function(Y) Y[, 1],    # voltage only
                      ncores = 2)
  
  FHNres_ats <- ODEmorris_ats(mod = FHNmod,
                              pars = names(FHNpars),
                              yini = FHNyini,
                              times = FHNtimes,
                              y_idx = 1,
                              seed = 2015,
                              binf = c(0.18, 0.18, 2.8),
                              bsup = c(0.22, 0.22, 3.2),
                              r = 4,
                              design =
                                list(type = "oat", levels = 100, grid.jump = 1),
                              scale = TRUE)
  
  FHNres_aos <- ODEmorris_aos(mod = FHNmod,
                              pars = names(FHNpars),
                              yini = FHNyini,
                              times = FHNtimes,
                              seed = 2015,
                              binf = c(0.18, 0.18, 2.8),
                              bsup = c(0.22, 0.22, 3.2),
                              r = 4,
                              design =
                                list(type = "oat", levels = 100, grid.jump = 1),
                              scale = TRUE)
  
  stopifnot(all.equal(FHNres$res, FHNres_ats$res,
                      tolerance = 1e-4, check.attributes = FALSE),
            all.equal(FHNres$res, FHNres_aos$res$Voltage,
                      tolerance = 1e-4, check.attributes = FALSE),
            all.equal(FHNres_ats$res, FHNres_aos$res$Voltage,
                      tolerance = 1e-16))
  
  FHNres_current <- ODEmorris(mod = FHNmod,
                              pars = names(FHNpars),
                              yini = FHNyini,
                              times = FHNtimes,
                              seed = 2015,
                              binf = c(0.18, 0.18, 2.8),
                              bsup = c(0.22, 0.22, 3.2),
                              r = 4,
                              design =
                                list(type = "oat", levels = 100, grid.jump = 1),
                              scale = TRUE,
                              trafo = function(Y) Y[, 2],    # current only
                              ncores = 2)
  
  FHNres_ats_current <- ODEmorris_ats(mod = FHNmod,
                              pars = names(FHNpars),
                              yini = FHNyini,
                              times = FHNtimes,
                              y_idx = 2,
                              seed = 2015,
                              binf = c(0.18, 0.18, 2.8),
                              bsup = c(0.22, 0.22, 3.2),
                              r = 4,
                              design =
                                list(type = "oat", levels = 100, grid.jump = 1),
                              scale = TRUE)
  
  stopifnot(all.equal(FHNres_current$res, FHNres_ats_current$res,
                      tolerance = 1e-4, check.attributes = FALSE),
            all.equal(FHNres_current$res, FHNres_aos$res$Current,
                      tolerance = 1e-4, check.attributes = FALSE),
            all.equal(FHNres_ats_current$res, FHNres_aos$res$Current,
                      tolerance = 1e-16))
})
