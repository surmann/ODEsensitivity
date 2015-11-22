library(ODEnetwork)

masses <- c(1, 1)
dampers <- diag(c(1, 1))
springs <- diag(c(1, 1))
springs[1, 2] <- 1
distances <- diag(c(0, 2))
distances[1, 2] <- 1
times <- seq(0, 20, by = 0.01)

# analytisch
odenet <- ODEnetwork(masses, dampers, springs, distances=distances)
ODEmod <- createOscillators(odenet)
ODEpars  <- createParamVec(odenet)
ODEyini  <- createState(odenet)
ODEtimes <- seq(0.01, 40, by = 5) # seq(0.01, 40, by = 0.01)

devtools::load_all()

ODEres <- ODEmorris(mod = ODEmod,
                    pars = names(ODEpars),
                    yini = ODEyini,
                    times = ODEtimes,
                    seed = 2015,
                    binf = 0,
                    bsup = 1,
                    r = 10,
                    design =
                      list(type = "oat", levels = 50, grid.jump = 1),
                    trafo = function(Y) Y[, 1],    # position only
                    ncores = 2)

# odenet <- setState(odenet, c(0.5, 1), c(0, 0))
# odenet <- simuNetwork(odenet, times)
# plot(odenet, state = "1")

# # nummerisch
# eventdata <- data.frame(var = c("v.1"), time = c(0), value = c(0))
# odenet <- setEvents(odenet, eventdata, type = "dirac")
# odenet <- simuNetwork(odenet, times)
# plot(odenet, state = "1")


