# Testing ODEmorris() for objects of class "ODEnetwork":

masses <- c(1, 1)
dampers <- diag(c(1, 1))
springs <- diag(c(1, 1))
springs[1, 2] <- 1
distances <- diag(c(0, 2))
distances[1, 2] <- 1
odenet <- ODEnetwork(masses, dampers, springs, 
                     cartesian = TRUE, distances = distances)
odenet <- setState(odenet, c(0.5, 1), c(0, 0))

ODEtimes <- seq(0.01, 20, by = 0.1)
ODEbinf <- c(rep(0.001, 9), -4, 0.001, -10)
ODEbsup <- c(2, 1.5, 6, 4, 6, 10, 2, 1.5, 6, -0.001, 6, -0.001)

devtools::load_all(".")

system.time(
  ODEres <- ODEmorris(odenet, ODEtimes, ode_method = "adams", seed = 2015, 
                      binf = ODEbinf, bsup = ODEbsup, r = 20)
)
# User      System verstrichen 
# 6.90        0.03        7.17

# pdf("test_ODEmorris.pdf", width = 12, height = 10)
plot(ODEres, y_idx = 1, type = "sep", legendPos = "topleft")
# dev.off()
