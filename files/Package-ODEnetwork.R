##--------------------------------------------------------------------##
## FOR1511: R-Code                                                    ##
##                                                                    ##
##   ODEnetwork Package                                               ##
##--------------------------------------------------------------------##

options(help_type = "html")
# library("devtools")                            # R developing
# devtools::install_github("surmann/ODEnetwork")
library("ODEnetwork")                            # Schwingungen


##----------------------------------------------------------------------
## Standardbeispiel
##----------------------------------------------------------------------

masses <- c(1, 2)
dampers <- diag(c(0.1, 0.5))
dampers[1, 2] <- 0.05
springs <- diag(c(4, 10))
springs[1, 2] <- 6
odenet <- ODEnetwork(masses, dampers, springs)
# Start definieren (Auslenkung(1, 2), Geschwindigkeit (1, 2)):
odenet <- setState(odenet, c(1, 3), c(0, 0))
# Simulation fuer Zeitpunkte starten:
odenet <- simuNetwork(odenet, seq(0, 60, by = 0.05))
plot(odenet)
plot(odenet, var = 2L)
plot(odenet, state = "1")
plot(odenet, state = "2")
plot(odenet, state = "1vs2")

