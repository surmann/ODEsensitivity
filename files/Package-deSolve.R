##--------------------------------------------------------------------##
## FOR1511: R-Code                                                    ##
##                                                                    ##
##   deSolve Package                                                  ##
##--------------------------------------------------------------------##

options(help_type = "html")
library("ODEnetwork")                          # Schwingungen
library("checkmate")                           # Plausibilitaet
library("deSolve")                             # numerische Lsg. DGLs


##----------------------------------------------------------------------
## Bsp. aus Paket
##----------------------------------------------------------------------

# Example1: Predator-Prey Lotka-Volterra model (with logistic prey)
LVmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Ingestion    <- rIng  * Prey * Predator
    GrowthPrey   <- rGrow * Prey * (1 - Prey/K)
    MortPredator <- rMort * Predator

    dPrey        <- GrowthPrey - Ingestion
    dPredator    <- Ingestion * assEff - MortPredator

    return(list(c(dPrey, dPredator)))
  })
}

LVpars  <- c(rIng   = 0.2,    # /day, rate of ingestion (= Nahrungsaufnahme)
             rGrow  = 1.0,    # /day, growth rate of prey
             rMort  = 0.2 ,   # /day, mortality rate of predator
             assEff = 0.5,    # -, assimilation efficiency
             K      = 10)     # mmol/m3, carrying capacity

LVyini  <- c(Prey = 1, Predator = 2)
LVtimes <- seq(0, 200, by = 1)
LVout   <- ode(LVyini, LVtimes, LVmod, LVpars)
summary(LVout)

## Default plot method
plot(LVout)

## User specified plotting
matplot(LVout[ , 1], LVout[ , 2:3], type = "l", xlab = "time",
        ylab = "Conc", main = "Lotka-Volterra", lwd = 2)
legend("topright", c("prey", "predator"), col = 1:2, lty = 1:2)


##----------------------------------------------------------------------
## Bsp. DGL einfach
##----------------------------------------------------------------------

# DGL definieren:
myMod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {

    dX <- alpha * X

    return(list(c(dX)))
  })
}

mypars  <- c(alpha = 2)     # some parameter alpha

myyini  <- c(X = 1)
mytimes <- seq(0, 10, by = 1)
myout   <- ode(myyini, mytimes, myMod, mypars)
summary(myout)

## Default plot method
plot(myout)    # genial!


##----------------------------------------------------------------------
## Bsp. aus ODEnetwork
##----------------------------------------------------------------------

masses    <- 1
dampers   <- as.matrix(0.5)
springs   <- as.matrix(4)
distances <- as.matrix(2)

odenet <- ODEnetwork(masses, dampers, springs, distances=distances)

# only state
odenet <- setState(odenet, 3, 0)
odenet <- simuNetwork(odenet, seq(0, 10, by = 0.1))
plot(odenet)
plot(odenet, state = "1vs2")

# Intera in "simuNetwork" [OS for oscillator]:
OSmod   <- ODEnetwork:::createOscillators(odenet)
OSpars  <- ODEnetwork:::createParamVec(odenet)
OSyini  <- ODEnetwork:::createState(odenet)
OStimes <- seq(0, 10, by = 0.01)

OSout <- ode(OSyini, OStimes, OSmod, OSpars)

# plot:
plot(OSout)


##----------------------------------------------------------------------
## Bsp. FitzHugh-Nagumo
##----------------------------------------------------------------------

FHNmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {

    dVoltage <- s * (Voltage - Voltage^3 / 3 + Current)
    dCurrent <- - 1 / s *(Voltage - a + b * Current)

    return(list(c(dVoltage, dCurrent)))
  })
}

FHNpars  <- c(a = 0.2,     # Paramter a
              b = 0.3,     # Paramter b
              s = 3)       # Paramter s (= c in der Originalnotation)

FHNyini  <- c(Voltage = -1, Current = 1)
FHNtimes <- seq(0, 20, by = 0.01)
FHNout   <- ode(FHNyini, FHNtimes, FHNmod, FHNpars)
head(FHNout)
summary(FHNout)

## Default plot method
plot(FHNout)  # sieht exakt so aus wie in dem Paper.

