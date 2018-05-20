Sensitivity Analysis of ODE Models
==============

[![CRAN status](https://www.r-pkg.org/badges/version/ODEsensitivity)](https://cran.r-project.org/package=ODEsensitivity)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/surmann/ODEsensitivity?branch=master&svg=true)](https://ci.appveyor.com/project/surmann/ODEsensitivity)
[![Travis build status](https://travis-ci.org/surmann/ODEsensitivity.svg?branch=master)](https://travis-ci.org/surmann/ODEsensitivity)
[![Coverage Status](https://coveralls.io/repos/github/surmann/ODEsensitivity/badge.svg?branch=master)](https://coveralls.io/github/surmann/ODEsensitivity?branch=master)

A package to perform sensitivity analysis for ordinary differential equation (ode)
models, e.g. the FitzHugh-Nagumo equations (cf. Ramsay et al., 2007) or low 
frequency oscillations (LFOs) in an electrical system (cf. Surmann et al., 
2014).
The package utilize the ode interface from `deSolve` and connects it with the 
sensitivity analysis from `sensitivity`. Additionally we add a method to
run the sensitivity analysis on variables with class `ODEnetwork`. A detailed
plotting function provides outputs on the calculations.
