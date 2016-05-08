#' @title Sobol' Sensitivity Analysis for ODE Models
#'
#' @description
#' \code{ODEsobol} is the generic function for performing a Sobol' sensitivity 
#' analysis of ODE models.
#'
#' @param mod
#'   either a model function supplied in the manner as needed for 
#'   \code{\link[deSolve]{ode}} (for \code{\link{ODEsobol.default}}) or an 
#'   object of class \code{ODEnetwork} (for \code{\link{ODEsobol.ODEnetwork}}).
#' @param ...
#'   further arguments passed to methods, see \code{\link{ODEsobol.default}} and
#'   \code{\link{ODEsobol.ODEnetwork}}.
#'
#' @details
#'   There are two methods for this generic function: 
#'   \code{\link{ODEsobol.default}} (for general ODE models) and
#'   \code{\link{ODEsobol.ODEnetwork}} (for objects of class \code{ODEnetwork},
#'   see package \code{ODEnetwork}).
#'
#' @author Frank Weber
#' @seealso \code{\link{ODEsobol.default}, \link{ODEsobol.ODEnetwork}}
#' 
#' @export
#'

ODEsobol <- function(mod, ...){
  UseMethod("ODEsobol", mod)
}
