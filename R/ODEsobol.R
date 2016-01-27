#' @title Sobol' SA for ODE Models
#'
#' @description
#' \code{ODEsobol} is the generic function for Sobol' sensitivity analysis in
#' ODE models.
#'
#' @param mod 
#'   either a function with arguments \code{Time, State, Pars} (for 
#'   \code{\link{ODEsobol.default}} or an object of class \code{ODEnetwork} (for
#'   \code{\link{ODEsobol.ODEnetwork}}.
#' @param ...
#'   further arguments passed to methods, see \code{\link{ODEsobol.default}} and
#'   \code{\link{ODEsobol.ODEnetwork}}.
#'
#' @details
#'   There are two methods for this generic function: 
#'   \code{\link{ODEsobol.default}} (for general ODE models) and
#'   \code{\link{ODEsobol.ODEnetwork}} (for objects of class \code{ODEnetwork},
#'   see the package \code{ODEnetwork}).
#'
#' @author Frank Weber
#' @seealso \code{\link{ODEsobol.default}, \link{ODEsobol.ODEnetwork}}
#' 
#' @export
#'

ODEsobol <- function(mod, ...){
  UseMethod("ODEsobol", mod)
}
