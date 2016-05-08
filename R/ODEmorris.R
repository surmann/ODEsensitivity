#' @title Morris Screening for ODE Models
#'
#' @description
#'   \code{ODEmorris} is the generic function for performing a sensitivity 
#'   analysis of ODE models using Morris's elementary effects screening method.
#'
#' @param mod
#'   either a model function supplied in the manner as needed for 
#'   \code{\link[deSolve]{ode}} (for \code{\link{ODEmorris.default}}) or an 
#'   object of class \code{ODEnetwork} (for \code{\link{ODEmorris.ODEnetwork}}).
#' @param ...
#'   further arguments passed to methods, see \code{\link{ODEmorris.default}} 
#'   and \code{\link{ODEmorris.ODEnetwork}}.
#'
#' @details
#'   There are two methods for this generic function: 
#'   \code{\link{ODEmorris.default}} (for general ODE models) and
#'   \code{\link{ODEmorris.ODEnetwork}} (for objects of class \code{ODEnetwork},
#'   see package \code{ODEnetwork}).
#'
#' @author Frank Weber
#' @seealso \code{\link{ODEmorris.default}, \link{ODEmorris.ODEnetwork}}
#' 
#' @export
#'

ODEmorris <- function(mod, ...){
  UseMethod("ODEmorris", mod)
}
