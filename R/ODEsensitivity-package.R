#' @title
#' Performing Sensitivity Analysis in ODE Models
#'
#' @description
#' \code{ODEsensitivity} provides methods to perform sensitivity analysis (SA)
#' in ordinary differential equation (ODE) models. Most of the functions are 
#' based on the implementations of Morris and Sobol' SA in the
#' \code{\link{sensitivity}} package (Pujol et al., 2015). However, a modified 
#' version of the \code{\link{sensitivity}}-package is needed for this package 
#' to work.
#'
#' @details
#' The main functions are \code{\link{ODEmorris}} and \code{\link{ODEsobol}},
#' which are generic functions and have a default method for general ODE models
#' (\code{\link{ODEmorris.default}, \link{ODEsobol.default}}) as well as a 
#' method for objects of class \code{ODEnetwork} 
#' (\code{\link{ODEmorris.ODEnetwork}, \link{ODEsobol.ODEnetwork}}). For the
#' latter two methods, the package \code{ODEnetwork} is required.
#' 
#' See the \code{\link{sensitivity}} package and its 
#' \code{\link[sensitivity]{morris}, \link[sensitivity]{soboljansen}} and 
#' \code{\link[sensitivity]{sobolmartinez}} implementations for further 
#' information on sensitivity analysis in \code{R}.
#'
#' @docType package
#' @name ODEsensitivity
NULL
