#' @title
#' Performing Sensitivity Analysis in ODE Models
#'
#' @description
#' \code{ODEsensitivity} provides methods to perform sensitivity analysis (SA)
#' in ordinary differential equation (ODE) models. Most of the functions are 
#' based on the implementations of Morris and Sobol' SA in the
#' \code{\link{sensitivity}} package (Pujol et al., 2015).
#'
#' @details
#' The main functions are \code{\link{ODEmorris}} and \code{\link{ODEsobol}},
#' which are built for objects of class \code{ODEnetwork}. The package 
#' \code{ODEnetwork} is required for these two functions to work. To perform
#' a sensitivity analysis for a more general type of ODEs, please use 
#' \code{\link{ODEmorris_trafo}, \link{ODEmorris_ats}, \link{ODEmorris_aos},
#' \link{ODEsobol_trafo}, \link{ODEsobol_ats}} or \code{\link{ODEsobol_aos}}.
#' However, these functions are included rather for historical reasons.
#' 
#' See the \code{\link{sensitivity}} package and its 
#' \code{\link[sensitivity]{morris}, \link[sensitivity]{soboljansen}} and 
#' \code{\link[sensitivity]{sobolmartinez}} implementations for further 
#' information.
#'
#' @docType package
#' @name ODEsensitivity
NULL
