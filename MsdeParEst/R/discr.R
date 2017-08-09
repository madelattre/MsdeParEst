# MsdeParEst R package ; file discr.r (last modified: 2017-08-09)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017
# UMR 518 AgroParisTech/INRA, 16 rue Claude Bernard, 75 231 Paris Cedex 05

#' Simulation Of Random Variables
#' 
#' @description Simulation of (discrete) random variables from a vector of probability (the nonparametrically estimated values
#'  of the density renormalised to sum at 1) and a vector of real values (the grid of estimation)
#' @param x n real numbers
#' @param p vector of probability, length n
#' @return
#' y a simulated value from the discrete distribution
#' @importFrom stats runif



discr <- function(x, p) {
    cp <- cumsum(p)
    U <- runif(1, 0, max(cp))
    x[which(cp >= U)[1]]
}
