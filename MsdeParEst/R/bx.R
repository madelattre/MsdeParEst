# MsdeParEst R package ; bx.r (last modified: 2017-08-09)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017
# UMR 518 AgroParisTech/INRA, 16 rue Claude Bernard, 75 231 Paris Cedex 05

#' Computation Of The Drift Coefficient
#' 
#' @description Computation of the drift coefficient
#' @param x vector of data
#' @return
#' \item{b}{The drift is \eqn{b(x,\phi) = \phi_1 b_1(x) + \phi_2 b_2(x)}, the output is the vector \eqn{(b_1,b_2)^t}}
#' @keywords drift



bx <- function(x) {
        return(matrix(c(1, -x), 2, 1))
}
