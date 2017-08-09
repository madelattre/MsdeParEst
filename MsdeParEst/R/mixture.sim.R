# MsdeParEst R package ; file mixture.sim.r (last modified: 2017-08-09)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017
# UMR 518 AgroParisTech/INRA, 16 rue Claude Bernard, 75 231 Paris Cedex 05

#' Simulation Of A Mixture Of Normal Distributions
#' 
#' @description Simulation of M random variables from a mixture of Gaussian distributions
#' @param M number of simulated variables 
#' @param param vector of parameters with the means and standard-deviations of the normal distributions 
#' @param prob mixture components probabilities
#' @return
#' \item{Y}{vector of simulated variables}
#' @importFrom stats rnorm
#' @details
#' If the distribution is \eqn{p1 N(\mu1,\sigma1^2) + (1-p1)N(\mu2, \sigma2^2)}
#' 
#' param=c(\eqn{\mu1, \sigma1, \mu2, \sigma2}) and prob=c(p1,1-p1)
#' 




mixture.sim <- function(M, param, prob) {
    
    N <- length(prob)
    
    index <- sample(1:N, M, replace = TRUE, prob)
    
  
        if (dim(param)[2] == 4) { #two series of random variables are simulated
            Y <- matrix(NA, dim(param)[2]/2, M)
            for (j in 1:M) {
                for (n in 1:2) {
                  Y[n, j] <- rnorm(1, param[index[j], 2 * n - 1], param[index[j], 2 * n])
                }
            }
        }
        if (dim(param)[2] == 2) { #one series of random variables is simulated
            Y <- rep(NA, M)
            for (j in 1:M) {
                Y[j] <- rnorm(1, param[index[j],1], param[index[j],2])
            }
        }
    
    res = list(Y=Y,index=index)
    
    return(res)
}


