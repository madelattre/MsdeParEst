#' Simulation Of A Mixture Of Two Normal Or Gamma Distributions
#' 
#' @description Simulation of M random variables from a mixture of two Gaussian or Gamma distributions
#' @param M number of simulated variables 
#' @param density.phi name of the chosen density 'mixture.normal' or 'mixture.gamma'
#' @param param vector of parameters with the proportion of mixture of the two distributions and means and standard-deviations of the two normal or 
#'   shapes and scales of the two Gamma distribution 
#' @param prob mixture components probabilities
#' @return
#' \item{Y}{vector of simulated variables}
#' @importFrom stats rnorm
#' @details
#' If 'mixture.normal', the distribution is \eqn{p N(\mu1,\sigma1^2) + (1-p)N(\mu2, \sigma2^2)}
#' 
#' and param=c(p, \eqn{\mu1, \sigma1, \mu2, \sigma2})
#' 
#' If 'mixture.gamma', the distribution is \eqn{p Gamma(shape1,scale1) + (1-p)Gamma(shape2,scale2)}
#' 
#' and param=c(p, shape1, scale1, shape2, scale2)





mixture.sim <- function(M, density.phi, param, prob) {
    
    N <- length(prob)
    
    index <- sample(1:N, M, replace = TRUE, prob)
    
    if (density.phi == "mixture.normal") {
        if (is.matrix(param) == 1) {
            Y <- matrix(NA, 2, M)
            for (j in 1:M) {
                for (n in 1:2) {
                  Y[n, j] <- rnorm(1, param[n, 2 * index[j] - 1], param[n, 2 * index[j]])
                }
            }
        }
        if (is.matrix(param) == 0) {
            Y <- rep(NA, M)
            for (j in 1:M) {
                Y[j] <- rnorm(1, param[2 * index[j] - 1], param[2 * index[j]])
            }
        }
        
    }
    
    # if (density.phi == "mixture.gamma") {
    #     for (j in 1:M) {
    #         Y[j] <- rgamma(1, shape = param[2 * index[j] - 1], scale = param[2 * index[j]])
    #     }
    # }
    
    return(Y)
}
