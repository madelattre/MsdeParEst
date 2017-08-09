# MsdeParEst R package ; file likelihoodNormalGamma.r (last modified: 2017-08-09)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017
# UMR 518 AgroParisTech/INRA, 16 rue Claude Bernard, 75 231 Paris Cedex 05

#' Computation Of The Log Likelihood of the Euer scheme in Mixed Stochastic Differential Equations with random effects in the drift and in the diffusion.
#' 
#' @description Computation of -2 log-likelihood of the Euler scheme the mixed SDE 
#' 
#'  \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j a(X_j(t)) dW_j(t)},
#'  
#' with Normal conditional distribution of the random effects in the drift \eqn{\alpha_j,\beta_j|\sigma_j \sim N(\mu,\sigma_j^2\Omega)}
#' and inverse Gamma distribution of the random effect in the diffusion \eqn{1/\sigma_j^2 \sim Gamma(a,\lambda)}.
#' @param a current value of the shape of the Gamma distribution.
#' @param lambda current value of the scape of the Gamma distribution.
#' @param mu current value of the mean of the normal distribution.
#' @param omega current value of the standard deviation of the normal distribution.
#' @param U vector of the M sufficient statistics U (see \code{\link{UVS}}).
#' @param V list of the M sufficient statistics V (see \code{\link{UVS}}).
#' @param S vector of the M sufficient statistics S (see \code{\link{UVS}}).
#' @param SigDelta vector of the M constant terms of the individual likelihood (see \code{\link{UVS}}). 
#' Required only if discrete = 1. Defaults to 0. 
#' @param K number of times of observations.
#' @param drift.random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects.
#' @return
#' \item{L}{value of -2 x loglikelihood}
#' @references
#' Estimation of the joint distribution of random effects for a discretely observed diffusion with random effects, M. Delattre, V. Genon-Catalot and C. Laredo, \emph{Preprint}, hal-01446063.


likelihoodNormalGamma <- function(a, lambda, mu, omega, U, V, S, SigDelta, K, drift.random) {
    
    if (is.infinite(suppressWarnings(log(gamma(a + K/2))))) {
        L <- Inf
    } else {
        invlambda <- 1/lambda
        ## Ajout
        I2 <- diag(c(1,1))
        ## Ajout
        M <- length(S)
        
        ## Ajout
        estimphi <- matrix(NA, 2, M)
        
        for (j in 1:M) {
          estimphi[, j] <- solve(V[[j]]) %*% U[, j]
        }
        ## Ajout
        
        loglik <- vector(length = M)
        LL <- vector(length = M)
        #V2 <- vector(length = M)
        
        ST <- vector(length = M)
        
        ## Ajout
        if (length(drift.random) == 2) {
          Omega <- matrix(c(omega[1]^2, 0, 0, omega[2]^2), 2, 2, byrow = TRUE)
        }
        
        if (sum(drift.random) == 1) {
          Omega <- matrix(c(omega[1]^2, 0, 0, 0), 2, 2, byrow = TRUE)
        }
        
        if (sum(drift.random) == 2) {
          Omega <- matrix(c(0, 0, 0, omega[2]^2), 2, 2, byrow = TRUE)
        }
        ## Ajout
        
        for (j in 1:M) {
            A <- (I2 + V[[j]] %*% Omega)
            Rinv <- solve(A) %*% V[[j]]
            b <- mu - estimphi[, j]
            ST[j] <- S[j] + t(b) %*% Rinv %*% b - t(U[, j]) %*% estimphi[, j]
            
            LL[j] <- a * log(invlambda) + log(gamma(a + K/2)) - log(gamma(a)) - (a + K/2) * 
                log(invlambda + S[j]/2) -log(det(A))/2

            loglik[j] <- K * log(2 * pi) + sum(SigDelta[j]) - 2 * LL[j]
            
        }
        
        
        L <- sum(loglik)
        
    }
    
    
    return(L)
    
}
