#' Computation Of The Log Likelihood In Mixed Stochastic Differential Equations.
#' 
#' @description Computation of -2 loglikelihood of the mixed SDE with Normal distribution of the random effects
#' and fixed effect in the diffusion coefficient:
#'  \eqn{dXj(t)= (\alpha_j- \beta_j Xj(t))dt + \sigma a(Xj(t)) dWj(t)}.
#' @param mu current value of the mean of the normal distribution.
#' @param omega current value of the standard deviation of the normal distribution.
#' @param sigma current value of the diffusion coefficient.
#' @param U vector of the M sufficient statistics U (see \code{\link{UVS}}).
#' @param V vector of the M sufficient statistics V (see \code{\link{UVS}}).
#' @param S vector of the M sufficient statistics S (see \code{\link{UVS}}).
#' @param SigDelta vector of the M constant terms of the individual likelihood (see \code{\link{UVS}}). 
#' Required only if discrete = 1. Defaults to 0. 
#' @param K number of times of observations.
#' @param estimphi matrix of estimators of the random effects. 
#' @param drift.random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects.
#' @param discrete 1 for discrete observations, 0 otherwise. If discrete = 0, the exact likelihood associated with continuous observations is 
#' discretized. If discrete = 1, the likelihood of the Euler scheme of the mixed SDE is computed. Defaults to 1.  
#' @return
#' \item{L}{value of -2 x loglikelihood}
#' @references
#' Maximum likelihood estimation for stochastic differential equations with random effects, M. Delattre, V. Genon-Catalot and A. Samson, \emph{Scandinavian Journal of Statistics 2012}, Vol 40, \bold{322--343}
#' Parametric inference for discrete observations of diffusion processes with mixed effects, M. Delattre, V. Genon-Catalot and C. Larédo, \emph{Preprint}, hal-01332630
#' 

likelihoodNormal <- function(mu, omega, sigma, U, V, S, SigDelta = 0, K, estimphi, drift.random, 
    discrete = 1) {
    
    if (length(drift.random) == 2) {
        Omega <- matrix(c(omega[1]^2, 0, 0, omega[2]^2), 2, 2, byrow = TRUE)
    }
    
    if (sum(drift.random) == 1) {
        Omega <- matrix(c(omega[1]^2, 0, 0, 0), 2, 2, byrow = TRUE)
    }
    
    if (sum(drift.random) == 2) {
        Omega <- matrix(c(omega[1]^2, 0, 0, omega[2]^2), 2, 2, byrow = TRUE)
    }
    
    M <- dim(U)[2]
    
    loglik <- vector(length = M)
    
    I2 <- diag(c(1, 1))
    
    if (discrete == 1) {
        
        for (j in 1:M) {
            A <- (I2 + V[[j]] %*% Omega)
            Rinv <- solve(A) %*% V[[j]]
            b <- mu - estimphi[, j]
            loglik[j] <- 2 * K * log(sigma) + log(det(A)) + 1/sigma^2 * (t(b) %*% Rinv %*% 
                b - t(U[, j]) %*% estimphi[, j]) + S[j]/sigma^2 + K * log(2 * pi) + sum(SigDelta[j])
        }
        
    }
    
    if (discrete == 0) {
        
        for (j in 1:M) {
            A <- (I2 + V[[j]] %*% Omega)
            Rinv <- solve(A) %*% V[[j]]
            b <- mu - estimphi[, j]
            loglik[j] <- log(det(A)) + 1/sigma^2 * (t(b) %*% Rinv %*% b - t(U[, j]) %*% 
                estimphi[, j])
        }
    }
    
    L <- sum(loglik)
    
    return(L)
}

