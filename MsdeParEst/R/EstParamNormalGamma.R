# MsdeParEst R package ; file EstParamNormalGamma.r (last modified: 2017-08-09)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017
# UMR 518 AgroParisTech/INRA, 16 rue Claude Bernard, 75 231 Paris Cedex 05

#' Estimation In Mixed Stochastic Differential Equations with random effects in the drift and in the diffusion coefficient
#' 
#' @description Estimation of the parameters of the mixed SDE
#' 
#'  \eqn{dX_j(t)= (\alpha_j - \beta_j X_j(t))dt + \sigma_j a(X_j(t)) dW_j(t)}
#'  
#' with Normal conditional distribution of the random effects in the drift \eqn{\alpha_j,\beta_j | \sigma_j \sim N(\mu,\sigma_j^2 \Omega)}
#' and the square root of an inverse Gamma distributed random effect in the diffusion \eqn{\sigma_j^2 \sim Gamma(a,\lambda)}.
#' @param U matrix of M sufficient statistics U (see \code{\link{UVS}}).
#' @param V list of the M sufficient statistics matrix V (see \code{\link{UVS}}).
#' @param S vector of the M sufficient statistics S (see \code{\link{UVS}}).
#' @param SigDelta vector of the M constant terms of the individual likelihood (see \code{\link{UVS}}). 
#' @param K number of times of observations.
#' @param drift.random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects.
#' @param drift.fixed NULL if the fixed effect(s) in the drift is (are) estimated, value of the fixed effect(s) otherwise. Default to NULL.
#' @return
#' \item{mu}{estimated value of the mean of the Normal distribution}
#' \item{omega}{estimated value of the standard deviation of the Normal distribution}
#' \item{a}{estimated value of the shape of the Gamma distribution.}
#' \item{lambda}{estimated value of the scale of the Gamma distribution.}
#' \item{BIChere}{BIC indicator}
#' \item{AIChere}{AIC indicator}
#' @importFrom stats optim
#' @importFrom stats var
#' @references
#' Estimaton of the joint distribution of random effects for a discretely observed diffusion with random effects, M. Delattre, V. Genon-Catalot and C. Laredo, \emph{Preprint}, hal-01446063.


EstParamNormalGamma <- function(U, V, S, SigDelta, K, drift.random, drift.fixed = NULL) {
    
    M <- length(S)
    
    k <- 0.1
    
    estimGamma <- (K/S) * ((S/K) >= (k/sqrt(K)))
    
    init.a <- mean(estimGamma)^2/var(estimGamma)
    init.lambda <- var(estimGamma)/mean(estimGamma)
    
    V1 <- function(param) {
        contrastGamma(param[1], param[2], U, V, S, K, c(0, 0))
    }
    
    res.gamma <- suppressWarnings(optim(c(init.a, init.lambda), fn = V1, method = "Nelder-Mead"))
    a <- res.gamma$par[1]
    lambda <- res.gamma$par[2]
    
    estimphi <- matrix(NA, 2, M)
    
    for (j in 1:M) {
        estimphi[, j] <- solve(V[[j]]) %*% U[, j]
    }
    
    init.mu <- apply(estimphi, 1, mean)
    
    if (length(drift.random) == 2) {
        
        init.omega2 <- c(mean(estimGamma * (estimphi[1, ] - mean(estimphi[1, ]))^2), mean(estimGamma * 
            (estimphi[2, ] - mean(estimphi[2, ]))^2))
        
        V2 <- function(param) {
            contrastNormal(c(param[1], param[2]), c(param[3], param[4]), U, V, S, K, estimphi, 
                drift.random)
        }
        
        res.normal <- suppressWarnings(optim(c(init.mu, sqrt(init.omega2)), fn = V2, method = "Nelder-Mead"))
        mu <- c(res.normal$par[1], res.normal$par[2])
        omega <- c(abs(res.normal$par[3]), abs(res.normal$par[4]))
        
        BIChere <- likelihoodNormalGamma(a, lambda, mu, omega, U, V, S, SigDelta, K, drift.random) + 
            6 * log(M)
        AIChere <- likelihoodNormalGamma(a, lambda, mu, omega, U, V, S, SigDelta, K, drift.random) + 
            12
    }
    
    if (is.null(drift.fixed)) {
        
        if (sum(drift.random) == 1) {
            
            init.omega2 <- mean(estimGamma * (estimphi[1, ] - mean(estimphi[1, ]))^2)
            
            
            V2 <- function(param) {
                contrastNormal(c(param[1], param[2]), c(param[3], 0), U = U, V = V, S = S, 
                  K = K, estimphi = estimphi, drift.random = drift.random)
            }
            
            res.normal <- suppressWarnings(optim(c(init.mu, sqrt(init.omega2)), fn = V2, 
                method = "Nelder-Mead"))
            mu <- c(res.normal$par[1], res.normal$par[2])
            omega <- c(abs(res.normal$par[3]), 0)
            
            BIChere <- likelihoodNormalGamma(a, lambda, mu, omega, U, V, S, SigDelta, K, 
                drift.random) + 5 * log(M)
            AIChere <- likelihoodNormalGamma(a, lambda, mu, omega, U, V, S, SigDelta, K, 
                drift.random) + 10
            
        }
        
        if (sum(drift.random) == 2) {
            init.omega2 <- mean(estimGamma * (estimphi[2, ] - mean(estimphi[2, ]))^2)
            
            
            V2 <- function(param) {
                contrastNormal(mu = c(param[1], param[2]), omega = c(0, param[3]), U = U, V = V, 
                  S = S, K = K, estimphi = estimphi, drift.random = drift.random)
            }
            
            res.normal <- suppressWarnings(optim(c(init.mu, sqrt(init.omega2)), fn = V2, 
                method = "Nelder-Mead"))
            mu <- c(res.normal$par[1], res.normal$par[2])
            omega <- c(0, abs(res.normal$par[3]))
            
            BIChere <- likelihoodNormalGamma(a, lambda, mu, omega, U, V, S, SigDelta, K, 
                drift.random) + 5 * log(M)
            AIChere <- likelihoodNormalGamma(a, lambda, mu, omega, U, V, S, SigDelta, K, 
                drift.random) + 10
        }
        
    }
    
    if (!is.null(drift.fixed)) {
        
        if (sum(drift.random) == 1) {
            
            init.omega2 <- mean(estimGamma * (estimphi[1, ] - mean(estimphi[1, ]))^2)
            
            V2 <- function(param) {
                contrastNormal(c(param[1], drift.fixed), c(param[2], 0), U = U, V = V, S = S, 
                  K = K, estimphi = estimphi, drift.random = drift.random)
            }
            
            res.normal <- suppressWarnings(optim(c(init.mu[1], sqrt(init.omega2)), fn = V2, 
                method = "Nelder-Mead"))
            mu <- c(res.normal$par[1], drift.fixed)
            omega <- c(abs(res.normal$par[2]), 0)
            
            BIChere <- likelihoodNormalGamma(a, lambda, mu, omega, U, V, S, SigDelta, K, 
                drift.random) + 5 * log(M)
            AIChere <- likelihoodNormalGamma(a, lambda, mu, omega, U, V, S, SigDelta, K, 
                drift.random) + 10
            
        }
        
        if (sum(drift.random) == 2) {
            init.omega2 <- mean(estimGamma * (estimphi[2, ] - mean(estimphi[2, ]))^2)
            
            
            V2 <- function(param) {
                contrastNormal(mu = c(drift.fixed, param[1]), omega = c(0, param[2]), U = U, 
                  V = V, S = S, K = K, estimphi = estimphi, drift.random = drift.random)
            }
            
            res.normal <- suppressWarnings(optim(c(init.mu, sqrt(init.omega2)), fn = V2, 
                method = "Nelder-Mead"))
            mu <- c(drift.fixed, res.normal$par[1])
            omega <- c(0, abs(res.normal$par[2]))
            
            BIChere <- likelihoodNormalGamma(a, lambda, mu, omega, U, V, S, SigDelta, K, 
                drift.random) + 5 * log(M)
            AIChere <- likelihoodNormalGamma(a, lambda, mu, omega, U, V, S, SigDelta, K, 
                drift.random) + 10
        }
        
    }
    return(list(mu = mu, omega = omega, a = a, lambda = lambda, BIChere = BIChere, AIChere = AIChere))
    
}

