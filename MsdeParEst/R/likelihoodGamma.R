# MsdeParEst R package ; file likelihoodGamma.r (last modified: 2017-08-09)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017
# UMR 518 AgroParisTech/INRA, 16 rue Claude Bernard, 75 231 Paris Cedex 05

#' Contrast based on the Euler approximation of the likelihood for parameter estimation when there is one random 
#' effect in the diffusion coefficient and no random effect in the drift.
#' 
#' @description Computation of the contrast based on the Euler approximation of -2 log-likelihood for the estimation
#' of the mixed SDE:
#' 
#'  \eqn{dXj(t)= (\alpha- \beta X_j(t))dt + \sigma_j a(X_j(t)) dW_j(t)}
#'  
#' with Gamma distribution for \eqn{1/\sigma_j^2} and fixed parameters in the drift.
#' @param a value of the shape of the Gamma distribution.
#' @param lambda value of the scape of the Gamma distribution.
#' @param U matrix of M sufficient statistics U (see \code{\link{UVS}}).
#' @param V list of the M sufficient statistics matrix V (see \code{\link{UVS}}).
#' @param S vector of the M sufficient statistics S (see \code{\link{UVS}}).
#' @param K number of times of observations.
#' @param drift.fixed values of the fixed effects in the drift. 
#' @return
#' \item{L}{value of the contrast}
#' @references
#' Estimation of population parameters in stochastic differential equations with random effects in the diffusion coefficient, M. Delattre, V. Genon-Catalot and A. Samson, \emph{ESAIM: Probability and Statistics 2015}, Vol 19, \bold{671 -- 688}
#' 
#' Parametric inference for discrete observations of diffusion processes with mixed effects, M. Delattre, V. Genon-Catalot and C. Laredo, \emph{Preprint 2016, hal-01332630}
#' 

contrastGamma <- function(a, lambda, U, V, S, K, drift.fixed) {
    
    invlambda <- 1/lambda
    
    M <- length(S)
    
    loglik <- vector(length = M)
    
    phi1 <- drift.fixed[1]
    phi2 <- drift.fixed[2]
    
    U1 <- rep(0, M)
    U2 <- rep(0, M)
    V11 <- rep(0, M)
    V22 <- rep(0, M)
    V12 <- rep(0, M)
    
    ST <- rep(0, M)
    
    for (j in 1:M) {
        U1[j] <- U[1, j]
        U2[j] <- U[2, j]
        V11[j] <- V[[j]][1, 1]
        V22[j] <- V[[j]][2, 2]
        V12[j] <- V[[j]][1, 2]
    }
    
    for (j in 1:M) {
        ST[j] <- (S[j] - 2 * (phi1 * U1[j] + phi2 * U2[j]) + phi1^2 * V11[j] + phi2^2 * 
            V22[j] + 2 * phi1 * phi2 * V12[j])
    }
    
    for (j in 1:M) {
        loglik[j] <- -2 * ((1 - a) * log(2/K + ST[j]/K) - log(gamma(a)) + a * log(invlambda) - 
            (a + K/2) * log((2 * invlambda/K + ST[j]/K)/(2/K + S[j]/K)))
    }
    
    L <- sum(loglik)
    
    return(L)
    
}

#' Computation of the Euler approximation of -2 log-likelihood when there is one random 
#' effect in the diffusion coefficient and fixed effects in the drift.
#' 
#' @description Computation of the Euler approximation of -2 log-likelihood of the mixed SDE:
#' 
#'  \eqn{dX_jj(t)= (\alpha - \beta X_j(t))dt + \sigma_j a(X_j(t)) dW_j(t)}
#'  
#' with inverse Gamma distribution for \eqn{\sigma_j^2} and fixed parameters \eqn{\alpha,\beta} in the drift.
#' @param a value of the shape of the Gamma distribution.
#' @param lambda value of the scape of the Gamma distribution.
#' @param U matrix of M sufficient statistics U (see \code{\link{UVS}}).
#' @param V list of the M sufficient statistics matrix V (see \code{\link{UVS}}).
#' @param S vector of the M sufficient statistics S (see \code{\link{UVS}}).
#' @param SigDelta vector of the M constant terms of the individual likelihood (see \code{\link{UVS}}).
#' @param K number of times of observations.
#' @param drift.fixed values of the fixed effects in the drift. 
#' @return
#' \item{L}{value of -2 x log-likelihood}
#' @references
#' Estimation of population parameters in stochastic differential equations with random effects in the diffusion coefficient, M. Delattre, V. Genon-Catalot and A. Samson, \emph{ESAIM: Probability and Statistics 2015}, Vol 19, \bold{671 -- 688}
#' 
#' Parametric inference for discrete observations of diffusion processes with mixed effects, M. Delattre, V. Genon-Catalot and C. Laredo, \emph{Preprint 2016}, hal-01332630
#' 

likelihoodGamma <- function(a, lambda, U, V, S, SigDelta, K, drift.fixed) {
    
    if (is.infinite(log(gamma(a + K/2)))) {
        L <- Inf
    } else {
        invlambda <- 1/lambda
        
        M <- length(S)
        
        loglik <- vector(length = M)
        
        phi1 <- drift.fixed[1]
        phi2 <- drift.fixed[2]
        
        U1 <- rep(0, M)
        U2 <- rep(0, M)
        V11 <- rep(0, M)
        V22 <- rep(0, M)
        V12 <- rep(0, M)
        
        ST <- rep(0, M)
        
        for (j in 1:M) {
            U1[j] <- U[1, j]
            U2[j] <- U[2, j]
            V11[j] <- V[[j]][1, 1]
            V22[j] <- V[[j]][2, 2]
            V12[j] <- V[[j]][1, 2]
        }
        
        for (j in 1:M) {
            ST[j] <- (S[j] - 2 * (phi1 * U1[j] + phi2 * U2[j]) + phi1^2 * V11[j] + phi2^2 * 
                V22[j] + 2 * phi1 * phi2 * V12[j])
        }
        
        for (j in 1:M) {
            loglik[j] <- K * log(2 * pi) + sum(SigDelta[j]) - 2 * (a * log(invlambda) + 
                log(gamma(a + K/2)) - log(gamma(a)) + (a + K/2) * log(invlambda + ST[j]/2))
        }
        
        
        L <- sum(loglik)
        
    }
    
    
    return(L)
    
}
