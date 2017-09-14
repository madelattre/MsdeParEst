# MsdeParEst R package ; file likelihoodNormal.r (last modified: 2017-08-28)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017

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

