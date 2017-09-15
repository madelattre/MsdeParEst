# MsdeParEst R package ; file contrastNormal.r (last modified: 2017-09-15)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017

contrastNormal <- function(mu, omega, U, V, S, K, estimphi, drift.random) {
    
    if (length(drift.random) == 2) {
        Omega <- matrix(c(omega[1]^2, 0, 0, omega[2]^2), 2, 2, byrow = TRUE)
    }
    
    if (sum(drift.random) == 1) {
        Omega <- matrix(c(omega[1]^2, 0, 0, 0), 2, 2, byrow = TRUE)
    }
    
    if (sum(drift.random) == 2) {
        Omega <- matrix(c(0, 0, 0, omega[2]^2), 2, 2, byrow = TRUE)
    }
    
    M <- dim(U)[2]
    
    loglik <- vector(length = M)
    
    I2 <- diag(c(1, 1))
    
    for (j in 1:M) {
        A <- (I2 + V[[j]] %*% Omega)
        Rinv <- solve(A) %*% V[[j]]
        b <- mu - estimphi[, j]
        loglik[j] <- log(det(A))/2 + K/(2 * S[j]) * (t(b) %*% Rinv %*% b - t(U[, j]) %*% 
            estimphi[, j])
    }
    
    L <- sum(loglik)
    
    return(L)
}
