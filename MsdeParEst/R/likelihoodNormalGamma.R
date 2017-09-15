# MsdeParEst R package ; file likelihoodNormalGamma.r (last modified: 2017-09-15)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017

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
