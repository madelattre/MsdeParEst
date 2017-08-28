# MsdeParEst R package ; file likelihoodGamma.r (last modified: 2017-08-28)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017
# UMR 518 AgroParisTech/INRA, 16 rue Claude Bernard, 75 231 Paris Cedex 05

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
