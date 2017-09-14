# MsdeParEst R package ; file UVS.r (last modified: 2017-08-28)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017

UVS <- function(X, model, times) {
    
    M <- dim(X)[1]
    K <- dim(X)[2]
    Xm <- X[, -K]
    delta <- diff(times)
    Tend <- times[K]
    
    S <- rep(NA, M)
    
    if (model == "OU") {
        for (j in 1:M) {
            S[j] <- sum((X[j, 2:K] - X[j, 1:(K - 1)])^2/delta)
        }
    }
    
    if (model == "CIR") {
        for (j in 1:M) {
            S[j] <- sum((X[j, 2:K] - X[j, 1:(K - 1)])^2/(delta * X[j, 1:(K - 1)]))
        }
    }
    
    
    U <- matrix(0, 2, M)
    
    V <- as.list(1:M)
    b <- as.list(1:M)
    
    Int1 <- rowSums(Xm * matrix(delta, M, length(delta)))  #Int1 <- apply(Xm * delta, 1, sum)
    
    if (model == "OU") {
        
        SigDelta <- rep(sum(log(delta)), M)
        
        Int2 <- rowSums(Xm^2 * matrix(delta, M, length(delta)))
        
        
        for (j in 1:M) {
            b[[j]] <- matrix(apply(matrix(X[j, ], 1, K), 2, bx), 2, K)  # 2xK  matrix
            
            U[, j] <- rowSums((b[[j]][, 1:(K - 1)] * matrix((X[j, 2:K] - X[j, 1:(K - 1)]), 
                2, K - 1, byrow = TRUE)))
            
            V[[j]] <- matrix(c(Tend, -Int1[j], -Int1[j], Int2[j]), 2, 2)
            
        }
        
    }
    
    if (model == "CIR") {
        
        SigDelta <- rep(0, M)
        
        Int3 <- rowSums(1/Xm * matrix(delta, M, length(delta)))
        
        b <- as.list(1:M)
        bsig <- as.list(1:M)
        
        for (j in 1:M) {
            b[[j]] <- matrix(apply(matrix(X[j, ], 1, K), 2, bx), 2, K)
            
            bsig[[j]] <- matrix(-1, 2, K)
            bsig[[j]][1, ] <- 1/X[j, ]
            
            U[, j] = rowSums((bsig[[j]][, 1:(K - 1)] * matrix((X[j, 2:K] - X[j, 1:(K - 
                1)]), 2, K - 1, byrow = TRUE)))
            
            V[[j]] <- matrix(c(Int3[j], -Tend, -Tend, Int1[j]), 2, 2)
            
            SigDelta[j] <- sum(log(delta) + log(X[j, 2:K]))
            
        }
    }
    
    
    return(list(U = U, V = V, S = S, SigDelta = SigDelta))
}
