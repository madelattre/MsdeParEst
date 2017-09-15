# MsdeParEst R package ; file EstParamNormalGamma.r (last modified: 2017-09-15)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017

EstParamNormalGamma <- function(U, V, S, SigDelta, K, drift.random, drift.fixed = NULL) {
    
    M <- length(S)
    
    # k <- 0.1
    
    # estimGamma <- (K/S) * ((S/K) >= (k/sqrt(K)))
    
    estimGamma <- (K/S) * ((S/K) > 0)
    
    init.a <- mean(estimGamma)^2/var(estimGamma)
    init.lambda <- var(estimGamma)/mean(estimGamma)
    
    V1 <- function(param) {
        contrastGamma(exp(param[1]), exp(param[2]), U, V, S, K, c(0, 0))
    }
    
    res.gamma <- suppressWarnings(optim(c(log(init.a), log(init.lambda)), fn = V1, method = "Nelder-Mead"))
    a <- exp(res.gamma$par[1])
    lambda <- exp(res.gamma$par[2])
    
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

