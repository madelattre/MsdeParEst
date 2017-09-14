# MsdeParEst R package ; file EstParamNormal.r (last modified: 2017-08-28)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017

EstParamNormal <- function(U, V, S, SigDelta = 0, K, drift.fixed = NULL, 
    sigma = NULL, drift.random, discrete = 1) {
    
    
    M <- length(S)
    
    estimphi <- matrix(NA, 2, M)
    
    for (j in 1:M) {
        estimphi[, j] <- solve(V[[j]]) %*% U[, j]
    }
    
    init.mu <- apply(estimphi, 1, mean)
    
    if (discrete == 1) {
        
        if (is.null(sigma)) {
            
            #k <- 0.1
            
            # estimGamma <- (K/S) * ((S/K) >= (k/sqrt(K)))
            
            estimGamma <- (K/S) * ((S/K) > 0)
          
            # if (prod(estimGamma==0)==1){
            #   estimGamma <- (K/S) * ((S/K) > 0)
            # }
            
            init.gamma <- mean(estimGamma)
            init.sigma <- 1/sqrt(init.gamma)
            
            
            if (length(drift.random) == 2) {
                init.omega2 <- apply(estimphi, 1, var) * init.gamma
                
                
                ln = function(param) {
                  likelihoodNormal(c(param[1], param[2]), c(param[3], param[4]), exp(param[5]), 
                    U, V, S, SigDelta, K, estimphi, drift.random, discrete)
                }
                
                res <- optim(c(init.mu, sqrt(init.omega2), init.sigma), fn = ln, method = "Nelder-Mead")
                
                mu <- c(res$par[1], res$par[2])
                sigma <- exp(res$par[5])
                omega <- c(abs(res$par[3]), abs(res$par[4])) * sigma
                
                nbparam <- 5
            }
            
            if (sum(drift.random) == 1) {
                init.omega2 <- var(estimphi[1, ]) * init.gamma
                
                
                if (is.null(drift.fixed)) {
                  ln = function(param) {
                    likelihoodNormal(c(param[1], param[2]), c(param[3], 0), exp(param[4]), 
                      U, V, S, SigDelta, K, estimphi, drift.random, discrete)
                  }
                  
                  res <- optim(c(init.mu, sqrt(init.omega2), init.sigma), fn = ln, method = "Nelder-Mead")
                  
                  mu <- c(res$par[1], res$par[2])
                  sigma <- exp(res$par[4])
                  omega <- c(abs(res$par[3]), 0) * sigma
                  
                  nbparam <- 4
                }
                
                if (!is.null(drift.fixed)) {
                  ln = function(param) {
                    likelihoodNormal(c(param[1], drift.fixed), c(param[2], 0), exp(param[3]), 
                      U, V, S, SigDelta, K, estimphi, drift.random, discrete)
                  }
                  
                  res <- optim(c(init.mu[1], sqrt(init.omega2), init.sigma), fn = ln, method = "Nelder-Mead")
                  
                  mu <- c(res$par[1], drift.fixed)
                  sigma <- exp(res$par[3])
                  omega <- c(abs(res$par[2]), 0) * sigma
                  
                  nbparam <- 3
                }
                
            }
            
            if (sum(drift.random) == 2) {
                init.omega2 <- var(estimphi[2, ]) * init.gamma
                
                
                if (is.null(drift.fixed)) {
                  ln = function(param) {
                    likelihoodNormal(c(param[1], param[2]), c(0, param[3]), exp(param[4]), 
                      U, V, S, SigDelta, K, estimphi, drift.random, discrete)
                  }
                  
                  res <- optim(c(init.mu, sqrt(init.omega2), init.sigma), fn = ln, method = "Nelder-Mead")
                  
                  mu <- c(res$par[1], res$par[2])
                  sigma <- exp(res$par[4])
                  omega <- c(0, abs(res$par[3])) * sigma
                  
                  nbparam <- 4
                }
                
                if (!is.null(drift.fixed)) {
                  ln = function(param) {
                    likelihoodNormal(c(drift.fixed, param[1]), c(0, param[2]), exp(param[3]), 
                      U, V, S, SigDelta, K, estimphi, drift.random, discrete)
                  }
                  
                  res <- optim(c(init.mu[2], sqrt(init.omega2), init.sigma), fn = ln, method = "Nelder-Mead")
                  
                  mu <- c(drift.fixed, res$par[1])
                  sigma <- exp(res$par[3])
                  omega <- c(0, abs(res$par[2])) * sigma
                  
                  nbparam <- 3
                }
                
            }
            
        }
        
        if (!is.null(sigma)) {
            
            init.mu <- apply(estimphi, 1, mean)
            
            if (length(drift.random) == 2) {
                
                init.omega2 <- apply(estimphi, 1, var)/sigma^2
                
                
                ln = function(param) {
                  likelihoodNormal(c(param[1], param[2]), c(param[3], param[4]), sigma, 
                    U, V, S, SigDelta, K, estimphi, drift.random, discrete)
                }
                
                res <- optim(c(init.mu, sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
                
                mu <- c(res$par[1], res$par[2])
                omega <- c(abs(res$par[3]), abs(res$par[4])) * sigma
                
                nbparam <- 4
            }
            
            if (sum(drift.random) == 1) {
                init.omega2 <- var(estimphi[1, ])/sigma^2
                
                if (is.null(drift.fixed)) {
                  ln = function(param) {
                    likelihoodNormal(c(param[1], param[2]), c(param[3], 0), sigma, U, V, 
                      S, SigDelta, K, estimphi, drift.random, discrete)
                  }
                  
                  res <- optim(c(init.mu, sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
                  
                  mu <- c(res$par[1], res$par[2])
                  omega <- c(abs(res$par[3]), 0) * sigma
                  
                  nbparam <- 3
                }
                
                if (!is.null(drift.fixed)) {
                  ln = function(param) {
                    likelihoodNormal(c(param[1], drift.fixed), c(param[2], 0), sigma, U, 
                      V, S, SigDelta, K, estimphi, drift.random, discrete)
                  }
                  
                  res <- optim(c(init.mu[1], sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
                  
                  mu <- c(res$par[1], drift.fixed)
                  omega <- c(abs(res$par[2]), 0) * sigma
                  
                  nbparam <- 2
                }
                
            }
            
            if (sum(drift.random) == 2) {
                init.omega2 <- var(estimphi[2, ])/sigma^2
                
                if (is.null(drift.fixed)) {
                  ln = function(param) {
                    likelihoodNormal(c(param[1], param[2]), c(0, param[3]), sigma, U, V, 
                      S, SigDelta, K, estimphi, drift.random, discrete)
                  }
                  
                  res <- optim(c(init.mu, sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
                  
                  mu <- c(res$par[1], res$par[2])
                  omega <- c(0, abs(res$par[3])) * sigma
                  
                  nbparam <- 3
                }
                
                if (!is.null(drift.fixed)) {
                  ln = function(param) {
                    likelihoodNormal(c(drift.fixed, param[1]), c(0, param[2]), sigma, U, 
                      V, S, SigDelta, K, estimphi, drift.random, discrete)
                  }
                  
                  res <- optim(c(init.mu[2], sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
                  
                  mu <- c(drift.fixed, res$par[1])
                  omega <- c(0, abs(res$par[2])) * sigma
                  
                  nbparam <- 2
                }
                
            }
        }
    }
    
    
    if (discrete == 0) {
        
        if (is.null(sigma)) {
            sigma <- sqrt(mean(S/(K - 1)))
        }
        
        if (length(drift.random) == 2) {
            init.omega2 <- apply(estimphi, 1, var)/sigma^2
            
            
            ln = function(param) {
                likelihoodNormal(c(param[1], param[2]), c(param[3], param[4]), sigma, U, 
                  V, S, SigDelta, K, estimphi, drift.random, discrete)
            }
            
            res <- optim(c(init.mu, sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
            
            mu <- c(res$par[1], res$par[2])
            omega <- c(abs(res$par[3]), abs(res$par[4])) * sigma
            
            nbparam <- 4 + is.null(sigma)
        }
        
        if (sum(drift.random) == 1) {
            init.omega2 <- var(estimphi[1, ])/sigma^2
            
            if (is.null(drift.fixed)) {
                ln = function(param) {
                  likelihoodNormal(c(param[1], param[2]), c(param[3], 0), sigma, U, V, 
                    S, SigDelta, K, estimphi, drift.random, discrete)
                }
                
                res <- optim(c(init.mu, sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
                
                mu <- c(res$par[1], res$par[2])
                omega <- c(abs(res$par[3]), 0) * sigma
                
                nbparam <- 3 + is.null(sigma)
            }
            
            if (!is.null(drift.fixed)) {
                ln = function(param) {
                  likelihoodNormal(c(param[1], drift.fixed), c(param[2], 0), sigma, U, 
                    V, S, SigDelta, K, estimphi, drift.random, discrete)
                }
                
                res <- optim(c(init.mu[1], sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
                
                mu <- c(res$par[1], drift.fixed)
                omega <- c(abs(res$par[2]), 0) * sigma
                
                nbparam <- 2 + is.null(sigma)
            }
            
        }
        
        if (sum(drift.random) == 2) {
            init.omega2 <- var(estimphi[2, ])/sigma^2
            
            if (is.null(drift.fixed)) {
                ln = function(param) {
                  likelihoodNormal(c(param[1], param[2]), c(0, param[3]), sigma, U, V, 
                    S, SigDelta, K, estimphi, drift.random, discrete)
                }
                
                res <- optim(c(init.mu, sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
                
                mu <- c(res$par[1], res$par[2])
                omega <- c(0, abs(res$par[3])) * sigma
                
                nbparam <- 3 + is.null(sigma)
            }
            
            if (!is.null(drift.fixed)) {
                ln = function(param) {
                  likelihoodNormal(c(drift.fixed, param[1]), c(0, param[2]), sigma, U, 
                    V, S, SigDelta, K, estimphi, drift.random, discrete)
                }
                
                res <- optim(c(init.mu[2], sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
                
                mu <- c(drift.fixed, res$par[1])
                omega <- c(0, abs(res$par[2])) * sigma
                
                nbparam <- 2 + is.null(sigma)
            }
            
        }
    }
    
    
    
    BIChere <- likelihoodNormal(mu, omega, sigma, U, V, S, SigDelta, K, estimphi, drift.random, 
        discrete) + nbparam * log(M)
    AIChere <- likelihoodNormal(mu, omega, sigma, U, V, S, SigDelta, K, estimphi, drift.random, 
        discrete) + nbparam * 2
    
    return(list(mu = mu, omega = omega, sigma = sigma, BIChere = BIChere, AIChere = AIChere))
    
}
