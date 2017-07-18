#' Estimation In Mixed Stochastic Differential Equations with random effects in the drift and fixed effect in the diffusion coefficient
#' 
#' @description Estimation of the parameters of the mixed SDE with Normal distribution of the random effects in the drift and
#' fixed parameter in the diffusion:
#'  \eqn{dXj(t)= (\alpha_j - \beta_j Xj(t))dt + \sigma a(Xj(t)) dWj(t)}, done with \code{\link{likelihoodNormal}}.
#' @param U matrix of M sufficient statistics U (see \code{\link{UVS}}).
#' @param V list of the M sufficient statistics matrix V (see \code{\link{UVS}}).
#' @param S vector of the M sufficient statistics S (see \code{\link{UVS}}).
#' @param SigDelta vector of the M constant terms of the individual likelihood (see \code{\link{UVS}}). 
#' Required only if discrete = 1. Defaults to 0. 
#' @param K number of times of observations.
#' @param drift.fixed value of the fixed effect in the drift if it is not estimated. Default to 0.
#' @param estim.drift.fix 1 if the fixed effect in the drift is estimated, 0 otherwise. Default to 1.
#' @param sigma value of the fixed effect in the diffusion if known (not estimated). Defaults to 0.
#' @param drift.random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects.
#' @param diffusion.estim 1 if sigma is estimated, 0 otherwise.
#' @param discrete 1 for discrete observations, 0 otherwise. If discrete = 0, the exact likelihood associated with continuous observations is 
#' discretized. If discrete = 1, the likelihood of the Euler scheme of the mixed SDE is computed. Defaults to 1.
#' @return
#' \item{mu}{estimated value of the mean of the Normal distribution}
#' \item{omega}{estimated value of the standard deviation of the Normal distribution}
#' \item{sigma}{value of the diffusion coefficient}
#' \item{BIChere}{BIC indicator}
#' \item{AIChere}{AIC indicator}

EstParamNormal <- function(U, V, S, SigDelta = 0, K, drift.fixed = 0, estim.drift.fix = 1, 
    sigma = 0, drift.random, diffusion.estim = 1, discrete = 1) {
    
    
    M <- length(S)
    
    estimphi <- matrix(NA, 2, M)
    
    for (j in 1:M) {
        estimphi[, j] <- solve(V[[j]]) %*% U[, j]
    }
    
    init.mu <- apply(estimphi, 1, mean)
    
    if (discrete == 1) {
        
        if (diffusion.estim == 1) {
            
            k <- 0.1
            
            estimGamma <- (K/S) * ((S/K) >= (k/sqrt(K)))
            
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
                
                
                if (estim.drift.fix == 1) {
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
                
                if (estim.drift.fix == 0) {
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
                
                
                if (estim.drift.fix == 1) {
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
                
                if (estim.drift.fix == 0) {
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
        
        if (diffusion.estim == 0) {
            
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
                
                if (estim.drift.fix == 1) {
                  ln = function(param) {
                    likelihoodNormal(c(param[1], param[2]), c(param[3], 0), sigma, U, V, 
                      S, SigDelta, K, estimphi, drift.random, discrete)
                  }
                  
                  res <- optim(c(init.mu, sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
                  
                  mu <- c(res$par[1], res$par[2])
                  omega <- c(abs(res$par[3]), 0) * sigma
                  
                  nbparam <- 3
                }
                
                if (estim.drift.fix == 0) {
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
                
                if (estim.drift.fix == 1) {
                  ln = function(param) {
                    likelihoodNormal(c(param[1], param[2]), c(0, param[3]), sigma, U, V, 
                      S, SigDelta, K, estimphi, drift.random, discrete)
                  }
                  
                  res <- optim(c(init.mu, sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
                  
                  mu <- c(res$par[1], res$par[2])
                  omega <- c(0, abs(res$par[3])) * sigma
                  
                  nbparam <- 3
                }
                
                if (estim.drift.fix == 0) {
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
        
        if (diffusion.estim == 1) {
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
            
            nbparam <- 4 + diffusion.estim
        }
        
        if (sum(drift.random) == 1) {
            init.omega2 <- var(estimphi[1, ])/sigma^2
            
            if (estim.drift.fix == 1) {
                ln = function(param) {
                  likelihoodNormal(c(param[1], param[2]), c(param[3], 0), sigma, U, V, 
                    S, SigDelta, K, estimphi, drift.random, discrete)
                }
                
                res <- optim(c(init.mu, sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
                
                mu <- c(res$par[1], res$par[2])
                omega <- c(abs(res$par[3]), 0) * sigma
                
                nbparam <- 3 + diffusion.estim
            }
            
            if (estim.drift.fix == 0) {
                ln = function(param) {
                  likelihoodNormal(c(param[1], drift.fixed), c(param[2], 0), sigma, U, 
                    V, S, SigDelta, K, estimphi, drift.random, discrete)
                }
                
                res <- optim(c(init.mu[1], sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
                
                mu <- c(res$par[1], drift.fixed)
                omega <- c(abs(res$par[2]), 0) * sigma
                
                nbparam <- 2 + diffusion.estim
            }
            
        }
        
        if (sum(drift.random) == 2) {
            init.omega2 <- var(estimphi[2, ])/sigma^2
            
            if (estim.drift.fix == 1) {
                ln = function(param) {
                  likelihoodNormal(c(param[1], param[2]), c(0, param[3]), sigma, U, V, 
                    S, SigDelta, K, estimphi, drift.random, discrete)
                }
                
                res <- optim(c(init.mu, sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
                
                mu <- c(res$par[1], res$par[2])
                omega <- c(0, abs(res$par[3])) * sigma
                
                nbparam <- 3 + diffusion.estim
            }
            
            if (estim.drift.fix == 0) {
                ln = function(param) {
                  likelihoodNormal(c(drift.fixed, param[1]), c(0, param[2]), sigma, U, 
                    V, S, SigDelta, K, estimphi, drift.random, discrete)
                }
                
                res <- optim(c(init.mu[2], sqrt(init.omega2)), fn = ln, method = "Nelder-Mead")
                
                mu <- c(drift.fixed, res$par[1])
                omega <- c(0, abs(res$par[2])) * sigma
                
                nbparam <- 2 + diffusion.estim
            }
            
        }
    }
    
    
    
    BIChere <- likelihoodNormal(mu, omega, sigma, U, V, S, SigDelta, K, estimphi, drift.random, 
        discrete) + nbparam * log(M)
    AIChere <- likelihoodNormal(mu, omega, sigma, U, V, S, SigDelta, K, estimphi, drift.random, 
        discrete) + nbparam * 2
    
    return(list(mu = mu, omega = omega, sigma = sigma, BIChere = BIChere, AIChere = AIChere))
    
}
