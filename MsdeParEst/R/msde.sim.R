#' Simulation Of A Mixed Stochastic Differential Equation
#' 
#' @description Simulation of M independent trajectories of a mixed stochastic differential equation (SDE) with linear drift
#'  \deqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j a(X_j(t)) dW_j(t)  ,  j=1, ..., M.}
#' There may be two random effects \eqn{(\alpha_j, \beta_j)} in the drift and one random effect \eqn{\sigma_j} in the diffusion coefficient.
#' @param M number of trajectories.
#' @param T horizon of simulation.
#' @param N number of simulation steps, default Tx100. 
#' @param model name of the SDE: 'OU' (Ornstein-Uhlenbeck) or 'CIR' (Cox-Ingersoll-Ross).
#' @param drift.random random effects in the drift: 0 if no random effect, 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects.
#' @param diffusion.random random effect in the diffusion coefficient: 0 if no random effect, 1 if one multiplicative random effect.
#' @param mixture 1 if the random effects in the drift follow a mixture of Normal distributions, 0 otherwise. Default to 0.
#' @param drift.param fixed effects in the drift: value of the fixed effect when there is only one random effect, 0 otherwise. 
#' If drift.random =2, fixed can be 0 but \eqn{\beta} has to be a non negative random variable for the estimation.
#' vector (not mixture) or matrix (mixture) of parameters of the distribution of the random effects in the drift.
#' @param diffusion.param diffusion parameter if the diffusion coefficient is fixed, vector of parameters of the distribution of the diffusion random effect otherwise. 
#' @param nb.mixt number of mixture components if the drift random effects follow a mixture distribution, default nb.mixt=1.
#' @param mixt.prop vector of mixture proportions if the drift random effects follow a mixture distribution, default mixt.prop=1.
#' @param t0 time origin, default 0.
#' @param X0 initial value of the process, default X0=0.001.
#' @param delta time step of the simulation (T/N).
#' @param op.plot 1 if a plot of the trajectories is required, default 0. 
#' @param add.plot 1 for add trajectories to an existing plot
#' @return
#' \item{X}{matrix (M x (N+1)) of the M trajectories. }
#' \item{phi}{vector (or matrix) of the M simulated random effects.}
#' 
#' @importFrom sde sde.sim
#' 
#' @details
#' Simulation of M independent trajectories of the SDE (the Brownian motions \eqn{W_j} are independent), with linear drift. There may be one or two
#' random effects in the drift: 
#' 
#' If drift.random = 0, \eqn{\alpha} and \eqn{\beta} are fixed effects (\eqn{\alpha_j \equiv \alpha} and \eqn{\beta_j \equiv \beta})
#' 
#' If drift.random = 1,  \eqn{\beta} is a fixed effect (\eqn{\beta_j \equiv \beta}), and the drift function is written \eqn{(\alpha_j- \beta X_j(t))}
#' 
#' If drift.random = 2, \eqn{\alpha} is a fixed effect (\eqn{\alpha_j \equiv \alpha}), and the drift function is written \eqn{(\alpha- \beta_j X_j(t))}
#' 
#' If drift.random = c(1,2), both effects are random, and the drift function is written \eqn{(\alpha_j- \beta_j X_j(t))}
#' 
#' Two diffusions are implemented:
#' 
#' Ornstein-Uhlenbeck model (OU): \eqn{a(X_j(t))=1}
#' 
#' Cox-Ingersoll-Ross model (CIR): \eqn{a(X_j(t))=\sqrt{X_j(t)}}
#' 
#' There may be either a fixed or a random effect in the diffusion coefficient:
#' 
#' If diffusion.random = 0, \eqn{\sigma} is a fixed effect (\eqn{\sigma_j \equiv \sigma}). In that case, the random effects in the drift follow a Normal distribution 
#' or a mixture of Normal distributions.
#' 
#' If diffusion.random = 1, \eqn{\sigma_j} is a random effect. In that case, \eqn{\sigma_j^2} follow an inverse Gamma distribution with parameters diffusion.param=c(shape,scale),
#' and conditional on \eqn{\sigma_j}, the random effects in the drift follow a Normal distribution with mean mu and variance \eqn{Omega*\sigma_j^2}
#'
#' @references This function mixedsde.sim is based on the package sde, function sde.sim. See Simulation and Inference for stochastic differential equation, S.Iacus, \emph{Springer Series in Statistics 2008}
#' Chapter 2
#' @seealso \url{http://cran.r-project.org/package=sde}
#' 



msde.sim <- function(M, T, N = 100, model, drift.random, diffusion.random, mixture = 0, 
    drift.param, diffusion.param, nb.mixt = 1, mixt.prop = 1, t0 = 0, 
    X0 = 0.01, delta = T/N, op.plot = 0, add.plot = FALSE) {
    
    
    ## local sde.sim to sink undesired output away into a tempfile
    con <- file(tempfile(), open = "w")
    on.exit(close(con))
    sde.sim <- function(...) {
        sink(con)
        res <- sde::sde.sim(...)
        sink(NULL)
        res
    }
    
    if (missing(X0)) {
        message("Be careful, X0 is missing thus the initial value X0=0.01 is used")
    }
    
    if (((sum(drift.random)) == 0) && (diffusion.random == 0)) {
      stop("There should be at least one random effect either in the drift or the diffusion coefficient.")
    }
    
    if (((diffusion.random == 1) && (mixture == 1))) {
      stop("If there is one random effect in the diffusion coefficient, the random effects in the drift can't follow a mixture of Normal distributions. Try mixture = 0.")
    }
    
    delta <- T/N
    times <- seq(t0, T, length = N + 1)
    
    X <- matrix(0, M, N + 1)
    
    index <- NA
    
    if (diffusion.random == 0) {
        
        sig <- diffusion.param[1]
        psi <- NA
        
        if ((sum(drift.random) > 2)) {
            # simulation of the random effects
            phi <- matrix(0, 2, M)
            
            if (mixture == 0) {
                phi[1, ] <- rnorm(M, drift.param[1], drift.param[2])
                phi[2, ] <- rnorm(M, drift.param[3], drift.param[4])
            }
            if (mixture == 1) {
                phi[1, ] <- mixture.sim(M, drift.param[,c(1,2)], mixt.prop)
                phi[2, ] <- mixture.sim(M, drift.param[,c(3,4)], mixt.prop)
            }
            
            
            # simulation of the series
            if (model == "OU") {
                for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(phi[, j], sig), model = "OU"))
                }
            }
            
            if (model == "CIR") {
                for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "milstein", theta = c(phi[, j], sig), sigma.x = expression(sig/(2 * 
                        sqrt(x))), sigma = expression(sig * sqrt(x)), model = "CIR"))
                }
                
            }
        }
        
        if (sum(drift.random) == 1) {
            
            # simulation of the random effects
            phi <- rep(0, M)
            if (mixture == 1) {
                sim <- mixture.sim(M, drift.param[,c(1,2)], mixt.prop)
                phi <- sim$Y
                index <- sim$index
                
                # simulation of the series
                if (model == "OU") {
                  for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                                                       method = "EA", theta = c(phi[j], drift.param[index[j],3], sig), model = "OU"))
                  }
                }
                
                if (model == "CIR") {
                  for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, N = N, X0 = X0, delta = delta, 
                                                       method = "milstein", theta = c(phi[j], drift.param[index[j],3], sig), sigma.x = expression(sig/(2 * 
                                                                                                                                                sqrt(x))), sigma = expression(sig * sqrt(x)), model = "CIR"))
                  }
                }
            }
            
            
            if (mixture == 0) {
                phi <- drift.param[1] + drift.param[2] * rnorm(M, mean = 0, sd = 1)
                
                # simulation of the series
                if (model == "OU") {
                  for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                                                       method = "EA", theta = c(phi[j], drift.param[3], sig), model = "OU"))
                  }
                }
                
                if (model == "CIR") {
                  for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, N = N, X0 = X0, delta = delta, 
                                                       method = "milstein", theta = c(phi[j], drift.param[3], sig), sigma.x = expression(sig/(2 * 
                                                                                                                                                sqrt(x))), sigma = expression(sig * sqrt(x)), model = "CIR"))
                  }
                }
            }
            
            
            

            
        }
        
        
        if (sum(drift.random) == 2) {
            # simulation of the random effects
            phi <- rep(0, M)
            if (mixture == 1) {
                sim <- mixture.sim(M, drift.param[,c(2,3)], mixt.prop)
                phi <- sim$Y
                index <- sim$index
                
                # simulation of the series
                if (model == "OU") {
                  for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                                                       method = "EA", theta = c(drift.param[index[j],1], phi[j], sig), model = "OU"))
                  }
                }
                
                if (model == "CIR") {
                  for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(t0, T, X0, N, delta, method = "milstein", 
                                                       theta = c(drift.param[index[j],1], phi[j], sig), sigma.x = expression(sig/(2 * 
                                                                                                                           sqrt(x))), sigma = expression(sig * sqrt(x)), model = "CIR"))
                  }
                  
                }
            }
            
            if (mixture == 0) {
                phi <- drift.param[2] + drift.param[3] * rnorm(M, mean = 0, sd = 1)
                
                # simulation of the series
                if (model == "OU") {
                  for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                                                       method = "EA", theta = c(drift.param[1], phi[j], sig), model = "OU"))
                  }
                }
                
                if (model == "CIR") {
                  for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(t0, T, X0, N, delta, method = "milstein", 
                                                       theta = c(drift.param[1], phi[j], sig), sigma.x = expression(sig/(2 * 
                                                                                                                           sqrt(x))), sigma = expression(sig * sqrt(x)), model = "CIR"))
                  }
                  
                }
            }
            
           
        }
        
    }
    
    if (diffusion.random == 1) {
        # simulation of the random effect of the diffusion coefficient
        
        psi <- 1/sqrt(rgamma(M, shape = diffusion.param[1], rate = 1/diffusion.param[2]))
        
        if ((sum(drift.random) > 2)) {
            # simulation of the random effects
            phi <- matrix(0, 2, M)
            
            if (mixture == 0) {
                for (j in 1:M) {
                  phi[1, j] <- rnorm(1, drift.param[1], sd = drift.param[2] * psi[j])
                  phi[2, j] <- rnorm(1, drift.param[3], sd = drift.param[4] * psi[j])
                }
            }
            
            
            # simulation of the series
            if (model == "OU") {
                for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(phi[, j], psi[j]), model = "OU"))
                }
            }
            
            if (model == "CIR") {
                
                message("to consider a random effect in the diffusion coefficient, alphaj should not be random")
                
            }
        }
        
        if (sum(drift.random) == 1) {
            
            # simulation of the random effects
            phi <- rep(0, M)
            
            if (mixture == 0) {
                for (j in 1:M) {
                  phi[j] <- drift.param[1] + drift.param[2] * psi[j] * rnorm(1, mean = 0, 
                    sd = 1)
                }
            }
            
            # simulation of the series
            if (model == "OU") {
                for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(phi[j], drift.param[3], psi[j]), model = "OU"))
                }
            }
            
            if (model == "CIR") {
                
                message("to consider a random effect in the diffusion coefficient, alphaj should not be random")
                
            }
            
        }
        
        if (sum(drift.random) == 2) {
            # simulation of the random effects
            phi <- rep(0, M)
            if (mixture == 0) {
                for (j in 1:M) {
                  phi[j] <- drift.param[2] + drift.param[3] * psi[j] * rnorm(1, mean = 0, 
                    sd = 1)
                }
            }
            
            # simulation of the series
            if (model == "OU") {
                for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(drift.param[1], phi[j], psi[j]), model = "OU"))
                }
            }
            
            if (model == "CIR") {
                for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(t0, T, X0, N, delta, method = "milstein", 
                      theta = c(drift.param[1], phi[j], psi[j]), sigma.x = expression(psi[j]/(2 * 
                        sqrt(x))), sigma = expression(psi[j] * sqrt(x)), model = "CIR"))
                }
                
            }
        }
        if (sum(drift.random) == 0) {
            phi <- NA
            # simulation of the series
            if (model == "OU") {
                for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(drift.param[1], drift.param[2], psi[j]), 
                      model = "OU"))
                }
            }
            
            if (model == "CIR") {
                for (j in 1:M) {
                    suppressMessages(X[j, ] <- sde.sim(t0, T, X0, N, delta, method = "milstein", 
                      theta = c(drift.param[1], drift.param[2], psi[j]), sigma.x = expression(psi[j]/(2 * 
                        sqrt(x))), sigma = expression(psi[j] * sqrt(x)), model = "CIR"))
                }
                
            }
        }
    }
    
    output <- X
    if (op.plot) {
        
        if (add.plot) {
            for (j in 1:M) {
                lines(delta * (0:N), X[j, ], col = j)
            }
        } else {
            plot(delta * (0:N), X[1, ], type = "l", ylim = c(min(X), max(X)), xlab = "time", 
                ylab = "", col = 1)
            for (j in 2:M) {
                lines(delta * (0:N), X[j, ], col = j)
            }
        }
    }
    return(list(index = index, phi = phi, psi = psi, X = X, times = times))
}
