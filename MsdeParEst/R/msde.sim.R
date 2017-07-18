#' Simulation Of A Mixed Stochastic Differential Equation
#' 
#' @description Simulation of M independent trajectories of a mixed stochastic differential equation (SDE) with linear drift  and two random effects \eqn{(\alpha_j, \beta_j)}
#'  \deqn{dX_j(t)= (\alpha_j- \beta_j X_(t))dt + \sigma a(X_j(t)) dW_j(t)  ,  j=1, ..., M.}
#' @param M number of trajectories.
#' @param T horizon of simulation.
#' @param N number of simulation steps, default Tx100. 
#' @param model name of the SDE: 'OU' (Ornstein-Uhlenbeck) or 'CIR' (Cox-Ingersoll-Ross).
#' @param drift.random random effects in the drift: 0 if no random effect, 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects.
#' @param diffusion.random random effect in the diffusion coefficient: 0 if no random effect, 1 if one multiplicative random effect.
#' @param density.phi name of the density of the random effects.
#' @param drift.fixed fixed effects in the drift: value of the fixed effect when there is only one random effect, 0 otherwise. 
#' If drift.random =2, fixed can be 0 but \eqn{\beta} has to be a non negative random variable for the estimation.
#' @param drift.param vector of parameters of the distribution of the random effects in the drift.
#' @param diffusion.param diffusion parameter if the diffusion coefficient is fixed, vector of parameters of the distribution of the diffusion random effect otherwise. 
#' @param nb.mixt number of mixture components if the drift random effects follow a mixture distribution, default nb.mixt=1.
#' @param mixt.prop vector of mixture proportions if the drift random effects follow a mixture distribution, default mixt.prop=1.
#' @param t0 time origin, default 0.
#' @param X0 initial value of the process, default X0=0.
#' @param invariant 1 if the initial value is simulated from the invariant distribution, default 0.01 and X0 is fixed.
#' @param delta time step of the simulation (T/N).
#' @param op.plot 1 if a plot of the trajectories is required, default 0. 
#' @param add.plot 1 for add trajectories to an existing plot
#' @return
#' \item{X}{matrix (M x (N+1)) of the M trajectories. }
#' \item{phi}{vector (or matrix) of the M simulated random effects.}
#' @details
#' Simulation of M independent trajectories of the SDE (the Brownian motions \eqn{Wj} are independent), with linear drift. Two diffusions are implemented, with one or two random effects:
#' \subsection{Ornstein-Uhlenbeck model (OU)}{
#' If random = 1, \eqn{\beta} is a fixed effect: \eqn{dX_j(t)= (\alpha_j- \beta X_j(t))dt + \sigma_j dW_j(t)  } 
#' 
#' If random = 2, \eqn{\alpha} is a fixed effect: \eqn{dX_j(t)= (\alpha - \beta_j X_j(t))dt + \sigma_j dW_j(t)  }
#' 
#' If random = c(1,2), \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j dW_j(t)  } 
#' }
#' \subsection{Cox-Ingersoll-Ross model (CIR)}{
#' If random = 1, \eqn{\beta} is a fixed effect: \eqn{dX_j(t)= (\alpha_j- \beta X_j(t))dt + \sigma_j \sqrt{X_j(t)} dW_j(t)  } 
#' 
#' If random = 2, \eqn{\alpha} is a fixed effect: \eqn{dX_j(t)= (\alpha - \beta_j X_j(t))dt + \sigma_j \sqrt{X_j(t)} dW_j(t)  } 
#' 
#' If random = c(1,2), \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j \sqrt{X_j(t)}  dW_j(t)  } 
#'}
#' The initial value of each trajectory can be simulated from the invariant distribution of the process:
#'  Normal distribution with mean \eqn{\alpha/\beta} and variance \eqn{\sigma^2/(2 \beta)} for the OU,  a gamma distribution
#'  \eqn{\Gamma(2\alpha/\sigma^2,  \sigma^2/(2\beta))} for the C-I-R model.
#' 
#' \subsection{Density of the diffusion random effect}{
#' We only consider inverse Gamma distributions for \eqn{\sigma_j^2} and diffusion.param=c(shape,scale). 
#' 
#' }
#' 
#' \subsection{Density of the drift random effects}{
#' Several densities are implemented for the random effects, depending on the number of random effects. 
#' 
#' \emph{If two random effects, choice between} 
#' 
#' 'normixtnormixt': Mixture of \eqn{N} normal distributions for both  \eqn{\alpha} \eqn{\beta} and param=c(mean_\eqn{\alpha}^\eqn{1}, sd_\eqn{\alpha}^\eqn{1}, ..., mean_\eqn{\alpha}^\eqn{N}, sd_\eqn{\alpha}^\eqn{N},mean_\eqn{\beta}^\eqn{1}, sd_\eqn{\beta}^\eqn{1}, ..., mean_\eqn{\beta}^\eqn{N}, sd_\eqn{\beta}^\eqn{N})
#' 
#' 'normalnormal': Normal distributions for both  \eqn{\alpha} \eqn{\beta} and param=c(mean_\eqn{\alpha}, sd_\eqn{\alpha}, mean_\eqn{\beta}, sd_\eqn{\beta})
#' 
#' 'gammagamma': Gamma distributions for both  \eqn{\alpha} \eqn{\beta} and param=c(shape_\eqn{\alpha}, scale_\eqn{\alpha}, shape_\eqn{\beta}, scale_\eqn{\beta})
#' 
#' 'gammainvgamma': Gamma for \eqn{\alpha}, Inverse Gamma for \eqn{\beta} and param=c(shape_\eqn{\alpha}, scale_\eqn{\alpha}, shape_\eqn{\beta}, scale_\eqn{\beta})
#' 
#' 'normalgamma':  Normal for \eqn{\alpha}, Gamma for \eqn{\beta} and param=c(mean_\eqn{\alpha}, sd_\eqn{\alpha}, shape_\eqn{\beta}, scale_\eqn{\beta})
#' 
#' 'normalinvgamma':  Normal for \eqn{\alpha}, Inverse Gamma for \eqn{\beta} and param=c(mean_\eqn{\alpha}, sd_\eqn{\alpha}, shape_\eqn{\beta}, scale_\eqn{\beta})
#' 
#' 'gammagamma2':  Gamma \eqn{+2 * \sigma^2} for \eqn{\alpha},  Gamma \eqn{+ 1} for \eqn{\beta}  and param=c(shape_\eqn{\alpha}, scale_\eqn{\alpha}, shape_\eqn{\beta}, scale_\eqn{\beta})
#' 
#' 'gammainvgamma2':   Gamma \eqn{+2 * \sigma^2} for \eqn{\alpha}, Inverse Gamma for \eqn{\beta} and param=c(shape_\eqn{\alpha}, scale_\eqn{\alpha}, shape_\eqn{\beta}, scale_\eqn{\beta})
#' 
#' \emph{If only \eqn{\alpha} is random, choice between}
#' 
#' 'normixt': Mixture of \eqn{N} normal distributions and param=c(mean_\eqn{\alpha}^\eqn{1}, sd_\eqn{\alpha}^\eqn{1}, ..., mean_\eqn{\alpha}^\eqn{N}, sd_\eqn{\alpha}^\eqn{N})
#' 
#' 'normal': Normal distribution with param=c(mean, sd)
#' 
#' lognormal': logNormal distribution with param=c(mean, sd)
#' 
#' 'mixture.normal': mixture of normal distributions \eqn{p N(\mu1,\sigma1^2) + (1-p)N(\mu2, \sigma2^2)} with 
#' param=c(p, \eqn{\mu1, \sigma1, \mu2, \sigma2})
#' 
#' 'gamma': Gamma distribution with param=c(shape, scale) 
#' 
#' 'mixture.gamma': mixture of Gamma distribution \eqn{p \Gamma(shape1,scale1) + (1-p)\Gamma(shape2,scale2)}
#' with param=c(p, shape1, scale1, shape2, scale2)
#' 
#' 'gamma2': Gamma distribution \eqn{+2 * \sigma^2} with param=c(shape, scale) 
#' 
#' 'mixed.gamma2': mixture of Gamma distribution \eqn{p \Gamma(shape1,scale1) + (1-p) \Gamma(shape2,scale2)} + \eqn{+2 * \sigma^2}
#' with param=c(p, shape1, scale1, shape2, scale2)
#' 
#' \emph{If only \eqn{\beta} is random, choice between}
#' 
#'  'normixt': Mixture of \eqn{N} normal distributions and param=c(mean_\eqn{\beta}^\eqn{1}, sd_\eqn{\beta}^\eqn{1}, ..., mean_\eqn{\beta}^\eqn{N}, sd_\eqn{\beta}^\eqn{N})
#' 
#'  'normal': Normal distribution with param=c(mean, sd)
#'  
#'  'gamma': Gamma distribution with param=c(shape, scale) 
#'  
#'  'mixture.gamma': mixture of Gamma distribution \eqn{p \Gamma(shape1,scale1) + (1-p) \Gamma(shape2,scale2)}
#' with param=c(p, shape1, scale1, shape2, scale2)
#' 
#' }
#'
#' @examples 
#' #Simulation of 5 trajectories of the OU SDE with random =1, and a Gamma distribution.
#' 
#' simuOU <- mixedsde.sim(M=5, T=10,N=1000,model='OU', random=1,fixed=0.5,
#' density.phi='gamma', param=c(1.8, 0.8) , sigma=0.1,op.plot=1)
#' X <- simuOU$X ; 
#' phi <- simuOU$phi
#' hist(phi)
#' @references This function mixedsde.sim is based on the package sde, function sde.sim. See Simulation and Inference for stochastic differential equation, S.Iacus, \emph{Springer Series in Statistics 2008}
#' Chapter 2
#' @seealso \url{http://cran.r-project.org/package=sde}
#' 



msde.sim <- function(M, T, N = 100, model, drift.random, diffusion.random, density.phi, 
    drift.fixed = 0, drift.param, diffusion.param, nb.mixt = 1, mixt.prop = 1, t0 = 0, 
    X0 = 0.01, invariant = 0, delta = T/N, op.plot = 0, add.plot = FALSE) {
    
    
    ## local sde.sim to sink undesired output away into a tempfile
    con <- file(tempfile(), open = "w")
    on.exit(close(con))
    sde.sim <- function(...) {
        sink(con)
        res <- sde::sde.sim(...)
        sink(NULL)
        res
    }
    
    if (missing(X0) && missing(invariant)) {
        message("be careful, X0 and invariant are missing thus the initial value X0=0.01 is used")
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
            
            if (density.phi == "normalnormal") {
                phi[1, ] <- rnorm(M, drift.param[1], drift.param[2])
                phi[2, ] <- rnorm(M, drift.param[3], drift.param[4])
            }
            if (density.phi == "mixture.normal") {
                phi <- mixture.sim(M, "mixture.normal", drift.param, mixt.prop)
                # phi[2, ] <- mixture.sim(M, 'mixture.normal', drift.param[2,], mixt.prop)
            }
            if (density.phi == "gammagamma") {
                phi[1, ] <- rgamma(M, drift.param[1], scale = drift.param[2])
                phi[2, ] <- rgamma(M, drift.param[3], scale = drift.param[4])
            }
            if (density.phi == "normalgamma") {
                phi[1, ] <- drift.param[1] + drift.param[2] * rnorm(M, mean = 0, sd = 1)
                phi[2, ] <- rgamma(M, drift.param[3], scale = drift.param[4])
            }
            if (density.phi == "normalinvgamma") {
                phi[1, ] <- drift.param[1] + drift.param[2] * rnorm(M, mean = 0, sd = 1)
                phi[2, ] <- 1/rgamma(M, drift.param[3], scale = drift.param[4])
            }
            if (density.phi == "gammainvgamma") {
                phi[1, ] <- rgamma(M, drift.param[1], scale = drift.param[2])
                phi[2, ] <- 1/rgamma(M, drift.param[3], scale = drift.param[4])
            }
            if (density.phi == "gammagamma2") {
                phi[1, ] <- 2 * sig^2 + rgamma(M, drift.param[1], scale = drift.param[2])
                phi[2, ] <- 1 + rgamma(M, drift.param[3], scale = drift.param[4])
            }
            if (density.phi == "gammainvgamma2") {
                phi[1, ] <- 2 * sig^2 + rgamma(M, drift.param[1], scale = drift.param[2])
                phi[2, ] <- 1/rgamma(M, drift.param[3], scale = drift.param[4])
            }
            
            # simulation of the series
            if (model == "OU") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- phi[1, j]/phi[2, j] + (sig/(sqrt(2 * phi[2, j]))) * rnorm(1)
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(phi[, j], sig), model = "OU"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(phi[, j], sig), model = "OU"))
                  }
                }
            }
            
            if (model == "CIR") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- rgamma(1, 2 * phi[1, j]/sig^2, scale = sig^2/(2 * phi[2, j]))
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "milstein", theta = c(phi[, j], sig), model = "CIR", sigma.x = expression(sig/(2 * 
                        sqrt(x))), sigma = expression(sig * sqrt(x))))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "milstein", theta = c(phi[, j], sig), sigma.x = expression(sig/(2 * 
                        sqrt(x))), sigma = expression(sig * sqrt(x)), model = "CIR"))
                  }
                }
                
            }
        }
        
        if (sum(drift.random) == 1) {
            
            # simulation of the random effects
            phi <- rep(0, M)
            if (density.phi == "mixture.normal") {
                phi <- mixture.sim(M, "mixture.normal", drift.param, mixt.prop)
            }
            if (density.phi == "normal") {
                phi <- drift.param[1] + drift.param[2] * rnorm(M, mean = 0, sd = 1)
            }
            if (density.phi == "lognormal") {
                phi <- drift.param[1] + drift.param[2] * rnorm(M, mean = 0, sd = 1)
                phi <- exp(phi)
            }
            if (density.phi == "gamma") {
                phi <- rgamma(M, shape = drift.param[1], scale = drift.param[2])
            }
            if (density.phi == "gamma2") {
                phi <- 2 * sig^2 + rgamma(M, shape = drift.param[1], scale = drift.param[2])
            }
            if (density.phi == "mixture.gamma") {
                phi <- mixture.sim(M, density.phi, drift.param, mixt.prop)
            }
            if (density.phi == "mixture.gamma2") {
                phi <- 2 * sig^2 + mixture.sim(M, density.phi, drift.param)
            }
            
            
            # simulation of the series
            if (model == "OU") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- phi[j]/drift.fixed + (sig/(sqrt(2 * drift.fixed))) * rnorm(1)
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(phi[j], drift.fixed, sig), model = "OU"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(phi[j], drift.fixed, sig), model = "OU"))
                  }
                }
            }
            
            if (model == "CIR") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- rgamma(1, 2 * phi[j]/sig^2, scale = sig^2/(2 * drift.fixed))
                    suppressMessages(X[j, ] <- sde.sim(T = T, N = N, X0 = X0, delta = delta, 
                      method = "milstein", theta = c(phi[j], drift.fixed, sig), sigma.x = expression(sig/(2 * 
                        sqrt(x))), sigma = expression(sig * sqrt(x)), model = "CIR"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, N = N, X0 = X0, delta = delta, 
                      method = "milstein", theta = c(phi[j], drift.fixed, sig), sigma.x = expression(sig/(2 * 
                        sqrt(x))), sigma = expression(sig * sqrt(x)), model = "CIR"))
                  }
                }
            }
            
        }
        
        
        if (sum(drift.random) == 2) {
            # simulation of the random effects
            phi <- rep(0, M)
            if (density.phi == "mixture.normal") {
                phi <- mixture.sim(M, "mixture.normal", drift.param, mixt.prop)
            }
            if (density.phi == "normal") {
                phi <- drift.param[1] + drift.param[2] * rnorm(M, mean = 0, sd = 1)
            }
            if (density.phi == "gamma") {
                phi <- rgamma(M, drift.param[1], scale = drift.param[2])
            }
            if (density.phi == "mixture.gamma") {
                phi <- mixture.sim(M, density.phi, drift.param, mixt.prop)
            }
            # simulation of the series
            if (model == "OU") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- (sig/(sqrt(2 * phi[j]))) * rnorm(1)
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(drift.fixed, phi[j], sig), model = "OU"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(drift.fixed, phi[j], sig), model = "OU"))
                  }
                }
            }
            
            if (model == "CIR") {
                for (j in 1:M) {
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(t0, T, X0, N, delta, method = "milstein", 
                      theta = c(drift.fixed, phi[j], sig), sigma.x = expression(sig/(2 * 
                        sqrt(x))), sigma = expression(sig * sqrt(x)), model = "CIR"))
                  }
                  if (invariant == 1) {
                    if (drift.fixed == 0) {
                      message("no invariant distribution, please fix X0")
                    }
                    if (drift.fixed != 0) {
                      X0 <- rgamma(1, 2 * drift.fixed/sig^2, scale = sig^2/(2 * phi[j]))
                      suppressMessages(X[j, ] <- sde.sim(t0 = t0, T = T, X0 = X0, N, delta = delta, 
                        method = "milstein", theta = c(drift.fixed, phi[j], sig), sigma.x = expression(sig/(2 * 
                          sqrt(x))), sigma = expression(sig * sqrt(x)), model = "CIR"))
                    }
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
            
            if (density.phi == "normalnormal") {
                for (j in 1:M) {
                  phi[1, j] <- rnorm(1, drift.param[1], sd = drift.param[2] * psi[j])
                  phi[2, j] <- rnorm(1, drift.param[3], sd = drift.param[4] * psi[j])
                }
            }
            
            
            # simulation of the series
            if (model == "OU") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- phi[1, j]/phi[2, j] + (psi[j]/(sqrt(2 * phi[2, j]))) * rnorm(1)
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(phi[, j], psi[j]), model = "OU"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(phi[, j], psi[j]), model = "OU"))
                  }
                }
            }
            
            if (model == "CIR") {
                
                message("to consider a random effect in the diffusion coefficient, alphaj should not be random")
                
            }
        }
        
        if (sum(drift.random) == 1) {
            
            # simulation of the random effects
            phi <- rep(0, M)
            
            if (density.phi == "normal") {
                for (j in 1:M) {
                  phi[j] <- drift.param[1] + drift.param[2] * psi[j] * rnorm(1, mean = 0, 
                    sd = 1)
                }
            }
            
            # simulation of the series
            if (model == "OU") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- phi[j]/drift.fixed + (psi[j]/(sqrt(2 * drift.fixed))) * rnorm(1)
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(phi[j], drift.fixed, psi[j]), model = "OU"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(phi[j], drift.fixed, psi[j]), model = "OU"))
                  }
                }
            }
            
            if (model == "CIR") {
                
                message("to consider a random effect in the diffusion coefficient, alphaj should not be random")
                
            }
            
        }
        
        if (sum(drift.random) == 2) {
            # simulation of the random effects
            phi <- rep(0, M)
            if (density.phi == "normal") {
                for (j in 1:M) {
                  phi[j] <- drift.param[1] + drift.param[2] * psi[j] * rnorm(1, mean = 0, 
                    sd = 1)
                }
            }
            
            # simulation of the series
            if (model == "OU") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- (psi[j]/(sqrt(2 * phi[j]))) * rnorm(1)
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(drift.fixed, phi[j], psi[j]), model = "OU"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(drift.fixed, phi[j], psi[j]), model = "OU"))
                  }
                }
            }
            
            if (model == "CIR") {
                for (j in 1:M) {
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(t0, T, X0, N, delta, method = "milstein", 
                      theta = c(drift.fixed, phi[j], psi[j]), sigma.x = expression(psi[j]/(2 * 
                        sqrt(x))), sigma = expression(psi[j] * sqrt(x)), model = "CIR"))
                  }
                  if (invariant == 1) {
                    if (drift.fixed == 0) {
                      message("no invariant distribution, please fix X0")
                    }
                    if (drift.fixed != 0) {
                      X0 <- rgamma(1, 2 * drift.fixed/psi[j]^2, scale = psi[j]^2/(2 * phi[j]))
                      suppressMessages(X[j, ] <- sde.sim(t0 = t0, T = T, X0 = X0, N, delta = delta, 
                        method = "milstein", theta = c(drift.fixed, phi[j], psi[j]), sigma.x = expression(psi[j]/(2 * 
                          sqrt(x))), sigma = expression(psi[j] * sqrt(x)), model = "CIR"))
                    }
                  }
                }
                
            }
        }
        if (sum(drift.random) == 0) {
            phi <- NA
            # simulation of the series
            if (model == "OU") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- (psi[j]/(sqrt(2 * drift.fixed[2]))) * rnorm(1)
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(drift.fixed[1], drift.fixed[2], psi[j]), 
                      model = "OU"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = delta, 
                      method = "EA", theta = c(drift.fixed[1], drift.fixed[2], psi[j]), 
                      model = "OU"))
                  }
                }
            }
            
            if (model == "CIR") {
                for (j in 1:M) {
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(t0, T, X0, N, delta, method = "milstein", 
                      theta = c(drift.fixed[1], drift.fixed[2], psi[j]), sigma.x = expression(psi[j]/(2 * 
                        sqrt(x))), sigma = expression(psi[j] * sqrt(x)), model = "CIR"))
                  }
                  if (invariant == 1) {
                    if (drift.fixed[2] == 0) {
                      message("no invariant distribution, please fix X0")
                    }
                    if (drift.fixed[2] != 0) {
                      X0 <- rgamma(1, 2 * drift.fixed[1]/psi[j]^2, scale = psi[j]^2/(2 * 
                        drift.fixed[2]))
                      suppressMessages(X[j, ] <- sde.sim(t0 = t0, T = T, X0 = X0, N, delta = delta, 
                        method = "milstein", theta = c(drift.fixed[1], drift.fixed[2], 
                          psi[j]), sigma.x = expression(psi[j]/(2 * sqrt(x))), sigma = expression(psi[j] * 
                          sqrt(x)), model = "CIR"))
                    }
                  }
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
