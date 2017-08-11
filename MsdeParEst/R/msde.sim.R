# MsdeParEst R package ; file msde.sim.r (last modified: 2017-08-11)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017
# UMR 518 AgroParisTech/INRA, 16 rue Claude Bernard, 75 231 Paris Cedex 05

#' Simulation Of A Mixed Stochastic Differential Equation
#' 
#' @description Simulation of M independent trajectories of a mixed stochastic differential equation (SDE) with linear drift
#'  \deqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j a(X_j(t)) dW_j(t)  ,  j=1, ..., M.}
#' There can be up to two random effects \eqn{(\alpha_j, \beta_j)} in the drift and one random effect \eqn{\sigma_j} in the diffusion coefficient.
#' @param M number of trajectories.
#' @param T horizon of simulation.
#' @param N number of simulation steps, default Tx100. 
#' @param model name of the SDE: 'OU' (Ornstein-Uhlenbeck) or 'CIR' (Cox-Ingersoll-Ross).
#' @param drift.random random effects in the drift: 0 if no random effect, 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects.
#' @param diffusion.random random effect in the diffusion coefficient: 0 if no random effect, 1 if one multiplicative random effect.
#' @param drift.param vector (not mixture) or matrix (mixture) of values of the fixed effects and/or the parameters of the distribution of the random effects in the drift (see details).
#' @param diffusion.param diffusion parameter if the diffusion coefficient is fixed, vector of parameters \eqn{c(a,\lambda)} of the distribution of the diffusion random effect otherwise. 
#' @param nb.mixt number of mixture components if the drift random effects follow a mixture distribution, default nb.mixt=1.
#' @param mixt.prop vector of mixture proportions if the drift random effects follow a mixture distribution, default mixt.prop=1.
#' @param t0 time origin, default 0.
#' @param X0 initial value of the process, default X0=0.001.
#' @param delta time step of the simulation (T/N).
#' @param op.plot 1 if a plot of the trajectories is required, default 0. 
#' @param add.plot 1 for add trajectories to an existing plot
#' @return
#' \item{X}{matrix (M x (N+1)) of the M trajectories. }
#' \item{times}{vector of the N+1 simulated observation times from t0 to T.}
#' \item{phi}{vector (or matrix) of the M simulated random effects of the drift.}
#' \item{psi}{vector of the M simulated values of \eqn{\sigma_j}.}
#' 
#' @importFrom sde sde.sim
#' @importFrom stats rnorm
#' 
#' @examples 
#' 
#'  \dontrun{
#'  # Example : one random effect in the drift and one fixed effect in the diffusion coefficient
#'  sim <- msde.sim(M = 100, T = 5, N = 5000, model = 'OU', drift.random = 2,
#'                 diffusion.random = 0, drift.param = c(0,1,sqrt(0.4/4)), diffusion.param = 0.5)
#'                 
#'  # Example : two random effects in the drift and one random effect in the diffusion coefficient
#'
#'  sim <- msde.sim(M = 100, T = 5, N = 5000, model = 'OU', drift.random = c(1,2),
#'                 diffusion.random = 1, drift.param = c(1,0.5,0.5,0.5), diffusion.param = c(8,1/2))
#'                 
#'  # Example : one fixed effect and one mixture random effect in the drift, and one fixed effect in
#'  # the diffusion coefficient
#'  
#'  sim <- msde.sim(M = 100, T = 5, N = 5000, model = 'OU', 
#'                  drift.random = 1, drift.param = matrix(c(0.5,1.8,0.25,0.25,1,1),nrow=2,byrow=F),
#'                  diffusion.random = 0, diffusion.param = 0.1, 
#'                  nb.mixt = 2, mixt.prop = c(0.5,0.5))
#'  } 
#' 
#' @details
#' Simulation of N discrete observations on time interval [t0,T] of M independent trajectories of the SDE 
#' \deqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j \ a(X_j(t)) dW_j(t),}
#' \eqn{j=1,\ldots,M}, where the \eqn{(W_j(t))} are independant Wiener processes. 
#' 
#' \bold{Specification of \eqn{\alpha,\beta,\sigma}}
#' 
#' The diffusion includes either a fixed effect or a random effect:
#' \enumerate{
#' \item if diffusion.random = 0: \eqn{\sigma_j \equiv \sigma} is fixed, and diffusion.param = \eqn{\sigma}.
#' In this case, the drift includes no, one or two random effects: 
#' \enumerate{
#' \item if drift.random = 0: \eqn{\alpha_j \equiv \alpha} and \eqn{\beta_j \equiv \beta} are fixed, and drift.param=c\eqn{(\alpha,\beta)}
#' \item if drift.random = 1: \eqn{\alpha_j} is random with distribution \eqn{N(\mu_{\alpha},\omega^2_{\alpha})} whereas \eqn{\beta_j \equiv \beta} is fixed, and drift.param=c\eqn{(\mu_{\alpha},\omega^2_{\alpha},\beta)}
#' \item if drift.random = 2: \eqn{\alpha_j \equiv \alpha} is fixed and \eqn{\beta_j} is random with distribution \eqn{N(\mu_{\beta},\omega^2_{\beta})},  and drift.param=c\eqn{(\alpha, \mu_{\beta},\omega^2_{\beta})}
#' \item if drift.random = c(1,2): \eqn{\alpha_j} and \eqn{\beta_j} are random with distributions \eqn{N(\mu_{\alpha},\omega^2_{\alpha})} and \eqn{N(\mu_{\beta},\omega^2_{\beta})} respectively,
#' and drift.param = c\eqn{(\mu_{\alpha},\omega^2_{\alpha},\mu_{\beta},\omega^2_{\beta})}
#' }
#' \item if diffusion.random = 1: \eqn{\sigma_j} is random such that \eqn{1/\sigma_j^2 \sim \Gamma}, and drift.param=c\eqn{(a,\lambda)}.
#' In this case, the drift includes at least one random effect:
#' \enumerate{
#' \item if drift.random = 1: \eqn{\alpha_j} is random with distribution \eqn{N(\mu_{\alpha}, \sigma_j^2 \omega^2_{\alpha})} whereas \eqn{\beta_j \equiv \beta} is fixed, and drift.param=c\eqn{(\mu_{\alpha},\omega^2_{\alpha},\beta)}
#' \item if drift.random = 2: \eqn{\alpha_j \equiv \alpha} is fixed and \eqn{\beta_j} is random with distribution \eqn{N(\mu_{\beta},\sigma_j^2 \omega^2_{\beta})},  and drift.param=c\eqn{(\alpha, \mu_{\beta},\omega^2_{\beta})}
#' \item if drift.random = c(1,2): \eqn{\alpha_j} and \eqn{\beta_j} are random with distributions \eqn{N(\mu_{\alpha},\sigma_j^2 \omega^2_{\alpha})} and \eqn{N(\mu_{\beta}, \sigma_j^2 \omega^2_{\beta})} respectively,
#' and drift.param = c\eqn{(\mu_{\alpha},\omega^2_{\alpha},\mu_{\beta},\omega^2_{\beta})}
#' }
#' 
#' If the random effects in the drift follow a mixture distribution (nb.mixt=K, K>1), drift.param is a matrix instead of a vector. Each line of the matrix
#' contains, as above, the parameter values for each mixture component. 
#' }
#' 
#'
#' @references This function mixedsde.sim is based on the package sde, function sde.sim. See 
#' 
#' Simulation and Inference for stochastic differential equation, S.Iacus, \emph{Springer Series in Statistics 2008}
#' Chapter 2
#' @seealso \url{http://cran.r-project.org/package=sde}
#' 



msde.sim <- function(M, T, N = 100, model, drift.random, diffusion.random, 
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
    
    ## Stops and warnings
    
    if (((M - round(M)) != 0) || (M <= 0)){
      stop("Invalid value for M. The number of simulated trajectories should be a positive integer")
    }
    
    if (((N - round(N)) != 0) || (N <= 0)){
      stop("Invalid value for N. The number of simulated time points should be a positive integer")
    }
    
    if (t0 < 0){
      stop("Invalid value for t0. The starting time point should be positive")
    }
    
    if (T <= 0){
      stop("Invalid value for T. The maximum time point should be positive")
    }
    
    if (t0 > T){
      stop("The starting time point t0 should be smaller than the maximum time point T")
    }
    
    if ((delta < 0) || (delta > (T-t0))){
      stop("Invalid value for delta")    
    }

        if((model != 'OU')&(model != 'CIR')){stop("A model must be precised: OU or CIR")}
    
    if (((nb.mixt - round(nb.mixt)) != 0) || (nb.mixt <= 0)){
      stop("The number of mixture components (nb.mixt) should be a positive integer")
    } 
    
    if ((diffusion.random!=0) && (diffusion.random!=1)){
      stop("Invalid value for diffusion.random, should be either 0 or 1")
    }
    
    valid.drift = list(0,1,2,c(1,2))
    
    check <- 0
    for (i in 1:4){
      check <- check + prod(drift.random %in% valid.drift[[i]])
    }
    
    if (check == 0){
      stop("Invalid value for drift.random, should be either 0, or 1, or 2, or c(1,2)")  
    }
    
    if ((model == "CIR") && (X0 <= 0)){
      stop('For the CIR model, the initial values should be positive.')
    }
    
    
    if (missing(X0)) {
        message("Be careful, X0 is missing thus the initial value X0=0.01 is used")
    }
    
    if ((diffusion.random == 1) && (length(diffusion.param) != 2)){
      stop("Invalid diffusion.param, should be a vector of length 2")
    }
    
    if ((diffusion.random == 0) && (length(diffusion.param) > 1)){
      message("Only the first value of diffusion.param is considered")
    }
    
    if (prod(diffusion.param > 0) != 1){
      stop("Invalid diffusion.param, should contain positive values")
    }
    
    if (nb.mixt == 1){
      if ((drift.random == 0) && (length(drift.param) != 2)){
        stop("Invalid drift.param, should be a vector of length 2")      
      }
      
      if (((sum(drift.random) == 1) || (sum(drift.random) == 2)) && (length(drift.param) != 3)){
        stop("Invalid drift.param, should be a vector of length 3")      
      }
      
      if ((sum(drift.random) == 3) && (length(drift.param) != 4)) {
        stop("Invalid drift.param, should be a vector of length 4") 
      }
    }else {
      if (dim(drift.param)[1] != nb.mixt){
        stop("Invalid dimensions for drift.param, should have as many lines as mixture components")
      }
      if (((sum(drift.random) == 1) || (sum(drift.random) == 2)) && (dim(drift.param)[2] != 3)){
        stop("Invalid dimensions for drift.param, should have 3 columns")
      }
      if ((sum(drift.random) == 3) && (dim(drift.param)[2] != 4)) {
        stop("Invalid dimensions for drift.param, should have 4 columns")
      }
      
    }
    
    
    if (((sum(drift.random)) == 0) && (diffusion.random == 0)) {
      stop("There should be at least one random effect either in the drift or the diffusion coefficient.")
    }
    
    if (((diffusion.random == 1) && (nb.mixt > 1))) {
      stop("If there is one random effect in the diffusion coefficient, the random effects in the drift can't follow a mixture of Normal distributions. Try nb.mixt = 1.")
    }
    
    if (length(mixt.prop) != nb.mixt) {
      stop("There should be as many mixing proportions as mixture components")
    }
    
    if ((prod(mixt.prop >= 0) == 0)) {
      stop("Invalid values for the mixing proportions, should be positive")
    } else {
      if (prod(mixt.prop <= 1) == 0) {
        mixt.prop <- mixt.prop/sum(mixt.prop)
      }
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
            
            if (nb.mixt == 1) {
                phi[1, ] <- rnorm(M, drift.param[1], drift.param[2])
                phi[2, ] <- rnorm(M, drift.param[3], drift.param[4])
            }
            if (nb.mixt > 1) {
                phi <- mixture.sim(M, drift.param, mixt.prop)$Y
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
            if (nb.mixt > 1) {
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
            
            
            if (nb.mixt == 1) {
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
            if (nb.mixt > 1) {
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
            
            if (nb.mixt == 1) {
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
            
            if (nb.mixt == 1) {
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
            
            if (nb.mixt == 1) {
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
            if (nb.mixt == 1) {
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
