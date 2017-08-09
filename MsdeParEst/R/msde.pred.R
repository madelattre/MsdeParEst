# MsdeParEst R package ; file msde.pred.r (last modified: 2017-08-09)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017
# UMR 518 AgroParisTech/INRA, 16 rue Claude Bernard, 75 231 Paris Cedex 05


#' Prediction Of Mixed Stochastic Differential Equations Trajectories
#' 
#' 
#' @description This function proposes to keep two thirds of the data (randomly chosen) to do the parametric
#' estimation of the density of the parameters and of the fixed parameters of model 
#'  \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j a(X_j(t)) dW_j(t)}
#'using the same method as in the estimationr function \code{msde.fit}
#' and then to predict new trajectories from the estimated model. 
#' The plot reflect the adequation between the last third of the data and the simulated one.
#' 
#'  
#' @param times vector of observation times
#' @param X matrix of the M trajectories (each row is a trajectory with as much columns as observations)
#' @param model name of the SDE: 'OU' (Ornstein-Uhlenbeck) or 'CIR' (Cox-Ingersoll-Ross)
#' @param drift.random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects
#' @param drift.fixed default NULL, fixed effect in the drift: value of the fixed effect when there is only one random effect and it is not estimated, NULL otherwise
#' @param diffusion.random default 0, 1 if one random effect in the diffusion, 0 if there is no random effect in the diffusion
#' @param diffusion.fixed default NULL, fixed effect in the diffusion: value of the fixed effect when there is no random effect in the diffusion and it is not estimated, NULL otherwise
#' @param mixture 1 if the random effects in the drift follow a mixture distribution, 0 otherwise. Default to 0.
#' @param nb.mixt default 1, number of mixture components for the distribution of the random effects in the drift 
#' @param Niter default 10, number of iterations for the EM algorithm if mixture = 1
#' @param discrete default 1, 1 for discrete observations, 0 otherwise. If discrete = 0, and diffusion.random = 0, the exact likelihood associated with continuous observations is 
#' discretized. If discrete = 1, the likelihood of the Euler scheme of the mixed SDE is computed. 
#' @param level alpha for the predicion intervals, default 0.05
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param plot.pred logical(1), if TRUE, the results are depicted grafically

#' @return 
#'
#' \item{res}{is the vector of subscript in \eqn{1,...,M} where the estimation of \eqn{phi} has been done,  most of the time \eqn{index= 1:M}}
#' \item{Xpred}{is the vector of subscript in \eqn{1,...,M} where the estimation of \eqn{phi} has been done,  most of the time \eqn{index= 1:M}}
#' \item{indexpred}{matrix of estimators of \eqn{\phi=\alpha, or \beta, or (\alpha,\beta)} from the efficient statitics (see \code{\link{UVS}}), matrix of two lines if drift.random =c(1,2), numerical type otherwise}
#' \item{phipred}{matrix of estimators of \eqn{\psi^2=\sigma^2} from the efficient statistics (see \code{\link{UVS}}), matrix of one line}
#' 

#' @importFrom stats dnorm
#' @importFrom stats rgamma
#' @importFrom stats density
#' @importFrom stats quantile
#' @importFrom methods new
#' @importFrom methods slot
#' @importFrom methods slotNames
#' @importFrom stats kmeans
#' @importFrom stats sd
#' @importFrom graphics lines
#' @importFrom grDevices dev.new
#' 
#' @examples
#'
#' \dontrun{
#' # Example 1: one random effect in the drift and one random effect in the diffusion coefficient.

#'   }
#' 
#' @keywords estimation
#' @references See  
#' Maximum Likelihood Estimation for Stochastic Differential Equations with Random Effects, Delattre, M., Genon-Catalot, V. and Samson, A. \emph{Scandinavian Journal of Statistics 40(2) 2012} \bold{322-343} 
#' 
#' Estimation of population parameters in stochastic differential equations with random effects in the diffusion coefficient, Delattre, M., Genon-Catalot, V. and Samson, A. \emph{ESAIM:PS 19 2015} \bold{671-688}
#' 
#' Mixtures of stochastic differential equations with random effects: application to data clustering, Delattre, M., Genon-Catalot, V. and Samson, A. \emph{Journal of Statistical Planning and Inference 173 2016} \bold{109-124} 
#' 
#' Parametric inference for discrete observations of diffusion processes with mixed effects, Delattre, M., Genon-Catalot, V. and Laredo, C. \emph{hal-01332630 2016}
#' 
#' Estimation of the joint distribution of random effects for a discretely observed diffusion with random effects, Delattre, M., Genon-Catalot, V. and Laredo, C. \emph{hal-01446063 2017}


msde.pred <- function(times, X, model = c("OU", "CIR"), drift.random, drift.fixed = NULL, 
    diffusion.random = 0, diffusion.fixed = NULL, mixture = 0, nb.mixt = 1,
    Niter = 10, discrete = 1, plot.pred = TRUE, level = 0.05, newwindow = FALSE) {
  
    model <- match.arg(model)
    
    if (newwindow) {
         #x11(width = 10)
      dev.new(width = 10)
    }
    if (is.matrix(X)) {
        if (nrow(X) == length(times)) {
            X <- t(X)
        } else {
            if (ncol(X) != length(times)) {
                stop("Length of times has to be equal to the columns of X")
            }
        }
    }
    
    
    ## local sde.sim to sink undesired output away into a tempfile
    con <- file(tempfile(), open = "w")
    on.exit(close(con))
    sde.sim <- function(...) {
        sink(con)
        res <- sde::sde.sim(...)
        sink(NULL)
        res
    }
    
    M <- dim(X)[1]
    K <- dim(X)[2]
    N <- dim(X)[2] - 1
    
    Mpred <- floor(M * 2/3) + 1
    ipred <- sample(seq(1,M,1),Mpred,replace=F)
    Xtrue <- X[-ipred,]
    Xestim <- X[ipred,]
    
    timestrue <- times
    Tend <- timestrue[length(timestrue)]
    delta <- round(diff(timestrue), 10)[1]
    
    times <- seq(0, Tend, by = delta)
    
    res <- msde.fit(times = times, X = Xestim, model = model, drift.random = drift.random, 
        diffusion.random = diffusion.random, mixture = mixture, nb.mixt = nb.mixt, 
        Niter = Niter)

    
    
    if (mixture == 1) {
        
        if (diffusion.random == 1) {
            warning("For the considered mixtures of SDE, there should not be random effects in the diffusion coeffcicient")
            diffusion.random <- 0
        }
        
        sig <- sqrt(res@sigma2)
        
        if (sum(drift.random) == 0) {
            stop("There should be at least one random effect in the drift")
        }
        
        if (sum(drift.random) != 3) {
            
            if (drift.random == 1) {
                paramfixed <- res@mu[Niter, , 2] 
                
                phipred <- rep(0, M - Mpred)
                
                simu <- mixture.sim(M - Mpred, matrix(c(res@mu[Niter, , 1], res@omega[Niter, , 1]),nrow=2,byrow=F), 
                                    res@mixt.prop[Niter, ])
                
                phipred <- simu$Y
                index <- simu$index

                if (model == "OU") {
                  indexpred <- 1:(M - Mpred)
                  Xpred <- matrix(0, M - Mpred, N + 1)
                  for (j in 1:(M - Mpred)) {
                    suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[j, 1], N = N, delta = Tend/N, 
                      method = "EA", theta = c(phipred[j], paramfixed[index[j]], sig), model = "OU"))
                  }
                }
                
                if (model == "CIR") {
                  indexpred <- which(phipred > 0)
                  phipred <- phipred[indexpred]
                  Mpred <- length(phipred)
                  Xpred <- matrix(0, M - Mpred, N + 1)
                  
                  for (j in 1:(M - Mpred)) {
                    
                    suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[indexpred[j], 1], N = N, delta = Tend/N, 
                      method = "milstein", theta = c(phipred[j], paramfixed[indexpred[j]], sig), model = "CIR", sigma.x = expression(sig/(2 * 
                        sqrt(x))), sigma = expression(sig * sqrt(x))))
                    
                  }
                }
                phipred <- matrix(phipred, 1, length(phipred))
            }
            
            if (drift.random == 2) {
                paramfixed <- res@mu[Niter, , 1]
                
                phipred <- rep(0, M - Mpred)
                
                simu <- mixture.sim(M - Mpred, matrix(c(res@mu[Niter, , 2], res@omega[Niter, , 2])
                                                      ,nrow=2,byrow=F), res@mixt.prop[Niter, ])
                
                phipred <- simu$Y
                index <- simu$index

                if (model == "OU") {
                  indexpred <- which(phipred > 0)
                  phipred <- phipred[indexpred]
                  Mprednew <- length(indexpred)
                  Xpred <- matrix(0, Mprednew, N + 1)
                  
                  for (j in 1:Mprednew) {
                    
                    suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[j, 1], N = N, delta = Tend/N, 
                      method = "EA", theta = c(paramfixed[index[j]], phipred[j], sig), model = "OU"))
                  }
                }
                
                if (model == "CIR") {
                  indexpred <- which(phipred > 0)
                  phipred <- phipred[indexpred]
                  Mprednew <- length(phipred)
                  Xpred <- matrix(0, Mprednew, N + 1)
                  
                  for (j in 1:Mprednew) {
                    
                    suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[indexpred[j], 1], N = N, delta = Tend/N, 
                      method = "milstein", theta = c(paramfixed[indexpred[j]], phipred[j], sig), model = "CIR", sigma.x = expression(sig/(2 * 
                        sqrt(x))), sigma = expression(sig * sqrt(x))))
                    
                  }
                }
                phipred <- matrix(phipred, 1, length(phipred))
            }
            
            if (plot.pred == TRUE) {
                op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                  cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)

                plot(timestrue, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 
                  0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
                for (j in indexpred) {
                  lines(timestrue, Xtrue[j, ], col = j)
                }
                
                plot(times, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, max(Xtrue, 
                  Xpred) * 1.5), main = "Predicted trajectories")
                for (j in 1:length(indexpred)) {
                  lines(times, Xpred[j, ], col = 1)
                }
                
                l.bound = level/2
                u.bound = 1 - level/2
                PI <- matrix(0, 2, N + 1)
                for (k in 1:(N + 1)) {
                  PI[, k] <- quantile(Xpred[, k], c(l.bound, u.bound))
                }
                lines(times, PI[1, ], col = "green", lwd = 2)
                lines(times, PI[2, ], col = "green", lwd = 2)
            }
            
        }
        
        
        if (sum(drift.random) == 3) {
            
            phipred <- matrix(0, 2, M - Mpred)
            
            param <- cbind(res@mu[Niter,,1],res@omega[Niter,,1])
            
            for (k in 2:nb.mixt){
              param <- cbind(param,res@mu[Niter,,k],res@omega[Niter,,k])
            }
            
            
            phipred <- mixture.sim(M - Mpred, param, res@mixt.prop[Niter, ])$Y
            
            if (model == "OU") {
                indexpred <- which(phipred[2, ] > 0)
                phipred <- phipred[, indexpred]
                Mprednew <- length(indexpred)
                Xpred <- matrix(0, Mprednew, N + 1)
                for (j in 1:Mprednew) {
                  
                  suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[j, 1], N = N, delta = Tend/N, method = "EA", 
                    theta = c(phipred[, j], sig), model = "OU"))
                  
                  
                }
            }
            if (model == "CIR") {
                indexpred <- which((phipred[1, ] > 0) & (phipred[2, ] > 0))
                phipred <- phipred[, indexpred]
                Mprednew <- length(indexpred)
                Xpred <- matrix(0, Mprednew, N + 1)
                
                for (j in 1:Mprednew) {
                  
                  suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[indexpred[j], 1], N = N, delta = Tend/N, 
                    method = "milstein", theta = c(phipred[, j], sig), model = "CIR", sigma.x = expression(sig/(2 * 
                      sqrt(x))), sigma = expression(sig * sqrt(x))))
                  
                }
            }
            
            if (plot.pred == TRUE) {
                
                op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                  cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
             
                plot(times, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 
                  0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
                for (j in indexpred) {
                  lines(times, Xtrue[j, ], col = j)
                }
                
                plot(timestrue, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, 
                  max(Xtrue, Xpred) * 1.5), main = "Predicted trajectories")
                for (j in 1:length(indexpred)) {
                  lines(timestrue, Xpred[j, ], col = 1)
                }
                
                l.bound = level/2
                u.bound = 1 - level/2
                PI <- matrix(0, 2, N + 1)
                for (k in 1:(N + 1)) {
                  PI[, k] <- quantile(Xpred[, k], c(l.bound, u.bound))
                }
                lines(timestrue, PI[1, ], col = "green", lwd = 2)
                lines(timestrue, PI[2, ], col = "green", lwd = 2)
            }
        }
        
        return(new(Class = "class.mixture.pred", estim = res, phipred = as.matrix(phipred), Xpred = Xpred, indexpred = indexpred))
        
    }
    
    # ###########
    if (mixture == 0) {
        
        if (diffusion.random == 0) {
            sig <- sqrt(res@sigma2)
            
            if (sum(drift.random != 3)) {
                
                if (sum(drift.random) == 1) {
                  paramfixed <- res@mu[2]
                  
                  phipred <- rep(0, M - Mpred)
                  
                  
                  phipred <- rnorm(M - Mpred, res@mu[1], res@omega[1])
                  
                  
                  if (model == "OU") {
                    indexpred <- 1:(M - Mpred)
                    for (j in 1:(M - Mpred)) {
                      suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[j, 1], N = N, delta = Tend/N, 
                        method = "EA", theta = c(phipred[j], paramfixed, sig), model = "OU"))
                    }
                  }
                  
                  if (model == "CIR") {
                    indexpred <- which(phipred > 0)
                    phipred <- phipred[indexpred]
                    Mprednew <- length(phipred)
                    Xpred <- matrix(0, Mprednew, N + 1)
                    
                    for (j in 1:Mprednew) {
                      
                      suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[indexpred[j], 1], N = N, delta = Tend/N, 
                        method = "milstein", theta = c(phipred[j], paramfixed, sig), model = "CIR", sigma.x = expression(sig/(2 * 
                          sqrt(x))), sigma = expression(sig * sqrt(x))))
                      
                      
                    }
                  }
                }
                
                if (sum(drift.random) == 2) {
                  paramfixed <- res@mu[1]
                  
                  phipred <- rep(0, M - Mpred)
                  
                  phipred <- rnorm(M - Mpred, res@mu[2], res@omega[2])
                  
                  
                  if (model == "OU") {
                    indexpred <- which(phipred > 0)
                    phipred <- phipred[indexpred]
                    Mprednew <- length(indexpred)
                    Xpred <- matrix(0, Mprednew, N + 1)
                    
                    for (j in 1:Mprednew) {
                      
                      suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[j, 1], N = N, delta = Tend/N, 
                        method = "EA", theta = c(paramfixed, phipred[j], sig), model = "OU"))
                      
                    }
                  }
                  
                  if (model == "CIR") {
                    indexpred <- which(phipred > 0)
                    phipred <- phipred[indexpred]
                    Mprednew <- length(phipred)
                    Xpred <- matrix(0, Mprednew, N + 1)
                    
                    for (j in 1:Mprednew) {
                      
                      
                      suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[indexpred[j], 1], N = N, delta = Tend/N, 
                        method = "milstein", theta = c(paramfixed, phipred[j], sig), model = "CIR", sigma.x = expression(sig/(2 * 
                          sqrt(x))), sigma = expression(sig * sqrt(x))))
                      
                    }
                  }
                }
                
                if (plot.pred == TRUE) {
                  op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                    cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
                  
                  plot(timestrue, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, 
                    Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
                  for (j in indexpred) {
                    lines(timestrue, Xtrue[j, ], col = j)
                  }
                  
                  plot(times, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, 
                    max(Xtrue, Xpred) * 1.5), main = "Predicted trajectories")
                  for (j in 1:length(indexpred)) {
                    lines(times, Xpred[j, ], col = 1)
                  }
                  
                  l.bound = level/2
                  u.bound = 1 - level/2
                  PI <- matrix(0, 2, N + 1)
                  for (k in 1:(N + 1)) {
                    PI[, k] <- quantile(Xpred[, k], c(l.bound, u.bound))
                  }
                  lines(times, PI[1, ], col = "green", lwd = 2)
                  lines(times, PI[2, ], col = "green", lwd = 2)
                }
            }
            if (sum(drift.random) == 3) {
                
                phipred <- matrix(0, 2, M - Mpred)
                
                phipred[1, ] <- rnorm(M - Mpred, res@mu[1], res@omega[1])
                
                phipred[2, ] <- rnorm(M - Mpred, res@mu[2], res@omega[2])
                
                if (model == "OU") {
                  indexpred <- which(phipred[2, ] > 0)
                  phipred <- phipred[, indexpred]
                  Mprednew <- length(indexpred)
                  Xpred <- matrix(0, Mprednew, N + 1)
                  for (j in 1:Mprednew) {
                    
                    suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[j, 1], N = N, delta = Tend/N, 
                      method = "EA", theta = c(phipred[, j], sig), model = "OU"))
                  }
                }
                if (model == "CIR") {
                  indexpred <- which((phipred[1, ] > 0) & (phipred[2, ] > 0))
                  phipred <- phipred[, indexpred]
                  Mprednew <- length(indexpred)
                  Xpred <- matrix(0, Mprednew, N + 1)
                  
                  for (j in 1:Mprednew) {
                    
                    suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[indexpred[j], 1], N = N, delta = Tend/N, 
                      method = "milstein", theta = c(phipred[, j], sig), model = "CIR", sigma.x = expression(sig/(2 * 
                        sqrt(x))), sigma = expression(sig * sqrt(x))))
                    
                  }
                }
                
                if (plot.pred == TRUE) {
                  
                  op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                    cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)

                  
                  plot(times, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 
                    0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
                  for (j in indexpred) {
                    lines(times, Xtrue[j, ], col = j)
                  }
                  
                  plot(timestrue, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, 
                    max(Xtrue, Xpred) * 1.5), main = "Predicted trajectories")
                  for (j in 1:length(indexpred)) {
                    lines(timestrue, Xpred[j, ], col = 1)
                  }
                  
                  l.bound = level/2
                  u.bound = 1 - level/2
                  PI <- matrix(0, 2, N + 1)
                  for (k in 1:(N + 1)) {
                    PI[, k] <- quantile(Xpred[, k], c(l.bound, u.bound))
                  }
                  lines(timestrue, PI[1, ], col = "green", lwd = 2)
                  lines(timestrue, PI[2, ], col = "green", lwd = 2)
                }
            }
        }
        
        
        
        if (diffusion.random == 1) {
            
            psipred <- rep(0, M - Mpred)
            psipred <- 1/sqrt(rgamma(M - Mpred, shape = res@a, rate = 1/res@lambda))
            
            if (sum(drift.random) != 3) {
               
                if (sum(drift.random) == 0) {
                  phipred <- 0
                  paramfixed <- res@mu
                  
                  if (model == "OU") {
                    indexpred <- 1:(M - Mpred)
                    Xpred <- matrix(0, (M-Mpred), N + 1)
                    for (j in 1:(M - Mpred)) {
                      suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[j, 1], N = N, delta = Tend/N, 
                                                             method = "EA", theta = c(paramfixed, psipred[j]), model = "OU"))
                      
                    }
                  }
                  if (model == "CIR") {
                    indexpred <- which(phipred > 0)
                    phipred <- phipred[indexpred]
                    psipred <- psipred[indexpred]
                    Mprednew <- length(phipred)
                    Xpred <- matrix(0, Mprednew, N + 1)
                    
                    for (j in 1:Mprednew) {
                      
                      suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[indexpred[j], 1], N = N, delta = Tend/N, 
                                                             method = "milstein", theta = c(paramfixed, psipred[j]), model = "CIR", sigma.x = expression(psipred[j]/(2 * 
                                                                                                                                                                                   sqrt(x))), sigma = expression(psipred[j] * sqrt(x))))
                      
                      
                    }
                  }
                  
                }
              
                if (sum(drift.random) == 1) {
                  
                  paramfixed <- res@mu[2]
                  
                  phipred <- rep(0, M - Mpred)
                  
                  for (j in 1:(M-Mpred)) {
                    phipred[j] <- rnorm(1, mean = res@mu[1], sd = res@omega[1] * psipred[j])
                  }
                  
                  
                  if (model == "OU") {
                    indexpred <- 1:(M - Mpred)
                    Xpred <- matrix(0, (M-Mpred), N + 1)
                    for (j in 1:(M - Mpred)) {
                      suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[j, 1], N = N, delta = Tend/N, 
                        method = "EA", theta = c(phipred[j], paramfixed, psipred[j]), model = "OU"))
                      
                    }
                  }
                  if (model == "CIR") {
                    indexpred <- which(phipred > 0)
                    phipred <- phipred[indexpred]
                    psipred <- psipred[indexpred]
                    Mprednew <- length(phipred)
                    Xpred <- matrix(0, Mprednew, N + 1)
                    
                    for (j in 1:Mprednew) {
                      
                      suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[indexpred[j], 1], N = N, delta = Tend/N, 
                        method = "milstein", theta = c(phipred[j], paramfixed, psipred[j]), model = "CIR", sigma.x = expression(psipred[j]/(2 * 
                          sqrt(x))), sigma = expression(psipred[j] * sqrt(x))))
                      
                      
                    }
                  }
                }
                
                if (sum(drift.random) == 2) {
                  paramfixed <- res@mu[1]
                  
                  phipred <- rep(0, M - Mpred)
                  
                  for (j in 1:(M-Mpred)) {
                    phipred[j] <- rnorm(1, res@mu[2], res@omega[2] * psipred[j])
                  }
                  
                  
                  if (model == "OU") {
                    indexpred <- which(phipred > 0)
                    phipred <- phipred[indexpred]
                    psipred <- psipred[indexpred]
                    Mprednew <- length(indexpred)
                    Xpred <- matrix(0, Mprednew, N + 1)
                    
                    for (j in 1:Mprednew) {
                      
                      suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[j, 1], N = N, delta = Tend/N, 
                        method = "EA", theta = c(paramfixed, phipred[j], psipred[j]), model = "OU"))
                      
                    }
                  }
                  
                  if (model == "CIR") {
                    indexpred <- which(phipred > 0)
                    phipred <- phipred[indexpred]
                    psipred <- psipred[indexpred]
                    Mprednew <- length(phipred)
                    Xpred <- matrix(0, Mprednew, N + 1)
                    
                    for (j in 1:Mprednew) {
                      
                      suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[indexpred[j], 1], N = N, delta = Tend/N, 
                        method = "milstein", theta = c(paramfixed, phipred[j], psipred[j]), model = "CIR", sigma.x = expression(psipred[j]/(2 * 
                          sqrt(x))), sigma = expression(psipred[j] * sqrt(x))))
                      
                    }
                  }
                }
                
                if (plot.pred == TRUE) {
                  op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                    cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)

                  
                  plot(timestrue, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, 
                    Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
                  for (j in indexpred) {
                    lines(timestrue, Xtrue[j, ], col = j)
                  }
                  
                  plot(times, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, 
                    max(Xtrue, Xpred) * 1.5), main = "Predicted trajectories")
                  for (j in 1:length(indexpred)) {
                    lines(times, Xpred[j, ], col = 1)
                  }
                  l.bound = level/2
                  u.bound = 1 - level/2
                  PI <- matrix(0, 2, N + 1)
                  for (k in 1:(N + 1)) {
                    PI[, k] <- quantile(Xpred[, k], c(l.bound, u.bound))
                  }
                  lines(times, PI[1, ], col = "green", lwd = 2)
                  lines(times, PI[2, ], col = "green", lwd = 2)
                }
            }
            
            
            if (sum(drift.random) == 3) {
                
                phipred <- matrix(0, 2, M - Mpred)
                
                phipred[1, ] <- rnorm(M - Mpred, res@mu[1], res@omega[1])
                
                phipred[2, ] <- rnorm(M - Mpred, res@mu[2], res@omega[2])
                
                if (model == "OU") {
                  indexpred <- which(phipred[2, ] > 0)
                  phipred <- phipred[, indexpred]
                  Mprednew <- length(indexpred)
                  Xpred <- matrix(0, Mprednew, N + 1)
                  for (j in 1:Mprednew) {
                    
                    suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[j, 1], N = N, delta = Tend/N, 
                      method = "EA", theta = c(phipred[, j], psipred[j]), model = "OU"))
                    
                  }
                }
                if (model == "CIR") {
                  indexpred <- which((phipred[1, ] > 0) & (phipred[2, ] > 0))
                  phipred <- phipred[, indexpred]
                  Mprednew <- length(indexpred)
                  Xpred <- matrix(0, Mprednew, N + 1)
                  
                  for (j in 1:Mprednew) {
                    
                    suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[indexpred[j], 1], N = N, delta = Tend/N, 
                      method = "milstein", theta = c(phipred[, j], psipred[j]), model = "CIR", sigma.x = expression(psipred[j]/(2 * 
                        sqrt(x))), sigma = expression(psipred[j] * sqrt(x))))
                    
                  }
                }
                
                if (plot.pred == TRUE) {
                  
                  op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                    cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)

                  
                  plot(times, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 
                    0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
                  for (j in indexpred) {
                    lines(times, Xtrue[j, ], col = j)
                  }
                  
                  plot(timestrue, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, 
                    max(Xtrue, Xpred) * 1.5), main = "Predicted trajectories")
                  for (j in 1:length(indexpred)) {
                    lines(timestrue, Xpred[j, ], col = 1)
                  }
                  
                  l.bound = level/2
                  u.bound = 1 - level/2
                  PI <- matrix(0, 2, N + 1)
                  for (k in 1:(N + 1)) {
                    PI[, k] <- quantile(Xpred[, k], c(l.bound, u.bound))
                  }
                  lines(timestrue, PI[1, ], col = "green", lwd = 2)
                  lines(timestrue, PI[2, ], col = "green", lwd = 2)
                }
            }
            
            
        }
      
      ############# 
      #res <- out(res)
      return(new(Class = "class.pred", estim = res, phipred = as.matrix(phipred), Xpred = Xpred, indexpred = indexpred))
      
    }
    
}







#' S4 class for the estimation results in the mixed SDE with random effects in the drift, in the diffusion or both 
#'  
#' @slot estim object of class Fit.class containing the results of the model estimation
#' @slot phipred matrix of simulated values for the random effects in the drift that are used for prediction (dimensions)
#' @slot Xpred matrix of simulated trajectories used for prediction (dimensions)
#' @slot idexpred vector of indexes of the true trajectories that are used for prediction

setClass(Class = "class.pred", representation = representation(estim = "Fit.class", phipred = "matrix", Xpred = "matrix", 
    indexpred = "numeric"))

#' S4 class for the parametric estimation results when the random effects in the drift follow
#' mixture of normal distributions
#'
#' @slot estim object of class Mixture.fit.class containing the results of the model estimation
#' @slot phipred numeric 1, 2, or c(1,2)
#' @slot Xpred matrix of predicted trajectories (dimensions)
#' @slot idexpred matrix of values on which the estimation of the density of the random effects is done

setClass(Class = "class.mixture.pred", representation = representation(estim = "Mixture.fit.class", phipred = "matrix", Xpred = "matrix",
    indexpred = "numeric"))










 
