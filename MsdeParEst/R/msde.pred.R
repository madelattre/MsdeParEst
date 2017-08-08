
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
#' @param estim.drift.fix default 0, 1 if the fixed effect in the drift is estimated, 0 if the fixed effect in the drift is known
#' @param diffusion.random default 0, 1 if one random effect in the diffusion, 0 if there is no random effect in the diffusion
#' @param diffusion.fixed default NULL, fixed effect in the diffusion: value of the fixed effect when there is no random effect in the diffusion and it is not estimated, NULL otherwise
#' @param estim.diffusion.fix default 0, 1 if the fixed effect in the diffusion is estimated, 0 otherwise
#' @param mixture 1 if the random effects in the drift follow a mixture distribution, 0 otherwise. Default to 0.
#' @param nb.mixt default 1, number of mixture components for the distribution of the random effects in the drift 
#' @param drift.fixed.mixed default 1, 1 if the value of the fixed effect in the drift is different from one mixture component to
#' another, 0 otherwise. 
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
#' @importFrom methods new
#' @importFrom methods slot
#' @importFrom methods slotNames
#' @importFrom stats kmeans
#' @importFrom stats sd
#' @importFrom graphics lines
#' @importFrom grDevices x11
#' 
#' @examples
#'
#' \dontrun{
#' # Example 1: one random effect in the drift and one random effect in the diffusion coefficient.
#' 
#' # -- Simulation
#' M <- 100
#' Tmax <- 5
#' N <- 5000
#' model <- 'OU'
#' drift.random <- 2
#' diffusion.random <- 1
#' drift.fixed <- 0
#' drift.param <- c(0.5,0.5)
#' diffusion.param <- c(8,1/2)
#'
#' sim1 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random, 
#' diffusion.random = diffusion.random, drift.fixed = drift.fixed, 
#' mixture = 0, drift.param = drift.param, diffusion.param = diffusion.param)
#'                  
#' 
#' pred1 <- pred( times = sim1$times, X = sim1$X, model = model, drift.random = drift.random,
#'               estim.drift.fix = estim.drift.fix, diffusion.random = diffusion.random, 
#'               estim.diffusion.fix = estim.diffusion.fix,level = 0.05)
#'
#' 
#' 
#' # Example 5: one fixed effect and one mixture random effect in the drift, and one fixed effect in 
#' # the diffusion coefficient
#' 
#' # -- Simulation
#' M <- 100
#' Tmax <- 5
#' N <- 5000
#' diffusion.random <- 0
#' diffusion.param <- 0.1
#' model <- 'OU'
#' drift.random <- 1
#' drift.fixed <- 1
#' nb.mixt <- 2
#' mixt.prop <- c(0.5,0.5)
#' param.ea1 <- c(0.5, 0.25, 1.8, 0.25)
#' param.ea2 <- c(1, 0.25, 1, 0.25) 
#' drift.param <- param.ea1
#' 
#' sim5 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random,
#'                  diffusion.random = diffusion.random, drift.fixed = drift.fixed,
#'                  mixture=1, drift.param = drift.param,
#'                  diffusion.param = diffusion.param, nb.mixt = nb.mixt, mixt.prop = mixt.prop)
#'
#' 
#' res5 <- msde.pred(times = sim5$times, X = sim5$X, model = model, drift.random = drift.random,
#'                   estim.drift.fix = estim.drift.fix, diffusion.random = diffusion.random, estim.diffusion.fix = estim.diffusion.fix,
#'                   mixture = mixture, nb.mixt=nb.mixt, Niter = Niter)
#' 
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


msde.pred <- function(times, X, model = c("OU", "CIR"), drift.random, drift.fixed = NULL, estim.drift.fix = 0, 
    diffusion.random = 0, diffusion.fixed = NULL, estim.diffusion.fix = 0, mixture = 0, nb.mixt = 1, drift.fixed.mixed = 0, 
    Niter = 10, discrete = 1, plot.pred = TRUE, level = 0.05, newwindow = FALSE) {
    model <- match.arg(model)
    
    if (newwindow) {
        x11(width = 10)
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
    Xtrue <- X[(Mpred + 1):M, ]
    
    Xestim <- X[1:Mpred, ]
    
    timestrue <- times
    Tend <- timestrue[length(timestrue)]
    delta <- round(diff(timestrue), 10)[1]
    
    times <- seq(0, Tend, by = delta)
    
    res <- msde.fit(times = times, X = Xestim, model = model, drift.random = drift.random, estim.drift.fix = estim.drift.fix, 
        diffusion.random = diffusion.random, estim.diffusion.fix = estim.diffusion.fix, mixture = mixture, nb.mixt = nb.mixt, 
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
                paramfixed <- res@mu[Niter, 1, 2]
                
                phipred <- rep(0, M - Mpred)
                
                phipred <- mixture.sim(M - Mpred, c(res@mu[Niter, 1, 1], res@omega[Niter, 2, 1], res@mu[Niter, 
                  2, 1], res@omega[Niter, 1, 1]), res@mixt.prop[10, ])
                # A VERIFIER
                
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
                  Mpred <- length(phipred)
                  Xpred <- matrix(0, M - Mpred, N + 1)
                  
                  for (j in 1:(M - Mpred)) {
                    
                    suppressMessages(Xpred[j, ] <- sde.sim(T = Tend, X0 = Xtrue[indexpred[j], 1], N = N, delta = Tend/N, 
                      method = "milstein", theta = c(phipred[j], paramfixed, sig), model = "CIR", sigma.x = expression(sig/(2 * 
                        sqrt(x))), sigma = expression(sig * sqrt(x))))
                    
                  }
                }
                phipred <- matrix(phipred, 1, length(phipred))
            }
            
            if (drift.random == 2) {
                paramfixed <- res@mu[Niter, 1, 1]
                
                phipred <- rep(0, M - Mpred)
                phipred <- mixture.sim(M - Mpred, c(res@mu[Niter, 1, 2], res@omega[Niter, 2, 2], res@mu[Niter, 
                  2, 2], res@omega[Niter, 1, 2]), res@mixt.prop[10, ])
                
                
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
                phipred <- matrix(phipred, 1, length(phipred))
            }
            
            if (plot.pred == TRUE) {
                op <- par(mfrow = c(1, 3), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                  cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
                
                
                plot(sort(res@estimphi[indexpred]), sort(phipred), pch = 18, xlim = c(min(res@estimphi[indexpred], 
                  phipred) * 0.8, max(res@estimphi[indexpred], phipred) * 1.2), ylim = c(min(res@estimphi[indexpred], 
                  phipred) * 0.8, max(res@estimphi[indexpred], phipred) * 1.2), ylab = "", xlab = "", main = "Random effect in the drift")
                abline(0, 1)
                
                
                plot(timestrue, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 
                  0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
                for (j in indexpred) {
                  lines(timestrue, Xtrue[j, ], col = j)
                }
                
                plot(times, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, max(Xtrue, 
                  Xpred) * 1.5), main = "Predictive trajectories")
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
        
        
        if (drift.random == c(1, 2)) {
            
            phipred <- matrix(0, 2, M - Mpred)
            
            phipred[1, ] <- mixture.sim(M - Mpred, c(res@mu[Niter, 1, 1], res@omega[Niter, 2, 1], res@mu[Niter, 
                2, 1], res@omega[Niter, 1, 1]), res@mixt.prop[10, ])
            
            phipred[2, ] <- mixture.sim(M - Mpred, c(res@mu[Niter, 1, 2], res@omega[Niter, 2, 2], res@mu[Niter, 
                2, 2], res@omega[Niter, 1, 2]), res@mixt.prop[10, ])
            
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
                
                op <- par(mfrow = c(2, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                  cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
                
                
                plot(sort(res@estimphi[1, indexpred]), sort(phipred[1, ]), pch = 18, xlim = c(min(res@estimphi[1, 
                  ], phipred[1, ]) * 0.8, max(res@estimphi[1, ], phipred[1, ]) * 1.2), ylim = c(min(res@estimphi[1, 
                  ], phipred[1, ]) * 0.8, max(res@estimphi[1, ], phipred[1, ]) * 1.2), ylab = "", xlab = "", main = "First random effect in the drift")
                abline(0, 1)
                plot(sort(res@estimphi[2, indexpred]), sort(phipred[2, ]), pch = 18, xlim = c(min(res@estimphi[2, 
                  ], phipred[2, ]) * 0.8, max(res@estimphi[2, ], phipred[2, ]) * 1.2), ylim = c(min(res@estimphi[2, 
                  ], phipred[2, ]) * 0.8, max(res@estimphi[2, ], phipred[2, ]) * 1.2), ylab = "", xlab = "", main = "Second random effect in the drift")
                abline(0, 1)
                
                
                plot(times, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 
                  0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
                for (j in indexpred) {
                  lines(times, Xtrue[j, ], col = j)
                }
                
                plot(timestrue, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, 
                  max(Xtrue, Xpred) * 1.5), main = "Predictive trajectories")
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
                  op <- par(mfrow = c(1, 3), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                    cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
                  
                  ## MODIF SUR x@estimphi
                  plot(sort(res@estimphi[indexpred]), sort(phipred), pch = 18, xlim = c(min(res@estimphi[indexpred], 
                    phipred) * 0.8, max(res@estimphi[indexpred], phipred) * 1.2), ylim = c(min(res@estimphi[indexpred], 
                    phipred) * 0.8, max(res@estimphi[indexpred], phipred) * 1.2), ylab = "", xlab = "", main = "Random effect in the drift")
                  abline(0, 1)
                  ## FIN MODIF SUR x@estimphi
                  
                  plot(timestrue, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, 
                    Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
                  for (j in indexpred) {
                    lines(timestrue, Xtrue[j, ], col = j)
                  }
                  
                  plot(times, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, 
                    max(Xtrue, Xpred) * 1.5), main = "Predictive trajectories")
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
            if (drift.random == c(1, 2)) {
                
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
                  
                  op <- par(mfrow = c(2, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                    cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
                  
                  
                  plot(sort(res@estimphi[1, indexpred]), sort(phipred[1, ]), pch = 18, xlim = c(min(res@estimphi[1, 
                    ], phipred[1, ]) * 0.8, max(res@estimphi[1, ], phipred[1, ]) * 1.2), ylim = c(min(res@estimphi[1, 
                    ], phipred[1, ]) * 0.8, max(res@estimphi[1, ], phipred[1, ]) * 1.2), ylab = "", xlab = "", 
                    main = "First random effect in the drift")
                  abline(0, 1)
                  plot(sort(res@estimphi[2, indexpred]), sort(phipred[2, ]), pch = 18, xlim = c(min(res@estimphi[2, 
                    ], phipred[2, ]) * 0.8, max(res@estimphi[2, ], phipred[2, ]) * 1.2), ylim = c(min(res@estimphi[2, 
                    ], phipred[2, ]) * 0.8, max(res@estimphi[2, ], phipred[2, ]) * 1.2), ylab = "", xlab = "", 
                    main = "Second random effect in the drift")
                  abline(0, 1)
                  
                  
                  plot(times, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 
                    0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
                  for (j in indexpred) {
                    lines(times, Xtrue[j, ], col = j)
                  }
                  
                  plot(timestrue, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, 
                    max(Xtrue, Xpred) * 1.5), main = "Predictive trajectories")
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
            
            if (sum(drift.random != 3)) {
                if (sum(drift.random) == 1) {
                  
                  paramfixed <- res@mu[2]
                  
                  phipred <- rep(0, M)
                  
                  for (j in 1:M) {
                    phipred[j] <- rnorm(1, mean = res@mu[1], sd = res@omega[1] * psipred[j])
                  }
                  
                  
                  if (model == "OU") {
                    indexpred <- 1:M
                    for (j in 1:M) {
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
                  
                  phipred <- rep(0, M)
                  
                  for (j in 1:M) {
                    phipred <- rnorm(M, res@mu[2], res@omega[2] * psipred[j])
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
                  op <- par(mfrow = c(2, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                    cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
                  
                  plot(sort(res@estimphi[indexpred]), sort(phipred), pch = 18, xlim = c(min(res@estimphi[indexpred], 
                    phipred) * 0.8, max(res@estimphi[indexpred], phipred) * 1.2), ylim = c(min(res@estimphi[indexpred], 
                    phipred) * 0.8, max(res@estimphi[indexpred], phipred) * 1.2), ylab = "", xlab = "", main = "Random effect in the drift")
                  abline(0, 1)
                  
                  plot(sort(sqrt(res@estimpsi2[indexpred])), sort(psipred), pch = 18, xlim = c(min(sqrt(res@estimpsi2[indexpred]), 
                    psipred) * 0.8, max(sqrt(rse@estimpsi2[indexpred]), psipred) * 1.2), ylim = c(min(sqrt(res@estimpsi2[indexpred]), 
                    psipred) * 0.8, max(sqrt(res@estimpsi2[indexpred]), psipred) * 1.2), ylab = "", xlab = "", 
                    main = "Random effect in the diffusion")
                  abline(0, 1)
                  
                  plot(timestrue, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, 
                    Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
                  for (j in indexpred) {
                    lines(timestrue, Xtrue[j, ], col = j)
                  }
                  
                  plot(times, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, 
                    max(Xtrue, Xpred) * 1.5), main = "Predictive trajectories")
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
            
            
            if (drift.random == c(1, 2)) {
                
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
                  
                  op <- par(mfrow = c(2, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                    cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
                  
                  
                  plot(sort(res@estimphi[1, indexpred]), sort(phipred[1, ]), pch = 18, xlim = c(min(res@estimphi[1, 
                    ], phipred[1, ]) * 0.8, max(res@estimphi[1, ], phipred[1, ]) * 1.2), ylim = c(min(res@estimphi[1, 
                    ], phipred[1, ]) * 0.8, max(res@estimphi[1, ], phipred[1, ]) * 1.2), ylab = "", xlab = "", 
                    main = "First random effect in the drift")
                  abline(0, 1)
                  plot(sort(res@estimphi[2, indexpred]), sort(phipred[2, ]), pch = 18, xlim = c(min(res@estimphi[2, 
                    ], phipred[2, ]) * 0.8, max(res@estimphi[2, ], phipred[2, ]) * 1.2), ylim = c(min(res@estimphi[2, 
                    ], phipred[2, ]) * 0.8, max(res@estimphi[2, ], phipred[2, ]) * 1.2), ylab = "", xlab = "", 
                    main = "Second random effect in the drift")
                  abline(0, 1)
                  
                  
                  plot(times, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 
                    0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
                  for (j in indexpred) {
                    lines(times, Xtrue[j, ], col = j)
                  }
                  
                  plot(timestrue, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 0.8, 
                    max(Xtrue, Xpred) * 1.5), main = "Predictive trajectories")
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
    }
    
    ############# 
    res <- out(res)
    return(new(Class = "class.mixture.pred", res = res, phipred = phipred, Xpred = Xpred, indexpred = indexpred))
}







#' S4 class for the estimation results in the mixed SDE with random effects in the drift, in the diffusion or both 
#'  
#' @slot res character 'OU' or 'CIR'
#' @slot phipred numeric 1, 2, or c(1,2)
#' @slot Xpred
#' @slot idexpred matrix of values on which the estimation of the density of the random effects is done

setClass(Class = "class.pred", representation = representation(res = "list", phipred = "matrix", Xpred = "matrix", 
    indexpred = "numeric"))

#' S4 class for the parametric estimation results when the random effects in the drift follow 
#' mixture of normal distributions  
#'  
#' @slot res character 'OU' or 'CIR'
#' @slot phipred numeric 1, 2, or c(1,2)
#' @slot Xpred
#' @slot idexpred matrix of values on which the estimation of the density of the random effects is done

setClass(Class = "class.mixture.pred", representation = representation(res = "list", phipred = "matrix", Xpred = "matrix", 
    indexpred = "numeric"))










 
