# MsdeParEst R package ; file msde.fit.r (last modified: 2017-08-09)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017
# UMR 518 AgroParisTech/INRA, 16 rue Claude Bernard, 75 231 Paris Cedex 05


#' Estimation Of The Random Effects In Mixed Stochastic Differential Equations
#' 
#' 
#' @description Parametric estimation of the joint density of the random effects \eqn{(\alpha_j, \beta_j, \sigma_j)} in the mixed SDE
#' 
#'  \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j a(X_j(t)) dW_j(t)}.
#'  
#' @param times vector of observation times
#' @param X matrix of the M trajectories (each row is a trajectory with as much columns as observations)
#' @param model name of the SDE: 'OU' (Ornstein-Uhlenbeck) or 'CIR' (Cox-Ingersoll-Ross)
#' @param drift.random random effects in the drift: 0 if only fixed effects, 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects. Defaults to c(1,2).
#' @param drift.fixed NULL if the fixed effect(s) in the drift is (are) estimated, value of the fixed effect(s) otherwise. Default to NULL
#' @param diffusion.random default 0, 1 if one random effect in the diffusion, 0 if there is no random effect in the diffusion
#' @param diffusion.fixed NULL if the fixed effect in the diffusion is estimated, value of the fixed effect otherwise. Default to NULL
#' @param mixture 1 if the random effects in the drift follow a mixture distribution, 0 otherwise. Default to 0.
#' @param nb.mixt default 1, number of mixture components for the distribution of the random effects in the drift otherwise. 
#' @param Niter default 10, number of iterations for the EM algorithm if mixture = 1
#' @param discrete default 1, 1 for discrete observations, 0 otherwise. If discrete = 0, and diffusion.random = 0, the exact likelihood associated with continuous observations is 
#' discretized. If discrete = 1, the likelihood of the Euler scheme of the mixed SDE is computed. 
#' @return 
#'

#' \item{index}{is the vector of subscript in \eqn{1,...,M} where the estimation of \eqn{phi} has been done,  most of the time \eqn{index= 1:M}}
#' \item{estimphi}{matrix of estimators of \eqn{\phi=\alpha, or \beta, or (\alpha,\beta)} from the efficient statitics (see \code{\link{UVS}}), matrix of two lines if drift.random =c(1,2), numerical type otherwise}
#' \item{estimpsi2}{matrix of estimators of \eqn{\psi^2=\sigma^2} from the efficient statistics (see \code{\link{UVS}}), matrix of one line}
#' \item{gridf}{grid of values for the plots of the random effects distribution in the drift, matrix form}
#' \item{gridg}{grid of values for the plots of the random effects distribution in the diffusion, matrix form}
#' \item{estimf}{estimator of the density of \eqn{\phi} from a kernel estimator from package: stats, function: density. Matrix form: one line if one random effect or square matrix otherwise}
#' \item{estimg}{estimator of the density of \eqn{\psi^2}. Matrix form: one line if one random effect or square matrix otherwise}
#' \item{mu}{estimator of the mean of the random effects normal density}
#' \item{omega}{estimator of the standard deviation of the random effects normal density}
#' \item{a}{estimated value of the shape of the Gamma distribution}
#' \item{lambda}{estimated value of the scale of the Gamma distribution}
#' \item{sigma2}{value of the diffusion coefficient if it is fixed}
#' \item{bic}{BIC criterium}
#' \item{aic}{AIC criterium}
#' \item{model}{initial choice}
#' \item{drift.random}{initial choice}
#' \item{diffusion.random}{initial choice}
#' \item{drift.fixed}{initial choice}
#' \item{estim.drift.fix}{initial choice}
#' \item{estim.diffusion.fixed}{initial choice}
#' \item{discrete}{initial choice}
#' \item{times}{initial choice}
#' \item{X}{initial choice}
#' 
#' For the 'paramMLmixture' method:
#' \item{mu}{estimated value of the mean at each iteration of the algorithm. Niter x N x 2 array. }
#' \item{omega}{estimated value of the standard deviation at each iteration of the algorithm. Niter x N x 2 array.}
#' \item{mixt.prop}{estimated value of the mixture proportions at each iteration of the algorithm. Niter x N matrix.}
#' \item{probindi}{posterior component probabilites. M x N matrix.}
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
#' 
#' @details
#' Estimation of the random effects density from M independent trajectories of the SDE (the Brownian motions \eqn{W_j} are independent), with linear drift. 
#' The drift includes no, one or two random effects: 
#' 
#' if drift.random = 0: \eqn{\alpha_j \equiv \alpha} and \eqn{\beta_j \equiv \beta} are fixed
#' 
#' if drift.random = 1: \eqn{\beta_j \equiv \beta} is fixed and \eqn{\alpha_j} is random
#' 
#' if drift.random = 2: \eqn{\alpha_j \equiv \alpha} is fixed and \eqn{\beta_j} is random  
#' 
#' if drift.random = c(1,2): \eqn{\alpha_j} and \eqn{\beta_j} are random
#' 
#' The diffusion includes either a fixed effect or a random effect:
#' 
#' if diffusion.random = 0: \eqn{\sigma_j \equiv \sigma} is fixed
#' 
#' if diffusion.random = 1: \eqn{\sigma_j} is random
#' 
#' If there is no random effect in the diffusion (diffusion.random = 0), the drift random effect follow Gaussian distributions: 
#' \eqn{\alpha_j,\beta_j \sim N(\mu,\Omega)}.
#' If there is one random effect (\eqn{\sigma_j}) in the diffusion (diffusion.random = 1), \eqn{\sigma_j \sim Gamma(a,\lambda)}, and 
#' \eqn{\alpha_j,\beta_j|\sigma_j \sim N(\mu,\sigma_j^2 \Omega)}.
#' 
#' Two diffusions are implemented: 
#' 
#' the Ornstein-Uhlenbeck model (OU) \eqn{a(X_j(t))=1}
#' 
#' the Cox-Ingersoll-Ross model (CIR) \eqn{a(X_j(t))=\sqrt{X_j(t)}}
#' 
#'  Validation method:
#'  For a number of trajectory numj (fixed by the user or randomly chosen) this function simulates 
#'  Mrep =100 (by default) new trajectories with the value of the estimated random effect. 
#'  Then it plots on the left graph the Mrep new trajectories 
#'  \eqn{(Xnumj^{k}(t1), ... Xnumj^{k}(tN)), k= 1, ... Mrep} with in red the true trajectory 
#'  \eqn{(Xnumj(t1), ... Xnumj(tN))}. The right graph is a qq-plot of the quantiles of samples 
#'  \eqn{(Xnumj^{1}(ti), ... Xnumj^{Mrep}(ti))}
#'  for each time \eqn{ti} compared with the uniform quantiles. The outputs of the function  
#'  are: a matrix \code{Xnew} dimension Mrepx N+1, vector of quantiles \code{quantiles} length 
#'  N and the number of the trajectory for the plot \code{numj} 
#' 
#'  Prediction method: (A COMPLETER)
#'  
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
#' # -- Estimation
#' 
#' # -----Fixed effect in the drift estimated
#' res1 <- msde.fit(times = sim1$times, X = sim1$X, model = 'OU', drift.random = 2, 
#'                  diffusion.random = 1, estim.drift.fix = 1, mixture = 0)
#' summary(res1)
#' valid(res1)
#' plot(res1)
#' 
#' # ----- Fixed effect in the drift known and not estimated
#' res1bis <- msde.fit(times = sim1$times, X = sim1$X, model = 'OU', drift.random = 2, 
#'                     diffusion.random = 1, drift.fixed=0, mixture = 0)
#' summary(res1bis)
#' 
#' # Example 2: one random effect in the drift and one fixed effect in the diffusion coefficient
#' 
#' # -- Simulation
#' M <- 100
#' Tmax <- 5
#' N <- 5000
#' model <- 'OU'
#' diffusion.random <- 0
#' diffusion.param <- 0.5
#' drift.random <- 2
#' drift.fixed <- 10
#' drift.param <- c(1,sqrt(0.4/4))
#' 
#' sim2 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random,
#' diffusion.random = diffusion.random, drift.fixed = drift.fixed,
#' mixture=0, drift.param = drift.param,
#' diffusion.param = diffusion.param)
#' 
#' # -- Estimation
#' res2 <- msde.fit(times = sim2$times, X = sim2$X, model = 'OU', drift.random = 2, 
#'                  diffusion.random = 0, estim.drift.fix = 1, mixture = 0)
#' 
#' summary(res2)
#' plot(res2)
#' 
#' # Example 3: two random effects in the drift and one random effect in the diffusion coefficient
#' 
#' # -- Simulation
#' M <- 100
#' Tmax <- 5
#' N <- 5000
#' model <- 'OU'
#' drift.random <- c(1,2)
#' diffusion.random <- 1
#' density.phi <- 'normalnormal'
#' drift.param <- c(1,0.5,0.5,0.5)
#' diffusion.param <- c(8,1/2)
#' 
#' sim3 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random,
#' diffusion.random = diffusion.random, mixture = 0,
#' drift.param = drift.param, diffusion.param = diffusion.param)
#' 
#' # -- Estimation
#' 
#' res3 <- msde.fit(times = sim3$times, X = sim3$X, model = 'OU', drift.random = c(1,2), 
#'                  diffusion.random = 1, mixture = 0)
#' summary(res3)
#' plot(res3)
#' 
#' # Example 4: fixed effects in the drift and one random effect in the diffusion coefficient
#' 
#' # -- Simulation
#' M <- 100
#' Tmax <- 5
#' N <- 5000
#' model <- 'OU'
#' drift.random <- 0
#' diffusion.random <- 1
#' drift.fixed <- c(0,1)
#' diffusion.param <- c(5,3)
#' 
#' sim4 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random,
#' diffusion.random = diffusion.random, drift.fixed = drift.fixed,
#' diffusion.param = diffusion.param)
#'
#' # -- Estimation
#' res4 <- msde.fit(times = sim4$times, X = sim4$X, model = 'OU', drift.random = 0, 
#'                  diffusion.random = 1, mixture = 0, estim.drift.fix = 0, 
#'                  drift.fixed = c(0,0), discrete = 1)
#' 
#' summary(res4)
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
#' # -- Estimation
#' res5 <- msde.fit(times = sim5$times, X = sim5$X, model = 'OU', drift.random = 1, 
#'                  estim.drift.fix = 1, diffusion.random = 0, estim.diffusion.fix = 1, 
#'                  mixture = 1, nb.mixt=2, Niter = 25)
#' 
#' summary(res5)
#' plot(res5)
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


msde.fit <- function(times, X, model = c("OU", "CIR"), drift.random = c(1,2), drift.fixed = NULL, 
    diffusion.random = 0, diffusion.fixed = NULL, mixture = 0, nb.mixt = 1,  
    Niter = 10, discrete = 1) {
    model <- match.arg(model)
    
    if((model != 'OU')&(model != 'CIR')){stop("A model must be precised: OU or CIR")}
    
    if (is.matrix(X)) {
        if (nrow(X) == length(times)) {
            X <- t(X)
        } else {
            if (ncol(X) != length(times)) {
                stop("Length of times has to be equal to the columns of X")
            }
        }
    }
    
    M <- dim(X)[1]
    K <- dim(X)[2]
    delta <- round(diff(times), 10)
    
    Tend <- times[length(times)]
    
    if (mixture == 0) {
        
        a <- 0
        lambda <- 0
        mu <- 0
        omega <- 0
        sigma2 <- 0
        
        gridf <- as.matrix(0)
        estimf <- as.matrix(0)
        
        estimg <- as.matrix(0)
        gridg <- as.matrix(0)

        estimphi <- as.matrix(0)
        
        if (model == "OU") {
            Mindex <- M
            index <- 1:Mindex
        }
        
        if (model == "CIR") {
            
            index <- which(rowSums(X <= 0) == 0)
            Mindex <- length(index)
            if (Mindex == 0) {
                warning("All the trajectories have non positive values the model CIR cannot be used", call. = FALSE)
                estimf <- 0
                estimphi <- 0
                estimpsi2 <- 0
                bic <- 0
                aic <- 0
                gridf <- 0
                mu <- 0
                omega <- 0
                sigma2 <- 0
                a <- 0
                lambda <- 0
          }
            
        }
        
        
        if (Mindex > 0) {
            
            # -- Errors and warnings
            
            if (((sum(drift.random)) == 0) && (diffusion.random == 0)) {
                stop("There should be at least one random effect either in the drift or the diffusion coefficient.")
            }
            
            
            # -- Computation of the sufficient statistics
            
            U <- matrix(0, 2, Mindex)
            V <- as.list(1:Mindex)
            S <- rep(0, Mindex)
            SigDelta <- rep(0, Mindex)
            
            estimUV <- UVS(X[index, ], model, times)
            
            U <- estimUV$U
            V <- estimUV$V
            S <- estimUV$S
            SigDelta <- estimUV$SigDelta
            
            deter <- lapply(V, det)
            index2 <- which((deter != Inf) & (deter != 0))  # indexes in 1:Mindex 
            Mindex2 <- length(index2)
            V <- V[index2]
            U <- U[, index2]
            S <- S[index2]
            SigDelta <- SigDelta[index2]
            
            # -- Estimation of the random effects
            
            estimphi <- matrix(0, 2, Mindex2)
            
            for (j in 1:Mindex2) {
                estimphi[, j] <- solve(V[[j]]) %*% U[, j]
            }
            
            estimpsi2 <- rep(0, Mindex2)
            estimpsi2 <- S/K
            
            
            # -- Parameter estimation
            
            
            if (((sum(drift.random)) >= 1) && (diffusion.random == 0)) {
                
                res <- EstParamNormal(U, V, S, SigDelta, K, drift.fixed, diffusion.fixed, drift.random, 
                  discrete)
                
                bic <- res$BIChere
                aic <- res$AIChere
                mu <- res$mu
                omega <- res$omega
                sigma2 <- res$sigma^2

                # -- Estimator of the density
                
                estimpsi2 <- matrix(estimpsi2, 1, length(estimpsi2), byrow = TRUE)
                
                if (sum(drift.random) == 3) {
                  
                    gridf <- matrix(0, 2, 500)
                    gridf[1, ] <- seq(mu[1] - 3 * omega[1], mu[1] + 3 * omega[1], length = 500)
                    gridf[2, ] <- seq(mu[2] - 3 * omega[2], mu[2] + 3 * omega[2], length = 500)

                  
                  estimf1 <- dnorm(gridf[1, ], mean = mu[1], sd = abs(omega[1]))
                  estimf2 <- dnorm(gridf[2, ], mean = mu[2], sd = abs(omega[2]))
                  
                  estimf <- rbind(estimf1, estimf2)

                }
                
                if (sum(drift.random) == 2) {
                  
                  estimphi <- estimphi[2, ]
                  

                    gridf <- seq(mu[2] - 3 * omega[2], mu[2] + 3 * omega[2], length = 500)

                  
                  estimf <- matrix(dnorm(gridf, mean = mu[2], sd = omega[2]), 1, length(gridf), byrow = TRUE)
                  gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
                  estimphi <- matrix(estimphi, 1, length(estimphi), byrow = TRUE)

                }
                
                if (sum(drift.random) == 1) {
                  
                  estimphi <- estimphi[1, ]
                  
                    gridf <- seq(mu[1] - 3 * omega[1], mu[1] + 3 * omega[1], length = 500)

                  estimf <- matrix(dnorm(gridf, mean = mu[1], sd = omega[1]), 1, length(gridf), byrow = TRUE)
                  gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
                  estimphi <- matrix(estimphi, 1, length(estimphi), byrow = TRUE)

                }
                
            }
            
            if ((sum(drift.random) == 0) && (diffusion.random == 1)) {
                
                if (discrete == 0) {
                  warning("Estimation cannot be performed from the likelihood associated with a continuous observation of the trajectories.")
                }
                
                res <- EstParamGamma(U, V, S, SigDelta, K, drift.fixed)
                
                bic <- res$BIChere
                aic <- res$AIChere
                mu <- res$mu
                a <- res$a
                lambda <- res$lambda
                
                
                # Estimator of the density

                gridg <- seq(min(estimpsi2) * 0.8, max(estimpsi2) * 1.2, length = 500)
           
                simupsi2 <- 1/rgamma(500, a, rate = 1/lambda)
                
                testpsi <- suppressMessages(density(simupsi2, from = min(simupsi2), to = max(simupsi2), bw = "ucv", 
                  n = 500))
                if (testpsi$bw < 0.1) {
                  testpsi <- suppressMessages(density(simupsi2, from = min(simupsi2), to = max(simupsi2), n = 500))
                }
                
                gridg <- seq(min(estimpsi2) * 0.8, max(estimpsi2) * 1.2, length = 500)
                estimg <- matrix(testpsi$y, nrow = 1)
                gridg <- matrix(gridg, nrow = 1)
                
                estimpsi2 <- matrix(estimpsi2, 1, length(estimpsi2), byrow = TRUE)

            }
            
            if (((sum(drift.random)) >= 1) && (diffusion.random == 1)) {
                
                # -- L'estimation de la densité va nécessiter de simuler des réalisation des effets aléatoires dans la dérive
                
                res <- EstParamNormalGamma(U, V, S, SigDelta, K, drift.random, drift.fixed)
                bic <- res$BIChere
                aic <- res$AIChere
                mu <- res$mu
                omega <- res$omega
                a <- res$a
                lambda <- res$lambda
                
                # -- Estimator of the density
                
                estimpsi2 <- matrix(estimpsi2, 1, length(estimpsi2), byrow = TRUE)
                
                simupsi2 <- 1/rgamma(500, a, rate = 1/lambda)
                testpsi <- suppressMessages(density(simupsi2, from = min(simupsi2), to = max(simupsi2), bw = "ucv", 
                  n = 500))
                if (testpsi$bw < 0.1) {
                  testpsi <- suppressMessages(density(simupsi2, from = min(simupsi2), to = max(simupsi2), n = 500))
                }
                
                gridg <- seq(min(estimpsi2) * 0.8, max(estimpsi2) * 1.2, length = 500)
                estimg <- matrix(testpsi$y, nrow = 1)
                gridg <- matrix(gridg, nrow = 1)
                
                
                if (sum(drift.random) == 3) {
                  
                  simuphi <- matrix(NA, 2, length(simupsi2))
                  simuphi[1, ] <- sapply(1:length(simupsi2), function(s) {
                    rnorm(1, mean = mu[1], sd = abs(omega[1] * sqrt(simupsi2[s])))
                  })
                  simuphi[2, ] <- sapply(1:length(simupsi2), function(s) {
                    rnorm(1, mean = mu[2], sd = abs(omega[2] * sqrt(simupsi2[s])))
                  })
                  
                  test1 <- suppressMessages(density(simuphi[1, ], from = min(simuphi[1, ]), to = max(simuphi[1, 
                    ]), bw = "ucv", n = 500))
                  test2 <- suppressMessages(density(simuphi[2, ], from = min(simuphi[2, ]), to = max(simuphi[2, 
                    ]), bw = "ucv", n = 500))
                  
                  if (test1$bw < 0.1) {
                    test1 <- suppressMessages(density(simuphi[1, ], from = min(simuphi[1, ]), to = max(simuphi[1, 
                      ]), n = 500))
                  }
                  estimf1 <- test1$y
                  gridf1 <- test1$x
                  
                  if (test2$bw < 0.1) {
                    test2 <- suppressMessages(density(simuphi[2, ], from = min(simuphi[2, ]), to = max(simuphi[2, 
                      ]), n = 500))
                  }
                  estimf2 <- test2$y
                  gridf2 <- test2$x
                  
                  gridf <- matrix(0, 2, 500)
                  gridf[1, ] <- gridf1
                  gridf[2, ] <- gridf2
                  
                  estimf <- matrix(0, 2, 500)
                  estimf[1, ] <- estimf1
                  estimf[2, ] <- estimf2

                }
                
                if (sum(drift.random) == 2) {
                  simuphi <- rep(NA, length(simupsi2))
                  simuphi <- sapply(1:length(simupsi2), function(s) {
                    rnorm(1, mean = mu[2], sd = abs(omega[2] * sqrt(simupsi2[s])))
                  })
                  
                  
                  test <- suppressMessages(density(simuphi, from = min(simuphi), to = max(simuphi), bw = "ucv", 
                    n = 500))
                  if (test$bw < 0.1) {
                    test <- suppressMessages(density(simuphi, from = min(simuphi), to = max(simuphi), n = 500))
                  }
                  
                  gridf <- test$x
                  estimf <- test$y
                  
                  estimf <- matrix(estimf, nrow = 1)
                  gridf <- matrix(gridf, nrow = 1)
                  
                }
                
                if (sum(drift.random) == 1) {
                  
                  simuphi <- rep(NA, length(simupsi2))
                  simuphi <- sapply(1:length(simupsi2), function(s) {
                    rnorm(1, mean = mu[1], sd = abs(omega[1] * sqrt(simupsi2[s])))
                  })
                  
                  test <- suppressMessages(density(simuphi, from = min(simuphi), to = max(simuphi), n = length(simuphi), 
                    bw = "ucv"))
                  
                  if (test$bw < 0.1) {
                    test <- suppressMessages(density(simuphi, from = min(simuphi), to = max(simuphi), n = length(simuphi)))
                  }
                  
                  
                  gridf <- test$x
                  estimf <- test$y
                  
                  
                  
                  estimf <- matrix(estimf, nrow = 1)
                  gridf <- matrix(gridf, nrow = 1)
                  
                
                }
                
            }
            
        }
        
        
        return(new(Class = "Fit.class", model = model, drift.random = drift.random, diffusion.random = diffusion.random, 
            gridf = gridf, mu = mu, omega = omega, a = a, lambda = lambda, sigma2 = sigma2, index = index, 
            estimphi = estimphi, estimpsi2 = estimpsi2, estimf = estimf, estim.drift.fix = as.numeric(is.null(drift.fixed)), 
            estim.diffusion.fix = as.numeric(is.null(diffusion.fixed)), 
            discrete = discrete, bic = bic, aic = aic, times = times, X = X, gridg = gridg, estimg = estimg))
    }
    
    
    if (mixture == 1) {
        
        gridf <- NULL
        estimf <- NULL

        
        
        if (diffusion.random == 1) {
            warning("For the considered mixtures of SDE, there should not be random effects in the diffusion coefficient. diffusion.random is set to 0.")
            diffusion.random <- 0
        }
        
        if (sum(drift.random) == 0) {
            stop("There should be at least one random effect in the drift.")
        }
        
        if (model == "OU") {
            Mindex <- M
            index <- 1:Mindex
        }
        
        if (model == "CIR") {
            
            index <- which(rowSums(X <= 0) == 0)
            Mindex <- length(index)
            if (Mindex == 0) {
                warning("All the trajectories have non positive values the model CIR cannot be used", call. = FALSE)
                estimf <- 0
                estimphi <- 0
                estimpsi2 <- 0
                bic <- 0
                aic <- 0
                gridf <- 0
                mu <- 0
                omega <- 0
                sigma2 <- 0
                a <- 0
                lambda <- 0
            }
            
        }
        
        # -- Computation of the sufficient statistics
        
        U <- matrix(0, 2, Mindex)
        V <- as.list(1:Mindex)
        S <- rep(0, Mindex)
        SigDelta <- rep(0, Mindex)
        
        estimUV <- UVS(X[index, ], model, times)
        
        U <- estimUV$U
        V <- estimUV$V
        S <- estimUV$S
        SigDelta <- estimUV$SigDelta
        
        deter <- lapply(V, det)
        index2 <- which((deter != Inf) & (deter != 0))  # indexes in 1:Mindex 
        Mindex2 <- length(index2)
        V <- V[index2]
        U <- U[, index2]
        S <- S[index2]
        SigDelta <- SigDelta[index2]
        
        # -- Initialization of the EM algorithm
        
        estimphi <- matrix(0, 2, Mindex2)
        
        for (j in 1:Mindex2) {
            estimphi[, j] <- solve(V[[j]]) %*% U[, j]
        }
        
        muinit <- matrix(0, nb.mixt, 2)
        omegainit <- matrix(0, nb.mixt, 2)
        
        km <- kmeans(t(estimphi), 2)
        
        if (sum(drift.random) == 3) {
            for (n in 1:nb.mixt) {
                label <- which(km$cluster == n)
                muinit[n, ] <- apply(estimphi[, label], 1, mean)
                omegainit[n, ] <- apply(estimphi[, label], 1, sd)
            }
        }
        
        if (sum(drift.random) == 2) {
            for (n in 1:nb.mixt) {
                label <- which(km$cluster == n)
                muinit[n, 2] <- mean(estimphi[2, label])
                omegainit[n, 2] <- sd(estimphi[2, label])
            }
            muinit[, 1] <- mean(estimphi[1, ])
            
        }
        
        if (sum(drift.random) == 1) {
            for (n in 1:nb.mixt) {
                label <- which(km$cluster == n)
                muinit[n, 1] <- mean(estimphi[1, label])
                omegainit[n, 1] <- sd(estimphi[1, label])
            }
            muinit[, 2] <- mean(estimphi[2, ])
        }
        
        probinit <- table(km$cluster)/Mindex2
        
        start <- list(mu = muinit, omega = omegainit, mixt.prop = probinit)
        
        # -- Parameter estimation
        
        res <- EM(U, V, S, K, drift.random, start, Niter, drift.fixed, diffusion.fixed)
        
        bic <- res$BIChere
        aic <- res$AIChere
        mu <- res$mu
        omega <- res$omega
        sigma2 <- res$sigma^2
        mixt.prop <- res$mixt.prop
        probindi <- res$probindi
        
        # -- Estimator of the density
        
        if (sum(drift.random) == 3) {
            
                gridf <- matrix(0, 2, 500)
                bg1 <- max(abs(estimphi[1, ])) * 1.2
                bg2 <- max(abs(estimphi[2, ])) * 1.2
                gridf[1, ] <- seq(min(abs(estimphi[1, ])) * 0.8, max(abs(estimphi[1, ])) * 1.2, length = 500)
                gridf[2, ] <- seq(min(abs(estimphi[2, ])) * 0.8, max(abs(estimphi[2, ])) * 1.2, length = 500)

            
            estimf1 <- rep(0, length(gridf[1, ]), byrow = TRUE)
            estimf2 <- rep(0, length(gridf[2, ]), byrow = TRUE)
            
            for (n in 1:nb.mixt) {
                estimf1 <- estimf1 + mixt.prop[n] * dnorm(gridf[1, ], mean = mu[Niter, n, 1], sd = abs(omega[Niter, 
                  n, 1]))
                estimf2 <- estimf2 + mixt.prop[n] * dnorm(gridf[2, ], mean = mu[Niter, n, 2], sd = abs(omega[Niter, 
                  n, 2]))
            }
            estimf <- rbind(estimf1, estimf2)
        }
        
        if (sum(drift.random) == 1) {
            
                gridf <- matrix(0, 2, 500)
                gridf <- seq(min(abs(estimphi[1, ])) * 0.8, max(abs(estimphi[1, ])) * 1.2, length = 500)

            
            estimf <- rep(0, length(gridf), byrow = TRUE)
            
            for (n in 1:nb.mixt) {
                estimf <- estimf + mixt.prop[n] * dnorm(gridf, mean = mu[Niter, n, 1], sd = abs(omega[Niter, n, 
                  1]))
            }
            
            gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
            estimf <- matrix(estimf, 1, length(estimf), byrow = TRUE)
        }
        
        if (sum(drift.random) == 2) {
            
                gridf <- matrix(0, 2, 500)
                gridf <- seq(min(abs(estimphi[2, ])) * 0.8, max(abs(estimphi[2, ])) * 1.2, length = 500)

            
            estimf <- rep(0, length(gridf), byrow = TRUE)
            
            for (n in 1:nb.mixt) {
                estimf <- estimf + mixt.prop[n] * dnorm(gridf, mean = mu[Niter, n, 2], sd = abs(omega[Niter, n, 
                  2]))
            }
            
            gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
            estimf <- matrix(estimf, 1, length(estimf), byrow = TRUE)
        }
        
        return(new(Class = "Mixture.fit.class", model = model, drift.random = drift.random, gridf = gridf, mu = mu, 
            omega = omega, mixt.prop = mixt.prop, sigma2 = sigma2, index = index, probindi = probindi, 
            estimf = estimf, estimphi = estimphi, bic = bic, aic = aic, estim.drift.fix = as.numeric(is.null(drift.fixed)), times = times, 
            X = X))
        
        
    }
    
}


#' S4 class for the estimation results in the mixed SDE with random effects in the drift, in the diffusion or both 
#'  
#' @slot model character 'OU' or 'CIR'
#' @slot drift.random numeric 0, 1, 2, or c(1,2)
#' @slot diffusion.random numeric 0 or 1
#' @slot gridf matrix of values on which the estimation of the density of the random effects in the drift is done
#' @slot gridg matrix of values on which the estimation of the density of the random effects in the diffusion is done
#' @slot mu numeric estimator of the mean mu of the drift random effects
#' @slot omega numeric estimator of the variance of the drift random effects
#' @slot a numeric estimator of the shape of the Gamma distribution for the diffusion random effect
#' @slot lambda numeric estimator of the scale of the Gamma distribution for the diffusion random effect
#' @slot sigma2 numeric estimated value of \eqn{\sigma^2} if the diffusion coefficient is not random
#' @slot index index of the used trajectories
#' @slot estimphi matrix of the estimator of the drift random effects 
#' @slot estimpsi2 vector of the estimator of the diffusion random effects \eqn{\sigma_j^2}
#' @slot estimf estimator of the (conditional) density of \eqn{\phi}, matrix form
#' @slot estimg estimator of the density of \eqn{\phi}, matrix form
#' @slot estim.drift.fix 1 if the user asked for the estimation of fixed parameter in the drift
#' @slot estim.diffusion.fix 1 if the user asked for the estimation of fixed diffusion coefficient
#' @slot discrete 1 if the estimation is based on the likelihood of discrete observations, 0 otherwise 
#' @slot bic numeric bic 
#' @slot aic numeric aic
#' @slot times vector of observation times, storage of input variable
#' @slot X matrix of observations, storage of input variable


setClass(Class = "Fit.class", representation = representation(model = "character", drift.random = "numeric", diffusion.random = "numeric", 
    gridf = "matrix", gridg = "matrix", mu = "numeric", omega = "numeric", a = "numeric", lambda = "numeric", 
    sigma2 = "numeric", index = "numeric", estimphi = "matrix", estimpsi2 = "matrix", estimf = "matrix", 
    estimg = "matrix", estim.drift.fix = "numeric", estim.diffusion.fix = "numeric", discrete = "numeric", bic = "numeric", 
    aic = "numeric", times = "numeric", X = "matrix"))


#' S4 class for the estimation results when the random effects in the drift follow 
#' mixture of normal distributions  
#'  
#' @slot model character 'OU' or 'CIR'
#' @slot drift.random numeric 1, 2, or c(1,2)
#' @slot gridf matrix of values on which the estimation of the density of the random effects is done
#' @slot mu array estimated value of the mean of the drift random effects at each iteration of the EM algorithm (Niter x nb.mixt x 2)
#' @slot omega array estimated value of the standard deviation of the drift random effects at each iteration of the EM algorithm (Niter x nb.mixt x 2)
#' @slot mixt.prop matrix estimated value of the mixing proportions at each iteration of the EM algorithm (Niter x nb.mixt)
#' @slot sigma2 numeric estimated value of \eqn{\sigma^2}
#' @slot index index of the used trajectories
#' @slot estimphi matrix of the estimator of the drift random effects
#' @slot probindi matrix of posterior component probabilities
#' @slot estimf matrix estimator of the density of the drift random effects
#' @slot estim.drift.fix numeric 1 if the user asked for the estimation of fixed parameter in the drift
#' @slot bic numeric bic 
#' @slot aic numeric aic
#' @slot times vector of observation times, storage of input variable
#' @slot X matrix of observations, storage of input variable


setClass(Class = "Mixture.fit.class", representation = representation(model = "character", drift.random = "numeric", 
    gridf = "matrix", mu = "array", omega = "array", mixt.prop = "matrix", sigma2 = "numeric", 
    index = "numeric", estimphi = "matrix", probindi = "matrix", estimf = "matrix", estim.drift.fix = "numeric", 
    bic = "numeric", aic = "numeric", times = "numeric", X = "matrix"))


########################################################### OUTPUTS
#' Transfers the class object to a list
#' 
#' @description Method for the S4 classes
#' @param x Fit.class or Mixture.fit.class class


out <- function(x) {
    sN <- slotNames(x)
    res <- lapply(sN, function(name) slot(x, name))
    names(res) <- sN
    res
}



# Summary for Fit.class -----------------------------------------------------------------



############################################################### SUMMARY

#' Short summary of the results of class object Fit.class
#' @description Method for the S4 class Fit.class
#' @param object Fit.class class
#'
#'
setMethod(f = "summary", signature = "Fit.class", definition = function(object) {
    
    cat(c("Number of trajectories used for estimation: ",length(object@index),"\n"))
  
    if (sum(object@drift.random) > 2) {
        
        cat("\nRandom effects in the drift:\n")
        
        driftreff1 <- (matrix(c(object@mu[1], object@omega[1]), 2, 1, byrow = TRUE))
        
        rownames(driftreff1) <- c("MLE mean 1", "MLE sd 1")
        
        print(driftreff1, quote = FALSE, right = TRUE)

        driftreff2 <- (matrix(c(object@mu[2], object@omega[2]), 2, 1, byrow = TRUE))
        
        rownames(driftreff2) <- c("MLE mean 2", "MLE sd 2")
        
        print(driftreff2, quote = FALSE, right = TRUE)
        
    }
    
    if (sum(object@drift.random) == 2) {
        cat("\nFixed and random effects in the drift:\n")
        
        driftreff <- (matrix(c(object@mu[1], object@mu[2], object@omega[2]), 3, 1, byrow = TRUE))
        
        if (object@estim.drift.fix == 1) {
            rownames(driftreff) <- c("MLE fixed effect", "MLE mean", "MLE sd")
        }
        if (object@estim.drift.fix == 0) {
            rownames(driftreff) <- c("Fixed effect (not estimated)", "MLE mean", "MLE sd")
        }
        
        print(driftreff, quote = FALSE, right = TRUE)
       
    }
    
    if (sum(object@drift.random) == 1) {
        cat("\nFixed and random effects in the drift:\n")
        
        driftreff <- (matrix(c(object@mu[1], object@omega[1], object@mu[2]), 3, 1, byrow = TRUE))
        
        if (object@estim.drift.fix == 1) {
            rownames(driftreff) <- c("MLE mean", "MLE sd", "MLE fixed effect")
        }
        if (object@estim.drift.fix == 0) {
            rownames(driftreff) <- c("MLE mean", "MLE sd", "Fixed effect (not estimated)")
        }
        
        print(driftreff, quote = FALSE, right = TRUE)
        
    }
    
    if (sum(object@drift.random) == 0) {
        
        if (object@estim.drift.fix == 1) {
            cat("\nFixed effects in the drift:\n")
            
            fixedeff <- (matrix(c(object@mu[1], object@mu[2]), 2, 1, byrow = TRUE))
            
            rownames(fixedeff) <- c("MLE fixed effect 1", "MLE fixed effect 2")
            
            print(fixedeff, quote = FALSE, right = TRUE)
        }
        
        if (object@estim.drift.fix == 0) {
            cat("\nFixed effects in the drift:\n")
            
            fixedeff <- (matrix(c(object@mu[1], object@mu[2]), 2, 1, byrow = TRUE))
            
            rownames(fixedeff) <- c("Fixed effect 1 (not estimated)", "Fixed effect 2 (not estimated)")
            
            print(fixedeff, quote = FALSE, right = TRUE)
        }
        
        
    }
    
    
    if (object@diffusion.random == 1) {
        
        cat("\nRandom effects in the diffusion coefficient:\n")
        
        diffreff <- (matrix(c(object@a, object@lambda, object@a * object@lambda, digamma(object@a) + log(object@lambda)), 
            4, 1, byrow = TRUE))
        
        rownames(diffreff) <- c("a (shape)", "lambda (scale)", "m = a*lambda", "t = psi(a) + log(lambda)")
        
        print(diffreff, quote = FALSE, right = TRUE)
        
    }
    
    if (object@diffusion.random == 0) {
        
        cat("\nFixed effect in the diffusion coefficient:\n")
        
        if ((object@estim.diffusion.fix == 1) && (object@discrete == 1)) {
            print(matrix(c("MLE sigma", round(sqrt(object@sigma2), 6)), 1, 2, byrow = TRUE), quote = FALSE, right = TRUE)
        }
        
        if ((object@estim.diffusion.fix == 1) && (object@discrete == 0)) {
            print(matrix(c("Estimated sigma (quadratic variations)", round(sqrt(object@sigma2), 6)), 1, 2, byrow = TRUE), 
                quote = FALSE, right = TRUE)
        }
        
        if (object@estim.diffusion.fix == 0) {
            print(matrix(c("sigma (not estimated)", round(sqrt(object@sigma2), 6)), 1, 2, byrow = TRUE), quote = FALSE, 
                right = TRUE)
        }
        
        
    }
    
    
    cat("\nCriteria for model selection:\n")
    
    
    info.criteria <- matrix(c(round(object@bic, 6), round(object@aic, 6)), 2, 1, byrow = TRUE)
    rownames(info.criteria) <- c("BIC", "AIC")
    colnames(info.criteria) <- c("")
    print(info.criteria, quote = FALSE, right = TRUE)
    
    
})


############################################################### SUMMARY

#' Short summary of the results of class object Mixture.fit.class
#' @description Method for the S4 class Mixture.fit.class
#' @param object Mixture.fit.class class
#'
setMethod("summary", "Mixture.fit.class", function(object) {
    
    cat(c("Number of trajectories used for estimation: ",length(object@index),"\n"))
  
    cat(c("Number of iterations of the EM algorithm", dim(object@mixt.prop)[1], "\n"))
  
    Niter <- dim(object@mixt.prop)[1]
    nb.mixt <- dim(object@mixt.prop)[2]
    
    if (object@bic != 0) {
        
        cat("\nMLE mixing proportions:\n")
        
        prop <- matrix(object@mixt.prop[Niter, ], 1, nb.mixt, byrow = TRUE)
        col.name <- c()
        for (j in 1:nb.mixt) {
            col.name <- c(col.name, paste("Comp.", j))
        }
        colnames(prop) <- col.name
        rownames(prop) <- c("MLE prop.")
        
        print(prop, quote = FALSE, right = TRUE)
        
        if (sum(object@drift.random) > 2) {
            
            cat("\nRandom effects in the drift:\n")
            
            driftreff <- (matrix(c(object@mu[Niter, , 1], object@omega[Niter, , 1], object@mu[Niter, , 2], object@omega[Niter, 
                , 2]), 4, nb.mixt, byrow = TRUE))
            
            rownames(driftreff) <- c("MLE mean 1", "MLE sd 1", "MLE mean 2", "MLE sd 2")
            colnames(driftreff) <- col.name
            
            print(driftreff, quote = FALSE, right = TRUE)
            
            
            
        }
        
        
        if (sum(object@drift.random) == 2) {
            cat("\nFixed and random effects in the drift:\n")
            
            driftreff <- (matrix(c(object@mu[Niter, , 1], object@mu[Niter, , 2], object@omega[Niter, , 2]), 3, 
                nb.mixt, byrow = TRUE))
            
            if (object@estim.drift.fix == 1) {
                rownames(driftreff) <- c("MLE fixed effect", "MLE mean", "MLE sd")
            }
            if (object@estim.drift.fix == 0) {
                rownames(driftreff) <- c("Fixed effect (not estimated)", "MLE mean", "MLE sd")
            }
            
            colnames(driftreff) <- col.name
            
            print(driftreff, quote = FALSE, right = TRUE)
            
        }
        
        if (sum(object@drift.random) == 1) {
            cat("\nFixed and random effects in the drift:\n")
            
            driftreff <- (matrix(c(object@mu[Niter, , 1], object@omega[Niter, , 1], object@mu[Niter, , 2]), 3, 
                nb.mixt, byrow = TRUE))
            
            if (object@estim.drift.fix == 1) {
                rownames(driftreff) <- c("MLE mean", "MLE sd", "MLE fixed effect")
            }
            if (object@estim.drift.fix == 0) {
                rownames(driftreff) <- c("MLE mean", "MLE sd", "Fixed effect (not estimated)")
            }
            
            colnames(driftreff) <- col.name
            print(driftreff, quote = FALSE, right = TRUE)
            
        }
        
        
        cat("\nCriteria for model selection:\n")
        
        
        info.criteria <- matrix(c(round(object@bic, 6), round(object@aic, 6)), 2, 1, byrow = TRUE)
        rownames(info.criteria) <- c("BIC", "AIC")
        colnames(info.criteria) <- c("")
        print(info.criteria, quote = FALSE, right = TRUE)
        
    }
    
})


########################################################### PLOT
#' Plot method for the estimation class object
#'
#' @description Plot method for the S4 class Fit.class
#' @param x Fit.class class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#'
#' @importFrom grDevices x11
#' @importFrom graphics par
#' @importFrom graphics plot

setMethod(f = "plot", signature = "Fit.class", definition = function(x, newwindow = FALSE, ...) {
    if (newwindow) {
        x11(width = 14)
    }
    
    
    if (x@diffusion.random == 0) {
        if (dim(x@gridf)[1] == 1) {
            
            op <- par(mfrow = c(1, 1), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(1, 1, 1, 1), omi = c(0.2, 
                0.2, 0.2, 0.2), cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
            if (x@drift.random==1){
              plot(x@gridf, x@estimf, main = expression(alpha[j]), xlab = "Value", ylab = "Density", type = "l")
            }
            if (x@drift.random==2){
              plot(x@gridf, x@estimf, main = expression(beta[j]), xlab = "Value", ylab = "Density", type = "l")
            }    
        }
        
        if (dim(x@gridf)[1] == 2) {
            
            op <- par(mfrow = c(2, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(1, 1, 1, 1), omi = c(0.2, 
                0.2, 0.2, 0.2), cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
            plot(x@gridf[1, ], x@estimf[1, ], main = expression(alpha[j]), xlab = "Value", ylab = "Density", type = "l")
            plot(x@gridf[2, ], x@estimf[2, ], main = expression(beta[j]), xlab = "Value", ylab = "Density", type = "l")
            title("Estimated densities", outer = TRUE)
            
        }
        
    } else {
        
        # if (x@gridf==0) { op <- par(mfrow = c(1, 1), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(1, 1, 1,
        # 1), omi = c(0.2, 0.2, 0.2, 0.2), cex.main = 1, cex.lab = 0.9, cex.axis = 0.9) plot(x@gridg, x@estimg, main =
        # 'Estimated density of the random effect in the diffusion',xlab='Value',ylab='Density') }
        
        if (dim(x@gridf)[1] == 1) {
          
          if (dim(x@gridf)[2] == 1) {
            
            op <- par(mfrow = c(1, 1), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(1, 1, 1, 1), omi = c(0.2, 
                                                                                                                   0.2, 0.2, 0.2), cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
            plot(x@gridg, x@estimg, main = expression(sigma[j]), xlab = "Value", ylab = "Density", type = "l")
            title("Estimated density", outer = TRUE)
            #"Random effect in the diffusion"
          } else {
            
            op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(1, 1, 1, 1), omi = c(0.2, 
                0.2, 0.2, 0.2), cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
            if (x@drift.random==1){
              plot(x@gridf, x@estimf, main = expression(alpha[j]), xlab = "Value", ylab = "Density", type = "l")
            }
            if (x@drift.random==2){
              plot(x@gridf, x@estimf, main = expression(beta[j]), xlab = "Value", ylab = "Density", type = "l")
            }
            plot(x@gridg, x@estimg, main = expression(sigma[j]), xlab = "Value", ylab = "Density", type = "l")
            title("Estimated densities", outer = TRUE)
          }
        }
        
        if (dim(x@gridf)[1] == 2) {
            
            op <- par(mfrow = c(2, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(1, 1, 1, 1), omi = c(0.2, 
                0.2, 0.2, 0.2), cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
            plot(x@gridf[1, ], x@estimf[1, ], main = expression(alpha[j]), xlab = "Value", ylab = "Density", type = "l")
            plot(x@gridf[2, ], x@estimf[2, ], main = expression(beta[j]), xlab = "Value", ylab = "Density", type = "l")
            plot(x@gridg, x@estimg, main = expression(sigma[j]), xlab = "Value", ylab = "Density", type = "l")
            title("Estimated densities", outer = TRUE)
            
        }
        
    }
    
})


########################################################### PLOT
#' Plot method for the mixture estimation class object
#'
#' @description Plot method for the S4 class Mixture.fit.class
#' @param x Mixture.fit.class class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#'
#' @importFrom grDevices x11
#' @importFrom graphics par
#' @importFrom graphics title
#' @importFrom graphics layout

setMethod(f = "plot", signature = "Mixture.fit.class", definition = function(x, newwindow = FALSE, ...) {
    if (newwindow) {
        x11(width = 14)
    }
    
    # Convergence plot of the EM algorithm
    
    
    mat = matrix(c(1, 3, 5, 2, 4, 5), 3, 2)
    op <- par(mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), omi = c(1, 1, 1, 1), oma = c(1, 1, 1, 1), cex.main = 1, 
        cex.lab = 0.9, cex.axis = 0.9)
    # mfcol = c(3, 2),
    
    nb.mixt <- dim(x@mu)[2]
    for (k in 1:nb.mixt) {
        # plot.new()
        layout(mat)
        tit <- paste("Convergence plot for mixture component ", k)
        plot(x@mu[, k, 1], type = "l", xlab = "Iteration", ylab = "")
        title("mu1")
        plot(x@mu[, k, 2], type = "l", xlab = "Iteration", ylab = "")
        title("mu2")
        plot(x@omega[, k, 1], type = "l", xlab = "Iteration", ylab = "")
        title("omega1")
        plot(x@omega[, k, 2], type = "l", xlab = "Iteration", ylab = "")
        title("omega2")
        plot(x@mixt.prop[, k], type = "l", xlab = "Iteration", ylab = "")
        title("Proportion")
        title(tit, outer = TRUE)
    }
    
    # Plots of the estimated densities
    
    if (sum(x@drift.random) == 3) {
        op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), omi = c(1, 1, 1, 1), oma = c(1, 
            1, 1, 1), cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
        
        
        plot(x@gridf[1, ], x@estimf[1, ], ylab = "Density", xlab = "value", type = "l")
        title(expression(alpha[j]))
        plot(x@gridf[2, ], x@estimf[2, ], ylab = "Density", xlab = "value", type = "l")
        title(expression(beta[j]))
        title("Estimated densities", outer = TRUE)
    }
    
    if (sum(x@drift.random) <= 2) {
        op <- par(mfrow = c(1, 1), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), omi = c(1, 1, 1, 1), oma = c(1, 
            1, 1, 1), cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
        
        
        plot(x@gridf, x@estimf, ylab = "Density", xlab = "value", type = "l")
        
        if (x@drift.random==1){
          title(expression(alpha[j]), outer = TRUE)
        }
        if (x@drift.random==2){
          title(expression(beta[j]), outer = TRUE)
        }
        
    }
    
})

#' ########################################################### VALIDATION

#' Validation of the chosen model.
#'
#' @description Validation of the chosen model. For the index numj, Mrep=100 new trajectories are simulated
#' with the value of the estimated random effect number numj. Two plots are given: on the left the simulated trajectories and the true one (red)
#' and one the left the corresponding qq-plot for each time.
#' @param x Fit.class or Mixture.fit.class class
#' @param ... other optional parameters
#'
setGeneric("valid", function(x, ...) {
    standardGeneric("valid")
})


#' Validation of the chosen model.
#'
#' @description Validation of the chosen model. For the index numj, Mrep=100 new trajectories are simulated
#' with the value of the estimated random effect number numj. Two plots are given: on the left the simulated trajectories and the true one (red)
#' and one the left the corresponding qq-plot for each time.
#' @param x Fit.class class
#' @param Mrep number of trajectories to be drawn
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param plot.valid logical(1), if TRUE, the results are depicted grafically
#' @param numj optional number of series to be validated
#' @param ... optional plot parameters
#' @importFrom grDevices x11
#' @importFrom graphics par
#' @importFrom graphics abline



setMethod(f = "valid", signature = "Fit.class", definition = function(x, Mrep = 100, newwindow = FALSE, plot.valid = TRUE, 
    numj, ...) {
    
    if (newwindow) {
        x11(width = 10)
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
    
    Xtrue <- x@X
    times <- round(x@times, 10)
    Tend <- max(times)
    del <- round(min(diff(times)), 10)
    timessimu <- round(seq(del, Tend, by = del), 10)
    
    M <- dim(Xtrue)[1]
    
    
    if (missing(numj)) {
        numj <- floor(runif(1, 1, M))
    }
    if (!missing(numj)) {
        numj <- numj
    }
    
    
    if (dim(x@gridf)[1] == 2) {
        
        phihat <- x@estimphi[, numj]
        
        if (x@diffusion.random == 0) {
            sig <- sqrt(x@sigma2)
            
            
            if (x@model == "OU") {
                suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                  method = "EA", theta = c(phihat, sig), model = "OU", M = Mrep)))
            }
            if (x@model == "CIR") {
                suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                  method = "milstein", theta = c(phihat, sig), model = "CIR", M = Mrep, sigma.x = expression(sig/(2 * 
                    sqrt(x))), sigma = expression(sig * sqrt(x)))))
            }
        }
        
        if (x@diffusion.random == 1) {
            phihat <- x@estimphi[, numj]
            psihat <- sqrt(x@estimpsi2[numj])
            
            if (x@model == "OU") {
                suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                  method = "EA", theta = c(phihat, psihat), model = "OU", M = Mrep)))
            }
            if (x@model == "CIR") {
                suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                  method = "milstein", theta = c(phihat, psihat), model = "CIR", M = Mrep, sigma.x = expression(psihat/(2 * 
                    sqrt(x))), sigma = expression(psihat * sqrt(x)))))
            }
        }
    }
    
    if (dim(x@gridf)[1] == 1) {
        phihat <- x@estimphi[numj]
        
        if (x@diffusion.random == 0) {
            sig <- sqrt(x@sigma2)
            if (sum(x@drift.random) == 1) {
                paramfixed <- x@mu[2]
                if (x@model == "OU") {
                  suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                    method = "EA", theta = c(phihat, paramfixed, sig), model = "OU", M = Mrep)))
                }
                
                if (x@model == "CIR") {
                  suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                    method = "milstein", theta = c(phihat, paramfixed, sig), model = "CIR", M = Mrep, sigma.x = expression(sig/(2 * 
                      sqrt(x))), sigma = expression(sig * sqrt(x)))))
                }
            }
            if (sum(x@drift.random) == 2) {
                paramfixed <- x@mu[1]
                if (x@model == "OU") {
                  suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                    method = "EA", theta = c(paramfixed, phihat, sig), model = "OU", M = Mrep)))
                }
                if (x@model == "CIR") {
                  suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                    method = "milstein", theta = c(paramfixed, phihat, sig), model = "CIR", M = Mrep, sigma.x = expression(sig/(2 * 
                      sqrt(x))), sigma = expression(sig * sqrt(x)))))
                }
            }
        }
        
        if (x@diffusion.random == 1) {
            psihat <- sqrt(x@estimpsi2[numj])
            if (sum(x@drift.random) == 1) {
                paramfixed <- x@mu[2]
                if (x@model == "OU") {
                  suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                    method = "EA", theta = c(phihat, paramfixed, psihat), model = "OU", M = Mrep)))
                }
                
                if (x@model == "CIR") {
                  suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                    method = "milstein", theta = c(phihat, paramfixed, psihat), model = "CIR", M = Mrep, sigma.x = expression(psihat/(2 * 
                      sqrt(x))), sigma = expression(psihat * sqrt(x)))))
                }
            }
            if (sum(x@drift.random) == 2) {
                paramfixed <- x@mu[1]
                if (x@model == "OU") {
                  suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                    method = "EA", theta = c(paramfixed, phihat, psihat), model = "OU", M = Mrep)))
                }
                if (x@model == "CIR") {
                  suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                    method = "milstein", theta = c(paramfixed, phihat, psihat), model = "CIR", M = Mrep, sigma.x = expression(psihat/(2 * 
                      sqrt(x))), sigma = expression(psihat * sqrt(x)))))
                }
            }
        }
    }
    
    vecttimes <- intersect(round(timessimu, 10), round(times, 10))
    
    N <- length(vecttimes)
    
    q <- rep(0, N)
    for (i in 1:N) {
        q[i] <- sum(Xtrue[numj, which(times == vecttimes[i])[1]] > Xnew[, which(timessimu == vecttimes[i])[1]])/Mrep
    }
    
    if (plot.valid == 1) {
        numj <- numj
        op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, 
            cex.lab = 0.7, cex.axis = 0.7)
        
        plot(c(0, timessimu), Xnew[1, ], type = "l", ylim = c(min(Xnew) * 0.8, max(Xnew) * 1.2), xlab = "", ylab = "")
        for (k in 1:Mrep) {
            lines(c(0, timessimu), Xnew[k, ])
        }
        lines(times, Xtrue[numj, ], col = "red", lwd = 2)
        
        
        plot(1:N/N, sort(q), xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
        abline(0, 1)
    }
    
    return(list(quantiles = q, Xnew = Xnew, numj = numj))
})


#' Validation of the chosen model for the mixture class.
#'
#' @description Validation of the chosen model. For the index numj, Mrep=100 new trajectories are simulated
#' with the value of the estimated random effect number numj. Two plots are given: on the left the simulated trajectories and the true one (red)
#' and one the left the corresponding qq-plot for each time.
#' @param x Mixture.fit.class class
#' @param Mrep number of trajectories to be drawn
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param plot.valid logical(1), if TRUE, the results are depicted grafically
#' @param numj optional number of series to be validated
#' @param ... optional plot parameters
#' @importFrom grDevices x11
#' @importFrom graphics par
#' @importFrom graphics abline


setMethod(f = "valid", signature = "Mixture.fit.class", definition = function(x, Mrep = 100, newwindow = FALSE, 
    plot.valid = TRUE, numj, ...) {
    
    if (newwindow) {
        x11(width = 10)
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
    
    Niter <- dim(x@mu)[1]
    
    Xtrue <- x@X
    times <- round(x@times, 10)
    Tend <- max(times)
    del <- round(min(diff(times)), 10)
    timessimu <- round(seq(del, Tend, by = del), 10)
    
    M <- dim(Xtrue)[1]
    
    
    if (missing(numj)) {
        numj <- floor(runif(1, 1, M))
    }
    if (!missing(numj)) {
        numj <- numj
    }
    
    # if (x@diffusion.random == 1) {
    #     warning("For the considered mixtures of SDE, there should not be random effects in the diffusion coeffcicient")
    #     diffusion.random <- 0
    # }
    # if (sum(x@drift.random) == 0) {
    #     stop("There should be at least one random effect in the drift")
    # }
    
    sig <- sqrt(x@sigma2)
    
    if (dim(x@gridf)[1] == 2) {
        
        phihat <- x@estimphi[, numj]
        
        if (x@model == "OU") {
            suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                method = "EA", theta = c(phihat, sig), model = "OU", M = Mrep)))
        }
        if (x@model == "CIR") {
            suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                method = "milstein", theta = c(phihat, sig), model = "CIR", M = Mrep, sigma.x = expression(sig/(2 * 
                  sqrt(x))), sigma = expression(sig * sqrt(x)))))
        }
    }
    
    if (dim(x@gridf)[1] == 1) {
        phihat <- x@estimphi[numj]
        component <- which.max(x@probindi[numj,])
        
        sig <- sqrt(x@sigma2)
        if (sum(x@drift.random) == 1) {
            paramfixed <- x@mu[Niter,component,2]
            if (x@model == "OU") {
              
                suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                  method = "EA", theta = c(phihat[1], paramfixed, sig), model = "OU", M = Mrep)))
            }
            
            if (x@model == "CIR") {
                suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                  method = "milstein", theta = c(phihat[1], paramfixed, sig), model = "CIR", M = Mrep, sigma.x = expression(sig/(2 * 
                    sqrt(x))), sigma = expression(sig * sqrt(x)))))
            }
        }
        if (sum(x@drift.random) == 2) {
            paramfixed <- x@mu[Niter,component,1]
            if (x@model == "OU") {
                suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                  method = "EA", theta = c(paramfixed, phihat[2], sig), model = "OU", M = Mrep)))
            }
            if (x@model == "CIR") {
                suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), delta = del, 
                  method = "milstein", theta = c(paramfixed, phihat[2], sig), model = "CIR", M = Mrep, sigma.x = expression(sig/(2 * 
                    sqrt(x))), sigma = expression(sig * sqrt(x)))))
            }
        }
    }
    
    vecttimes <- intersect(round(timessimu, 10), round(times, 10))
    
    N <- length(vecttimes)
    
    q <- rep(0, N)
    for (i in 1:N) {
        q[i] <- sum(Xtrue[numj, which(times == vecttimes[i])[1]] > Xnew[, which(timessimu == vecttimes[i])[1]])/Mrep
    }
    
    if (plot.valid == 1) {
        numj <- numj
        op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), cex.main = 0.8, 
            cex.lab = 0.7, cex.axis = 0.7)
        
        plot(c(0, timessimu), Xnew[1, ], type = "l", ylim = c(min(Xnew) * 0.8, max(Xnew) * 1.2), xlab = "", ylab = "")
        for (k in 1:Mrep) {
            lines(c(0, timessimu), Xnew[k, ])
        }
        lines(times, Xtrue[numj, ], col = "red", lwd = 2)
        
        
        plot(1:N/N, sort(q), xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
        abline(0, 1)
    }
    
    return(list(quantiles = q, Xnew = Xnew, numj = numj))
})


 
