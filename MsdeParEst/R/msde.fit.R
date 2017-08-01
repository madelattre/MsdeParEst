#' Estimation Of The Random Effects In Mixed Stochastic Differential Equations
#' 
#' 
#' @description Parametric estimation of the random effects \eqn{(\alpha_j, \beta_j, \sigma_j)} joint density in the mixed SDE
#'  \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j a(X_j(t)) dW_j(t)}.
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
#' @param Niter default 10, number of iterations for the EM algorithm if mixture='paramMLmixture'
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
#' \item{estim.drift.fixed}{initial choice}
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
#' The drift includes one or two random effects and the diffusion includes either a fixed effect or a random effect. Two diffusions are implemented.
#' \subsection{Ornstein-Uhlenbeck model (OU), if diffusion.random = 0}{
#' 
#' If drift.random = 1, \eqn{\beta} is a fixed effect: \eqn{dX_j(t)= (\alpha_j- \beta X_j(t))dt + \sigma dW_j(t)  } 
#' 
#' If drift.random = 2, \eqn{\alpha} is a fixed effect: \eqn{dX_j(t)= (\alpha - \beta_j X_j(t))dt + \sigma dW_j(t)  }
#' 
#' If drift.random = c(1,2), \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma dW_j(t)  } 
#' }
#' \subsection{Ornstein-Uhlenbeck model (OU), if diffusion.random = 1}{
#' 
#' If drift.random = 1, \eqn{\beta} is a fixed effect: \eqn{dX_j(t)= (\alpha_j- \beta X_j(t))dt + \sigma_j dW_j(t)  } 
#' 
#' If drift.random = 2, \eqn{\alpha} is a fixed effect: \eqn{dX_j(t)= (\alpha - \beta_j X_j(t))dt + \sigma_j dW_j(t)  }
#' 
#' If drift.random = c(1,2), \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j dW_j(t)  } 
#' }
#' \subsection{Cox-Ingersoll-Ross model (CIR), if diffusion.random = 0}{
#' 
#' If drift.random = 1, \eqn{\beta} is a fixed effect: \eqn{dX_j(t)= (\alpha_j- \beta X_j(t))dt + \sigma \sqrt{X_(t)} dWj_(t)  } 
#' 
#' If drift.random = 2, \eqn{\alpha} is a fixed effect: \eqn{dX_j(t)= (\alpha - \beta_j X_j(t))dt + \sigma \sqrt{X_j(t)} dW_j(t)  } 
#' 
#' If drift.random = c(1,2), \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma \sqrt{X_j(t)}  dW_j(t)  } 
#'}
#' \subsection{Cox-Ingersoll-Ross model (CIR), if diffusion.random = 1}{
#' 
#' If drift.random = 1, \eqn{\beta} is a fixed effect: \eqn{dX_j(t)= (\alpha_j- \beta X_j(t))dt + \sigma_j \sqrt{X_(t)} dWj_(t)  } 
#' 
#' If drift.random = 2, \eqn{\alpha} is a fixed effect: \eqn{dX_j(t)= (\alpha - \beta_j X_j(t))dt + \sigma_j \sqrt{X_j(t)} dW_j(t)  } 
#' 
#' If drift.random = c(1,2), \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j \sqrt{X_j(t)}  dW_j(t)  } 
#'}
#'
#' If diffusion.random = 0, the random effects in the drift follow a Gaussian distribution (or a mixture of Gaussian distributions) with mean mu and standard deviation omega.
#' If diffusion.random = 1, the random effects in the diffusion \eqn{sigma_j^2} follow an Inverse Gamma distribution with shape a and scale lambda, and the random effects in the drift follow a Gaussian distribution conditional to the random effect in the drift, with mean mu and standard deviation omega*sigma^2.  
#'
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
#' model <- "OU"
#' drift.random <- 2
#' diffusion.random <- 1
#' drift.fixed <- 0
#' density.phi <- 'normal'
#' drift.param <- c(0.5,0.5)
#' diffusion.param <- c(8,1/2)
#'
#' sim1 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random, 
#'                  diffusion.random = diffusion.random, drift.fixed = drift.fixed, 
#'                  density.phi = density.phi, drift.param = drift.param, 
#'                  diffusion.param = diffusion.param)
#'                  
#' # -- Estimation
#' 
#' # -----Fixed effect in the drift estimated
#' res1 <- msde.fit(times = sim1$times, X = sim1$X, model = "OU", drift.random = 2, 
#'                  diffusion.random = 1, estim.drift.fix = 1, mixture = 0)
#' summary(res1)
#' print(res1)
#' #pred1 <- pred(res1, invariant = 0, level = 0.05, newwindow = FALSE, plot.pred = TRUE)
#' valid(res1)
#' plot(res1)
#' 
#' # ----- Fixed effect in the drift known and not estimated
#' res1bis <- msde.fit(times = sim1$times, X = sim1$X, model = "OU", drift.random = 2, 
#'                     diffusion.random = 1, drift.fixed=0, mixture = 0)
#' summary(res1bis)
#' 
#' # Example 2: one random effect in the drift and one fixed effect in the diffusion coefficient
#' 
#' # -- Simulation
#' M <- 100
#' Tmax <- 5
#' N <- 5000
#' model <- "OU"
#' diffusion.random <- 0
#' diffusion.param <- 0.5
#' drift.random <- 2
#' drift.fixed <- 10
#' density.phi <- 'normal'
#' drift.param <- c(1,sqrt(0.4/4))
#' 
#' sim2 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random, 
#'                  diffusion.random = diffusion.random, drift.fixed = drift.fixed, 
#'                  density.phi = density.phi, drift.param = drift.param, 
#'                  diffusion.param = diffusion.param)
#' 
#' # -- Estimation
#' res2 <- msde.fit(times = sim2$times, X = sim2$X, model = "OU", drift.random = 2, 
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
#' model <- "OU"
#' drift.random <- c(1,2)
#' diffusion.random <- 1
#' density.phi <- 'normalnormal'
#' drift.param <- c(1,0.5,0.5,0.5)
#' diffusion.param <- c(8,1/2)
#' 
#' sim3 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random, 
#'                  diffusion.random = diffusion.random, density.phi = density.phi, 
#'                  drift.param = drift.param, diffusion.param = diffusion.param)
#' 
#' # -- Estimation
#' 
#' res3 <- msde.fit(times = sim3$times, X = sim3$X, model = "OU", drift.random = c(1,2), 
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
#' 
#' model <- "OU"
#' drift.random <- 0
#' diffusion.random <- 1
#' drift.fixed <- c(0,1)
#' diffusion.param <- c(5,3)
#' sim4 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random, 
#'                  diffusion.random = diffusion.random, drift.fixed = drift.fixed, 
#'                  diffusion.param = diffusion.param)
#'
#' # -- Estimation
#' res4 <- msde.fit(times = sim4$times, X = sim4$X, model = "OU", drift.random = 0, 
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
#' model <- "OU"
#' drift.random <- 1
#' drift.fixed <- 1
#' density.phi <- 'mixture.normal'
#' nb.mixt <- 2
#' mixt.prop <- c(0.5,0.5)
#' param.ea1 <- c(0.5, 0.25, 1.8, 0.25)
#' param.ea2 <- c(1, 0.25, 1, 0.25) 
#' drift.param <- param.ea1
#' 
#' sim5 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random,
#'                  diffusion.random = diffusion.random, drift.fixed = drift.fixed,
#'                  density.phi = density.phi, drift.param = drift.param, 
#'                  diffusion.param = diffusion.param, nb.mixt = nb.mixt, mixt.prop = mixt.prop)
#'
#' # -- Estimation
#' res5 <- msde.fit(times = sim5$times, X = sim5$X, model = "OU", drift.random = 1, 
#'                  estim.drift.fix = 1, diffusion.random = 0, estim.diffusion.fix = 1, 
#'                  mixture = 1, nb.mixt=2, Niter = 25)
#' 
#' summary(res5)
#' print(res5)
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


msde.fit <- function(times, X, model = c("OU", "CIR"), drift.random, drift.fixed = NULL, 
                     estim.drift.fix = 0, diffusion.random = 0, diffusion.fixed = NULL, estim.diffusion.fix = 0, 
                     mixture = 0,  
                     nb.mixt = 1, drift.fixed.mixed = 0, Niter = 10, discrete = 1) {
  model <- match.arg(model)
  #mixture <- match.arg(mixture)
  
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
    
    gridf <- NULL
    
    estimg <- 0
    gridg <- 0
    estimf <- 0
    #gridf <- 0
    
    if (model == "OU") {
      Mindex <- M
      index <- 1:Mindex
    }
    
    if (model == "CIR") {
      
      index <- which(rowSums(X <= 0) == 0)
      Mindex <- length(index)
      if (Mindex == 0) {
        warning("All the trajectories have non positive values the model CIR cannot be used", 
                call. = FALSE)
        estimf <- 0
        estimphi <- 0
        estimpsi2 <- 0
        bic <- 0
        aic <- 0
        gridf <- 0
        mu <- 0
        omega <- 0
        cutoff <- 0
        sigma2 <- 0
        a <- 0
        lambda <- 0
        estimf.trunc <- 0
        estimphi.trunc <- 0
        estimpsi2.trunc <- 0
      }
      
    }
    
    
    if (Mindex > 0) {
      
      # -- Errors and warnings
      
      if (((sum(drift.random)) == 0) && (diffusion.random == 0)) {
        stop("There should be at least one random effect either in the drift or the diffusion coefficient.")
      }
      
      if (diffusion.random == 0) {
        if ((is.null(diffusion.fixed)) && (estim.diffusion.fix == 0)) {
          warning("estim.diffusion.fix == 0. No value is specified for the diffusion coefficient. It is therefore estimated.")
          estim.diffusion.fix <- 1
        }
        
        if ((!is.null(diffusion.fixed)) && (estim.diffusion.fix == 1)) {
          warning("estim.diffusion.fix == 1. The diffusion coefficient is estimated. The specified value is therefore not used.")
        }
      }
      
      
      if ((is.null(drift.fixed)) && (estim.drift.fix == 0)) {
        warning("estim.drift.fix == 0. No value is specified for the drift parameters. They are therefore estimated.")
        estim.drift.fix <- 1
      }
      
      if ((!is.null(drift.fixed)) && (estim.drift.fix == 1)) {
        warning("estim.drift.fix == 1. The drift parameters are estimated. The specified values are therefore not used.")
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
        
        res <- EstParamNormal(U, V, S, SigDelta, K, drift.fixed, estim.drift.fix, 
                              diffusion.fixed, drift.random, estim.diffusion.fix, discrete)
        
        bic <- res$BIChere
        aic <- res$AIChere
        mu <- res$mu
        omega <- res$omega
        sigma2 <- res$sigma^2
        
        
        gridg <- as.matrix(gridg)
        estimg <- as.matrix(estimg)
        
        # -- Estimator of the density
        
        estimpsi2 <- matrix(estimpsi2, 1, length(estimpsi2), byrow = TRUE)
        
        if (sum(drift.random) == 3) {
          
          if (is.null(gridf) == 1) {
            gridf <- matrix(0, 2, 500)
            gridf[1, ] <- seq(mu[1]-3*omega[1], mu[1]+3*omega[1], 
                              length = 500)
            gridf[2, ] <- seq(mu[2]-3*omega[2], mu[2]+3*omega[2], 
                              length = 500)
          }
          
          if (is.null(gridf) == 0) {
            gridf <- gridf
          }
          
          estimf1 <- dnorm(gridf[1, ], mean = mu[1], sd = abs(omega[1]))
          estimf2 <- dnorm(gridf[2, ], mean = mu[2], sd = abs(omega[2]))
          
          estimf <- rbind(estimf1,estimf2)
          
          estimf.trunc <- estimf
          estimphi.trunc <- estimphi
          
          cutoff <- FALSE
        }
        
        if (sum(drift.random) == 2) {
          
          estimphi <- estimphi[2, ]
          
          if (is.null(gridf) == 1) {
            gridf <- seq(mu[2]-3*omega[2], mu[2]+3*omega[2], 
                         length = 500)
          }
          
          if (is.null(gridf) == 0) {
            gridf <- gridf
          }
          
          estimf <- matrix(dnorm(gridf, mean = mu[2], sd = omega[2]), 1, length(gridf), 
                           byrow = TRUE)
          gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
          estimphi <- matrix(estimphi, 1, length(estimphi), byrow = TRUE)
          estimphi.trunc <- estimphi
          estimf.trunc <- estimf
          
          cutoff <- FALSE
          
        }
        
        if (sum(drift.random) == 1) {
          
          estimphi <- estimphi[1, ]
          
          if (is.null(gridf) == 1) {
            gridf <- seq(mu[1]-3*omega[1], mu[1]+3*omega[1], 
                         length = 500)  }
          
          if (is.null(gridf) == 0) {
            gridf <- gridf
          }
          
          estimf <- matrix(dnorm(gridf, mean = mu[1], sd = omega[1]), 1, length(gridf), 
                           byrow = TRUE)
          gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
          estimphi <- matrix(estimphi, 1, length(estimphi), byrow = TRUE)
          estimphi.trunc <- estimphi
          estimf.trunc <- estimf
          
          cutoff <- FALSE
          
        }
        
      }
      
      if ((sum(drift.random) == 0) && (diffusion.random == 1)) {
        
        if (discrete == 0) {
          warning("Estimation cannot be performed from the likelihood associated with a continuous observation of the trajectories.")
        }
        
        res <- EstParamGamma(U, V, S, SigDelta, K, drift.fixed, estim.drift.fix)
        
        bic <- res$BIChere
        aic <- res$AIChere
        mu <- res$mu
        a <- res$a
        lambda <- res$lambda
        
        
        # Estimator of the density
        
        estimphi <- 0
        
        if (is.null(gridg) == 1) {
          gridg <- seq(min(estimpsi2) * 0.8, max(estimpsi2) * 1.2, length = 500)
        }
        
        if (is.null(gridg) == 0) {
          gridg <- gridg
        }
        
        simupsi2 <- 1/rgamma(500, a, rate = 1/lambda)
        
        testpsi <- density(simupsi2, from = min(simupsi2), to = max(simupsi2), bw = 'ucv', n = 500)
        if (testpsi$bw < 0.1) {
          testpsi <- density(simupsi2, from = min(simupsi2), to = max(simupsi2), n = 500)
        }
        
        gridg <- seq(min(estimpsi2) * 0.8, max(estimpsi2) * 1.2, length = 500)
        estimg <- matrix(testpsi$y, nrow = 1)
        gridg <- matrix(gridg, nrow = 1)
        
        estimpsi2 <- matrix(estimpsi2, 1, length(estimpsi2), byrow = TRUE)
        ## ? troncation ??? estimpsi2.trunc <- estimpsi2 estimf.trunc <- estimf
        cutoff <- FALSE
        gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
        estimf <- matrix(estimf, 1, length(estimf), byrow = TRUE)
        estimphi <- matrix(estimphi, 1, length(estimphi), byrow = TRUE)
      }
      
      if (((sum(drift.random)) >= 1) && (diffusion.random == 1)) {
        
        # -- L'estimation de la densité va nécessiter de simuler des réalisation des effets
        # aléatoires dans la dérive
        
        res <- EstParamNormalGamma(U, V, S, SigDelta, K, drift.random, drift.fixed, 
                                   estim.drift.fix)
        bic <- res$BIChere
        aic <- res$AIChere
        mu <- res$mu
        omega <- res$omega
        a <- res$a
        lambda <- res$lambda
        
        # -- Estimator of the density
        
        estimpsi2 <- matrix(estimpsi2, 1, length(estimpsi2), byrow = TRUE)
        
        simupsi2 <- 1/rgamma(500, a, rate = 1/lambda)
        testpsi <- density(simupsi2, from = min(simupsi2), to = max(simupsi2), bw = 'ucv', n = 500)
        if (testpsi$bw < 0.1) {
          testpsi <- density(simupsi2, from = min(simupsi2), to = max(simupsi2), n = 500)
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
          
          test1 <- density(simuphi[1, ], from = min(simuphi[1, ]), to = max(simuphi[1, ]), bw = 'ucv', n = 500)
          test2 <- density(simuphi[2, ], from = min(simuphi[2, ]), to = max(simuphi[2, ]), bw = 'ucv', n = 500)
          
          if (test1$bw < 0.1) {
            test1 <- density(simuphi[1, ], from = min(simuphi[1, ]), to = max(simuphi[1, ]), n = 500)
          }
          estimf1 <- test1$y
          gridf1 <- test1$x
          
          if (test2$bw < 0.1) {
            test2 <- density(simuphi[2, ], from = min(simuphi[2, ]), to = max(simuphi[2, ]), n = 500)
          }
          estimf2 <- test2$y
          gridf2 <- test2$x
          
          gridf <- matrix(0, 2, 500)
          gridf[1, ] <- gridf1
          gridf[2, ] <- gridf2
          
          estimf <- matrix(0, 2, 500)
          estimf[1, ] <- estimf1
          estimf[2, ] <- estimf2
          
          estimf.trunc <- estimf
          estimphi.trunc <- estimphi
          
          cutoff <- FALSE
        }
        
        if (sum(drift.random) == 2) {
          simuphi <- rep(NA, length(simupsi2))
          simuphi <- sapply(1:length(simupsi2), function(s) {
            rnorm(1, mean = mu[2], sd = abs(omega[2] * sqrt(simupsi2[s])))
          })
          
          
          test <- density(simuphi, from = min(simuphi), to = max(simuphi), bw = 'ucv', n = 500)
          if (test$bw < 0.1) {
            test <- density(simuphi, from = min(simuphi), to = max(simuphi), n = 500)
          }
          
          gridf <- test$x
          estimf <- test$y
          
          estimf <- matrix(estimf, nrow = 1)
          gridf <- matrix(gridf, nrow = 1)
          
          estimf.trunc <- estimf
          estimphi.trunc <- estimphi
          
          cutoff <- FALSE
          
        }
        
        if (sum(drift.random) == 1) {
          
          simuphi <- rep(NA, length(simupsi2))
          simuphi <- sapply(1:length(simupsi2), function(s) {
            rnorm(1, mean = mu[1], sd = abs(omega[1] * sqrt(simupsi2[s])))
          })
          
          test <- density(simuphi, from = min(simuphi), to = max(simuphi), n = length(simuphi), bw = "ucv")
          
          if (test$bw < 0.1) {
            test <- density(simuphi, from = min(simuphi), to = max(simuphi), n = length(simuphi))
          }
          
          
          gridf <- test$x
          estimf <- test$y
          
          
          
          estimf <- matrix(estimf, nrow = 1)
          gridf <- matrix(gridf, nrow = 1)
          
          estimf.trunc <- estimf
          estimphi.trunc <- estimphi
          
          cutoff <- FALSE
          
        }
        
      }
      
    }
    
    
    return(new(Class = "Freq.fit", model = model, drift.random = drift.random, diffusion.random = diffusion.random, 
               gridf = gridf, mu = mu, omega = omega, a = a, lambda = lambda, cutoff = cutoff, 
               sigma2 = sigma2, index = index, estimphi = estimphi, estimpsi2 = estimpsi2, 
               estimf = estimf, estim.drift.fix = estim.drift.fix, estim.diffusion.fix = estim.diffusion.fix, 
               discrete = discrete, bic = bic, aic = aic, times = times, X = X, gridg = gridg, 
               estimg = estimg))
  }
  
  
  if (mixture == 1) {
    
    ## 
    cutoff <- FALSE
    gridf <- NULL
    estimf <- NULL
    ## 
    
    
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
        warning("All the trajectories have non positive values the model CIR cannot be used", 
                call. = FALSE)
        estimf <- 0
        estimphi <- 0
        estimpsi2 <- 0
        bic <- 0
        aic <- 0
        gridf <- 0
        mu <- 0
        omega <- 0
        cutoff <- 0
        sigma2 <- 0
        a <- 0
        lambda <- 0
        estimf.trunc <- 0
        estimphi.trunc <- 0
        estimpsi2.trunc <- 0
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
      muinit[,1] <- mean(estimphi[1,])
      
    }
    
    if (sum(drift.random) == 1) {
      for (n in 1:nb.mixt) {
        label <- which(km$cluster == n)
        muinit[n, 1] <- mean(estimphi[1, label])
        omegainit[n, 1] <- sd(estimphi[1, label])
      }
      muinit[,2] <- mean(estimphi[2,])
    }
    
    probinit <- table(km$cluster)/Mindex2
    
    start <- list(mu = muinit, omega = omegainit, mixt.prop = probinit)
    
    # -- Parameter estimation
    
    res <- EM(U, V, S, K, drift.random, start, Niter, drift.fixed, estim.drift.fix, 
              drift.fixed.mixed, diffusion.fixed, estim.diffusion.fix)
    
    bic <- res$BIChere
    aic <- res$AIChere
    mu <- res$mu
    omega <- res$omega
    sigma2 <- res$sigma^2
    mixt.prop <- res$mixt.prop
    probindi <- res$probindi
    
    # -- Estimator of the density
    
    if (sum(drift.random) == 3) {
      
      if (is.null(gridf) == 1) {
        gridf <- matrix(0, 2, 500)
        bg1 <- max(abs(estimphi[1,])) * 1.2
        bg2 <- max(abs(estimphi[2,])) * 1.2
        gridf[1, ] <- seq(min(abs(estimphi[1,])) * 0.8, max(abs(estimphi[1,])) * 1.2, length = 500)
        gridf[2, ] <- seq(min(abs(estimphi[2,])) * 0.8, max(abs(estimphi[2,])) * 1.2, length = 500)
      }
      
      if (is.null(gridf) == 0) {
        gridf <- gridf
      }
      
      estimf1 <- rep(0, length(gridf[1, ]), byrow = TRUE)
      estimf2 <- rep(0, length(gridf[2, ]), byrow = TRUE)
      
      for (n in 1:nb.mixt) {
        estimf1 <- estimf1 + mixt.prop[n] * dnorm(gridf[1, ], mean = mu[Niter, 
                                                                        n, 1], sd = abs(omega[Niter, n, 1]))
        estimf2 <- estimf2 + mixt.prop[n] * dnorm(gridf[2, ], mean = mu[Niter, 
                                                                        n, 2], sd = abs(omega[Niter, n, 2]))
      }
      #estimf <- estimf1 %*% t(estimf2)
      estimf <- rbind(estimf1,estimf2)
      #estimf.trunc <- estimf
      #estimphi.trunc <- estimphi
      
      cutoff <- FALSE
    }
    
    if (sum(drift.random) == 1 ) {
      
      if (is.null(gridf) == 1) {
        gridf <- matrix(0, 2, 500)
        gridf <- seq(min(abs(estimphi[1,])) * 0.8, max(abs(estimphi[1,])) * 1.2, length = 500)
      }
      
      if (is.null(gridf) == 0) {
        gridf <- gridf
      }
      
      estimf <- rep(0, length(gridf), byrow = TRUE)
      
      for (n in 1:nb.mixt) {
        estimf <- estimf + mixt.prop[n] * dnorm(gridf, mean = mu[Niter, n, 1], sd = abs(omega[Niter, n, 1]))
      }
      
      gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
      estimf <- matrix(estimf, 1, length(estimf), byrow = TRUE)
      cutoff <- FALSE
    }
    
    if (sum(drift.random) == 2 ) {
      
      if (is.null(gridf) == 1) {
        gridf <- matrix(0, 2, 500)
        gridf <- seq(min(abs(estimphi[2,])) * 0.8, max(abs(estimphi[2,])) * 1.2, length = 500)
      }
      
      if (is.null(gridf) == 0) {
        gridf <- gridf
      }
      
      estimf <- rep(0, length(gridf), byrow = TRUE)
      
      for (n in 1:nb.mixt) {
        estimf <- estimf + mixt.prop[n] * dnorm(gridf, mean = mu[Niter, n, 2], sd = abs(omega[Niter, n, 2]))
      }
      
      gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
      estimf <- matrix(estimf, 1, length(estimf), byrow = TRUE)
      cutoff <- FALSE
    }
    
    return(new(Class = "Freq.mixture.fit", model = model, drift.random = drift.random, 
               gridf = gridf, mu = mu, omega = omega, mixt.prop = mixt.prop, cutoff = cutoff, 
               sigma2 = sigma2, index = index, probindi = probindi, estimf = estimf, estimphi = estimphi, 
               bic = bic, aic = aic, estim.drift.fix = estim.drift.fix, times = times, X = X))
    
    
  }
  
}


# S4 class for the parametric estimation results -------------------------


#' S4 class for the parametric estimation results  
#'  
#' @slot model character 'OU' or 'CIR'
#' @slot drift.random numeric 1, 2, or c(1,2)
#' @slot diffusion.random numeric 0 or 1
#' @slot gridf matrix of values on which the estimation of the density of the random effects in the drift is done
#' @slot gridg matrix of values on which the estimation of the density of the random effects in the diffusion is done
#' @slot mu numeric MLE estimator for parametric approach
#' @slot omega numeric  MLE estimator for parametric approach
#' @slot a numeric MLE estimator for parametric approach
#' @slot lambda numeric MLE estimator for parametric approach
#' @slot cutoff value of the cutoff if there is one
#' @slot sigma2 numeric estimated value of \eqn{\sigma^2} if the diffusion coefficient is not random
#' @slot index index of the used trajectories
#' @slot estimphi matrix of the estimator of the drift random effects
#' @slot estimpsi2 vector of the estimator of the diffusion random effects
#' @slot estimf estimator of the (conditional) density of \eqn{\phi}, matrix form
#' @slot estimg estimator of the density of \eqn{\phi}, matrix form
#' @slot estim.drift.fix 1 if the user asked for the estimation of fixed parameter in the drift
#' @slot estim.diffusion.fix 1 if the user asked for the estimation of fixed diffusion coefficient
#' @slot discrete 1 if the estimation is based on the likelihood of discrete observations, 0 otherwise 
#' @slot bic numeric bic 
#' @slot aic numeric aic
#' @slot times vector of observation times, storage of input variable
#' @slot X matrix of observations, storage of input variable


setClass(Class = "Freq.fit", representation = representation(model = "character", drift.random = "numeric", 
                                                             diffusion.random = "numeric", gridf = "matrix", gridg = "matrix", mu = "numeric", omega = "numeric", 
                                                             a = "numeric", lambda = "numeric", cutoff = "logical", sigma2 = "numeric", index = "numeric", 
                                                             estimphi = "matrix", estimpsi2 = "matrix", estimf = "matrix", estimg = "matrix", estim.drift.fix = "numeric", 
                                                             estim.diffusion.fix = "numeric", discrete = "numeric", bic = "numeric", aic = "numeric", 
                                                             times = "numeric", X = "matrix"))

# S4 class for the parametric estimation results (mixture) -------------------------


#' S4 class for the parametric estimation results when the random effects in the drift follow 
#' mixture of normal distributions  
#'  
#' @slot model character 'OU' or 'CIR'
#' @slot drift.random numeric 1, 2, or c(1,2)
#' @slot gridf matrix of values on which the estimation of the density of the random effects is done
#' @slot mu array estimated value of the mean of phi at each iteration of the EM algorithm (Niter x nb.mixt x 2)
#' @slot omega array estimated value of the standard deviation of phi at each iteration of the EM algorithm (Niter x nb.mixt x 2)
#' @slot mixt.prop matrix estimated value of the mixing proportions at each iteration of the EM algorithm (Niter x nb.mixt)
#' @slot cutoff value of the cutoff if there is one
#' @slot sigma2 numeric estimated value of \eqn{\sigma^2} if the diffusion coefficient is not random
#' @slot index index of the used trajectories
#' @slot estimphi matrix of the estimator of the drift random effects
#' @slot probindi matrix of posterior component probabilities
#' @slot estimf matrix estimator of the density of \eqn{\phi} 
#' @slot estim.drift.fix numeric 1 if the user asked for the estimation of fixed parameter in the drift
#' @slot bic numeric bic 
#' @slot aic numeric aic
#' @slot times vector of observation times, storage of input variable
#' @slot X matrix of observations, storage of input variable


setClass(Class = "Freq.mixture.fit", representation = representation(model = "character", 
                                                                     drift.random = "numeric", gridf = "matrix", mu = "array", omega = "array", mixt.prop = "matrix", 
                                                                     cutoff = "logical", sigma2 = "numeric", index = "numeric", estimphi = "matrix", probindi = "matrix", 
                                                                     estimf = "matrix", estim.drift.fix = "numeric", bic = "numeric", aic = "numeric", times = "numeric", 
                                                                     X = "matrix"))


# Outputs (fonction out) --------------------------------------------------


########################################################### OUTPUTS
#' Transfers the class object to a list
#' 
#' @description Method for the S4 classes
#' @param x Freq.fit or Freq.mixture.fit class


out <- function(x) {
  sN <- slotNames(x)
  res <- lapply(sN, function(name) slot(x, name))
  names(res) <- sN
  res
}



# Summary for Freq.fit
# -----------------------------------------------------------------



############################################################### SUMMARY

#' Short summary of the results of class object Freq.fit
#' @description Method for the S4 class Freq.fit
#' @param object Freq.fit class
#'
#'
setMethod(f = "summary", signature = "Freq.fit", definition = function(object) {
  
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
    
    diffreff <- (matrix(c(object@a, object@lambda, object@a * object@lambda, digamma(object@a) + 
                            log(object@lambda)), 4, 1, byrow = TRUE))
    
    rownames(diffreff) <- c("a (shape)", "lambda (scale)", "m = a*lambda", "t = psi(a) + log(lambda)")
    
    print(diffreff, quote = FALSE, right = TRUE)
    
  }
  
  if (object@diffusion.random == 0) {
    
    cat("\nFixed effect in the diffusion coefficient:\n")
    
    if ((object@estim.diffusion.fix == 1) && (object@discrete == 1)) {
      print(matrix(c("MLE sigma", round(sqrt(object@sigma2), 6)), 1, 2, byrow = TRUE), 
            quote = FALSE, right = TRUE)
    }
    
    if ((object@estim.diffusion.fix == 1) && (object@discrete == 0)) {
      print(matrix(c("Estimated sigma (quadratic variations)", round(sqrt(object@sigma2), 
                                                                     6)), 1, 2, byrow = TRUE), quote = FALSE, right = TRUE)
    }
    
    if (object@estim.diffusion.fix == 0) {
      print(matrix(c("sigma (not estimated)", round(sqrt(object@sigma2), 6)), 1, 
                   2, byrow = TRUE), quote = FALSE, right = TRUE)
    }
    
    
  }
  
  
  cat("\nCriteria for model selection:\n")
  
  
  info.criteria <- matrix(c(round(object@bic, 6), round(object@aic, 6)), 2, 1, byrow = TRUE)
  rownames(info.criteria) <- c("BIC", "AIC")
  colnames(info.criteria) <- c("")
  print(info.criteria, quote = FALSE, right = TRUE)
  
  
})


# Summary pour Freq.mixture.fit -------------------------------------------



############################################################### SUMMARY

#' Short summary of the results of class object Freq.mixture.fit
#' @description Method for the S4 class Freq.fit
#' @param object Freq.mixture.fit class
#'
setMethod("summary", "Freq.mixture.fit", function(object) {
  
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
      
      driftreff <- (matrix(c(object@mu[Niter, , 1], object@omega[Niter, , 1], object@mu[Niter, 
                                                                                        , 2], object@omega[Niter, , 2]), 4, nb.mixt, byrow = TRUE))
      
      rownames(driftreff) <- c("MLE mean 1", "MLE sd 1", "MLE mean 2", "MLE sd 2")
      colnames(driftreff) <- col.name
      
      print(driftreff, quote = FALSE, right = TRUE)
      
      
      
    }
    
    
    if (sum(object@drift.random) == 2) {
      cat("\nFixed and random effects in the drift:\n")
      
      driftreff <- (matrix(c(object@mu[Niter, , 1], object@mu[Niter, , 2], object@omega[Niter, 
                                                                                        , 2]), 3, nb.mixt, byrow = TRUE))
      
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
      
      driftreff <- (matrix(c(object@mu[Niter, , 1], object@omega[Niter, , 1], object@mu[Niter, 
                                                                                        , 2]), 3, nb.mixt, byrow = TRUE))
      
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


########################################################### PRINT


# Print for Freq.fit ------------------------------------------------------


#' Description of print
#' @description Method for the S4 class Freq.fit
#' @param x Freq.fit class
#'
setMethod("print", "Freq.fit", function(x) {
  if (x@bic != 0) {
    
    print(c("number of used trajectories", length(x@index)))
  }
  
  if (x@bic == 0) {
    
    if (sum(x@cutoff) != 0) {
      print(matrix(c("number of used trajectories", length(x@index), "number of truncated values", 
                     length(x@cutoff) - sum(x@cutoff)), 2, 2, byrow = TRUE))
    }
    if (sum(x@cutoff) == 0) {
      print(c("number of used trajectories", length(x@index)))
    }
  }
})

#' Description of print
#' @description Method for the S4 class Freq.mixture.fit
#' @param x Freq.mixture.fit class
#'
setMethod("print", "Freq.mixture.fit", function(x) {
  if (x@bic != 0) {
    
    print(c("number of used trajectories", length(x@index)))
  }
  
  if (x@bic == 0) {
    
    if (sum(x@cutoff) != 0) {
      print(matrix(c("number of used trajectories", length(x@index), "number of truncated values", 
                     length(x@cutoff) - sum(x@cutoff)), 2, 2, byrow = TRUE))
    }
    if (sum(x@cutoff) == 0) {
      print(c("number of used trajectories", length(x@index)))
    }
  }
  
  print(c("Number of iterations of the EM algorithm", dim(x@mixt.prop)[1]))
})



# Plot for Freq.fit -------------------------------------------------------



########################################################### PLOT
#' Plot method for the frequentist estimation class object
#'
#' @description Plot method for the S4 class Freq.fit
#' @param x Freq.fit class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#'
#' @importFrom grDevices x11
#' @importFrom graphics par
#' @importFrom graphics plot

setMethod(f = "plot", signature = "Freq.fit", definition = function(x, newwindow = FALSE, 
                                                                    ...) {
  if (newwindow) {
    x11(width = 14)
  }
  
  
  if (x@diffusion.random == 0) {
    if (dim(x@gridf)[1] == 1) {
      
      op <- par(mfrow = c(1, 1), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), 
                oma = c(1, 1, 1, 1), omi = c(0.2, 0.2, 0.2, 0.2), 
                cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
      plot(x@gridf, x@estimf, main = "Estimated density of the random effect",xlab="Value",ylab="Density")
      
    }
    
    if (dim(x@gridf)[1] == 2) {
      
      op <- par(mfrow = c(2, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), 
                oma = c(1, 1, 1, 1), omi = c(0.2, 0.2, 0.2, 0.2), 
                cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
      plot(x@gridf[1,], x@estimf[1,], main='First random effect' ,xlab="Value",ylab="Density")
      plot(x@gridf[2,], x@estimf[2,], main='Second random effect',xlab="Value",ylab="Density")
      title('Estimated densities',outer=TRUE)
      
    }
    
  } else {
    
   # if (x@gridf==0) {
      # op <- par(mfrow = c(1, 1), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), 
      #           oma = c(1, 1, 1, 1), omi = c(0.2, 0.2, 0.2, 0.2), 
      #           cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
      # plot(x@gridg, x@estimg, main = "Estimated density of the random effect in the diffusion",xlab="Value",ylab="Density")
   #  }
    
    if (dim(x@gridf)[1] == 1) {
      
      op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), 
                oma = c(1, 1, 1, 1), omi = c(0.2, 0.2, 0.2, 0.2), 
                cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
      plot(x@gridf, x@estimf, main = "Random effect in the drift",xlab="Value",ylab="Density")
      plot(x@gridg, x@estimg, main = "Random effect in the diffusion",xlab="Value",ylab="Density")
      title("Estimated densities",outer=TRUE)
    }
    
    if (dim(x@gridf)[1] == 2) {
      
      op <- par(mfrow = c(2, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), 
                oma = c(1, 1, 1, 1), omi = c(0.2, 0.2, 0.2, 0.2), cex.main = 1, 
                cex.lab = 0.9, cex.axis = 0.9)
      plot(x@gridf[1,], x@estimf[1,], main='First drift random effect',xlab="Value",ylab="Density")
      plot(x@gridf[2,], x@estimf[2,], main='Second drift random effect',xlab="Value",ylab="Density")
      plot(x@gridg, x@estimg, main = "Random effect in the diffusion",xlab="Value",ylab="Density")
      title('Estimated densities',outer=TRUE)
      
    }
    
  }
  
})

# Plot for Freq.mixture.fit (A FAIRE) -------------------------------------------------------



########################################################### PLOT
#' Plot method for the frequentist estimation class object
#'
#' @description Plot method for the S4 class Freq.mixture.fit
#' @param x Freq.mixture.fit class
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param ... optional plot parameters
#'
#' @importFrom grDevices x11
#' @importFrom graphics par
#' @importFrom graphics title
#' @importFrom graphics layout

setMethod(f = 'plot', signature = 'Freq.mixture.fit', definition = function(x, newwindow = FALSE, ...) {
  if (newwindow) {
    x11(width = 14)
  }
  
  # Convergence plot of the EM algorithm
  
  
  mat = matrix(c(1,3,5,2,4,5), 3, 2)
  op <- par(mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), omi = c(1, 1, 1, 1),
            oma = c(1, 1, 1, 1), cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
  #mfcol = c(3, 2), 
  
  nb.mixt <- dim(x@mu)[2]
  for (k in 1:nb.mixt){
    #plot.new()
    layout(mat)
    tit <- paste("Convergence plot for mixture component ",k) 
    plot(x@mu[,k,1],type='l',xlab='Iteration',ylab='')
    title("mu1")
    plot(x@mu[,k,2],type='l',xlab='Iteration',ylab='')
    title('mu2')
    plot(x@omega[,k,1],type='l',xlab='Iteration',ylab='')
    title('omega1')
    plot(x@omega[,k,2],type='l',xlab='Iteration',ylab='')
    title('omega2') 
    plot(x@mixt.prop[,k],type='l',xlab='Iteration',ylab='')
    title('Proportion')
    title(tit , outer = TRUE)
  }
  
  # Plots of the estimated densities
  
  if (sum(x@drift.random)==3){
    op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), omi = c(1, 1, 1, 1),
              oma = c(1, 1, 1, 1), cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
    
    
    plot(x@gridf[1,],x@estimf[1,],ylab="Density",xlab="value")
    title('First random effect')
    plot(x@gridf[2,],x@estimf[2,],ylab="Density",xlab="value")
    title('Second random effect')
    title('Estimated densities', outer=TRUE)
  }
  
  if (sum(x@drift.random)<=2){
    op <- par(mfrow = c(1, 1), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), omi = c(1, 1, 1, 1),
              oma = c(1, 1, 1, 1), cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
    
    
    plot(x@gridf,x@estimf,ylab="Density",xlab="value")
    title('Estimated density', outer=TRUE)
  }
  
})

#' ########################################################### VALIDATION

#' Validation of the chosen model.
#'
#' @description Validation of the chosen model. For the index numj, Mrep=100 new trajectories are simulated
#' with the value of the estimated random effect number numj. Two plots are given: on the left the simulated trajectories and the true one (red)
#' and one the left the corresponding qq-plot for each time.
#' @param x Freq.fit or Freq.mixture.fit class
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
#' @param x Freq.fit class
#' @param Mrep number of trajectories to be drawn
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param plot.valid logical(1), if TRUE, the results are depicted grafically
#' @param numj optional number of series to be validated
#' @param ... optional plot parameters
#' @importFrom grDevices x11
#' @importFrom graphics par
#' @importFrom graphics abline
#' @references
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: a R package to fit mixed stochastic differential equations.
#'

setMethod(f = "valid", signature = "Freq.fit", definition = function(x, Mrep = 100, newwindow = FALSE, 
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
        suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), 
                                           delta = del, method = "EA", theta = c(phihat, sig), model = "OU", M = Mrep)))
      }
      if (x@model == "CIR") {
        suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), 
                                           delta = del, method = "milstein", theta = c(phihat, sig), model = "CIR", 
                                           M = Mrep, sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig * 
                                                                                                                   sqrt(x)))))
      }
    }
    
    if (x@diffusion.random == 1) {
      phihat <- x@estimphi[, numj]
      psihat <- sqrt(x@estimpsi2[numj])
      
      if (x@model == "OU") {
        suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), 
                                           delta = del, method = "EA", theta = c(phihat, psihat), model = "OU", 
                                           M = Mrep)))
      }
      if (x@model == "CIR") {
        suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), 
                                           delta = del, method = "milstein", theta = c(phihat, psihat), model = "CIR", 
                                           M = Mrep, sigma.x = expression(psihat/(2 * sqrt(x))), sigma = expression(psihat * 
                                                                                                                      sqrt(x)))))
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
          suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), 
                                             delta = del, method = "EA", theta = c(phihat, paramfixed, sig), model = "OU", 
                                             M = Mrep)))
        }
        
        if (x@model == "CIR") {
          suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), 
                                             delta = del, method = "milstein", theta = c(phihat, paramfixed, sig), 
                                             model = "CIR", M = Mrep, sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig * 
                                                                                                                                    sqrt(x)))))
        }
      }
      if (sum(x@drift.random) == 2) {
        paramfixed <- x@mu[1]
        if (x@model == "OU") {
          suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), 
                                             delta = del, method = "EA", theta = c(paramfixed, phihat, sig), model = "OU", 
                                             M = Mrep)))
        }
        if (x@model == "CIR") {
          suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), 
                                             delta = del, method = "milstein", theta = c(paramfixed, phihat, sig), 
                                             model = "CIR", M = Mrep, sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig * 
                                                                                                                                    sqrt(x)))))
        }
      }
    }
    
    if (x@diffusion.random == 1) {
      psihat <- sqrt(x@estimpsi2[numj])
      if (sum(x@drift.random) == 1) {
        paramfixed <- x@mu[2]
        if (x@model == "OU") {
          suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), 
                                             delta = del, method = "EA", theta = c(phihat, paramfixed, psihat), 
                                             model = "OU", M = Mrep)))
        }
        
        if (x@model == "CIR") {
          suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), 
                                             delta = del, method = "milstein", theta = c(phihat, paramfixed, psihat), 
                                             model = "CIR", M = Mrep, sigma.x = expression(psihat/(2 * sqrt(x))), 
                                             sigma = expression(psihat * sqrt(x)))))
        }
      }
      if (sum(x@drift.random) == 2) {
        paramfixed <- x@mu[1]
        if (x@model == "OU") {
          suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), 
                                             delta = del, method = "EA", theta = c(paramfixed, phihat, psihat), 
                                             model = "OU", M = Mrep)))
        }
        if (x@model == "CIR") {
          suppressMessages(Xnew <- t(sde.sim(T = Tend, X0 = Xtrue[numj, 1], N = length(timessimu), 
                                             delta = del, method = "milstein", theta = c(paramfixed, phihat, psihat), 
                                             model = "CIR", M = Mrep, sigma.x = expression(psihat/(2 * sqrt(x))), 
                                             sigma = expression(psihat * sqrt(x)))))
        }
      }
    }
  }
  
  vecttimes <- intersect(round(timessimu, 10), round(times, 10))
  
  N <- length(vecttimes)
  
  q <- rep(0, N)
  for (i in 1:N) {
    q[i] <- sum(Xtrue[numj, which(times == vecttimes[i])[1]] > Xnew[, which(timessimu == 
                                                                              vecttimes[i])[1]])/Mrep
  }
  
  if (plot.valid == 1) {
    numj <- numj
    op <- par(mfrow = c(1, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 
                                                                                      0, 0, 0), cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
    
    plot(c(0, timessimu), Xnew[1, ], type = "l", ylim = c(min(Xnew) * 0.8, max(Xnew) * 
                                                            1.2), xlab = "", ylab = "")
    for (k in 1:Mrep) {
      lines(c(0, timessimu), Xnew[k, ])
    }
    lines(times, Xtrue[numj, ], col = "red", lwd = 2)
    
    
    plot(1:N/N, sort(q), xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
    abline(0, 1)
  }
  
  return(list(quantiles = q, Xnew = Xnew, numj = numj))
})

########################################################### PREDICTION

#' Prediction method
#'
#' @description Prediction
#' @param x Freq.fit 
#' @param ... other optional parameters
#'
#' @importFrom graphics abline
#' @importFrom stats quantile
#' @importFrom sde sde.sim


setGeneric("pred", function(x, ...) {
  standardGeneric("pred")
})


#' Prediction method for the Freq.fit class object
#'
#' @description Frequentist prediction
#' @param x Freq.fit class
#' @param invariant 1 if the initial value is from the invariant distribution, default X0 is fixed from Xtrue
#' @param level alpha for the predicion intervals, default 0.05
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot
#' @param plot.pred logical(1), if TRUE, the results are depicted grafically
#' @param ... optional plot parameters
#' @importFrom grDevices x11
#' @importFrom graphics par
#' @references
#' Dion, C., Hermann, S. and Samson, A. (2016). Mixedsde: a R package to fit mixed stochastic differential equations.
#'
setMethod(f = "pred", signature = "Freq.fit", definition = function(x, invariant = 0, level = 0.05,
                                                                    newwindow = FALSE, plot.pred = TRUE, ...) {
  if (newwindow) {
    x11(width = 10)
  }
  
  Xtrue <- x@X
  timestrue <- x@times
  T <- timestrue[length(timestrue)]
  
  if (dim(x@gridf)[1] == 1) {
    if (x@diffusion.random == 0) {
      sig <- sqrt(x@sigma2)
      
      index <- x@index
      M <- length(index)
      N <- dim(Xtrue)[2] - 1
      delta <- T/N
      Xpred <- matrix(0, M, N + 1)
      times <- seq(0, T, by = delta)
      
      if (sum(x@drift.random) == 1) {
        paramfixed <- x@mu[2]
        
        phipred <- rep(0, M)
        
        if (x@bic == 0) {
          p <- x@estimf/sum(x@estimf)
          for (i in 1:M) {
            phipred[i] <- discr(x@gridf, p)
          }
        }
        if (x@bic != 0) {
          phipred <- rnorm(M, x@mu[1], x@omega[1])
        }
        
        if (x@model == "OU") {
          indexpred <- 1:M
          for (j in 1:M) {
            if (invariant == FALSE) {
              suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[j, 1], N = N,
                                                     delta = T/N, method = "EA", theta = c(phipred[j], paramfixed, sig),
                                                     model = "OU"))
            }
            if (invariant == TRUE) {
              X0 <- phipred[j]/paramfixed + (sig/(sqrt(2 * paramfixed))) * rnorm(1)
              suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N,
                                                     method = "EA", theta = c(phipred[j], paramfixed, sig), model = "OU"))
            }
          }
        }
        
        if (x@model == "CIR") {
          indexpred <- which(phipred > 0)
          phipred <- phipred[indexpred]
          Mpred <- length(phipred)
          Xpred <- matrix(0, Mpred, N + 1)
          
          for (j in 1:Mpred) {
            if (invariant == 0) {
              suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[indexpred[j],
                                                                       1], N = N, delta = T/N, method = "milstein", theta = c(phipred[j],
                                                                                                                              paramfixed, sig), model = "CIR", sigma.x = expression(sig/(2 *
                                                                                                                                                                                           sqrt(x))), sigma = expression(sig * sqrt(x))))
            }
            if (invariant == 1) {
              X0 <- rgamma(1, 2 * phipred[j]/sig^2, scale = sig^2/(2 * paramfixed))
              suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N,
                                                     method = "milstein", theta = c(phipred[j], paramfixed, sig), model = "CIR",
                                                     sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig *
                                                                                                                   sqrt(x))))
            }
          }
        }
      }
      
      if (sum(x@drift.random) == 2) {
        paramfixed <- x@mu[1]
        
        phipred <- rep(0, M)
        
        if (x@bic == 0) {
          p <- x@estimf/sum(x@estimf)
          for (i in 1:M) {
            phipred[i] <- discr(x@gridf, p)
          }
        }
        
        if (x@bic != 0) {
          phipred <- rnorm(M, x@mu[2], x@omega[2])
        }
        
        if (x@model == "OU") {
          indexpred <- which(phipred > 0)
          phipred <- phipred[indexpred]
          Mpred <- length(indexpred)
          Xpred <- matrix(0, Mpred, N + 1)
          
          for (j in 1:Mpred) {
            if (invariant == 0) {
              suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[j, 1], N = N,
                                                     delta = T/N, method = "EA", theta = c(paramfixed, phipred[j], sig),
                                                     model = "OU"))
            }
            if (invariant == 1) {
              X0 <- (sig/(sqrt(2 * phipred[j]))) * rnorm(1)
              suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N,
                                                     method = "EA", theta = c(paramfixed, phipred[j], sig), model = "OU"))
            }
          }
        }
        
        if (x@model == "CIR") {
          indexpred <- which(phipred > 0)
          phipred <- phipred[indexpred]
          Mpred <- length(phipred)
          Xpred <- matrix(0, Mpred, N + 1)
          
          for (j in 1:Mpred) {
            if (invariant == 1) {
              if (paramfixed == 0) {
                print("no invariant distribution, please fix X0")
              }
              if (paramfixed != 0) {
                X0 <- rgamma(1, 2 * paramfixed/sig^2, scale = sig^2/(2 * phipred[j]))
                suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N,
                                                       method = "milstein", theta = c(paramfixed, phipred[j], sig),
                                                       model = "CIR", sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig *
                                                                                                                                    sqrt(x))))
              }
            }
            if (invariant == 0) {
              suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[indexpred[j],
                                                                       1], N = N, delta = T/N, method = "milstein", theta = c(paramfixed,
                                                                                                                              phipred[j], sig), model = "CIR", sigma.x = expression(sig/(2 *
                                                                                                                                                                                           sqrt(x))), sigma = expression(sig * sqrt(x))))
            }
          }
        }
      }
      
      if (plot.pred == TRUE) {
        op <- par(mfrow = c(1, 3), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0),
                  oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
        
        ## MODIF SUR x@estimphi
        plot(sort(x@estimphi[indexpred]), sort(phipred), pch = 18, xlim = c(min(x@estimphi[indexpred],
                                                                                phipred) * 0.8, max(x@estimphi[indexpred], phipred) * 1.2), ylim = c(min(x@estimphi[indexpred],
                                                                                                                                                         phipred) * 0.8, max(x@estimphi[indexpred], phipred) * 1.2), ylab = "",
             xlab = "", main = "Random effect in the drift")
        abline(0, 1)
        ## FIN MODIF SUR x@estimphi
        
        plot(timestrue, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "",
             ylim = c(min(Xtrue, Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
        for (j in indexpred) {
          lines(timestrue, Xtrue[j, ], col = j)
        }
        
        plot(times, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue,
                                                                               Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), main = "Predictive trajectories")
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
      
      return(list(phipred = phipred, Xpred = Xpred, indexpred = indexpred))
      
    }
    
    ## NOUVEAU, CAS D'UN EFFET ALEATOIRE DANS LA DIFFUSION
    if (x@diffusion.random == 1) {
      
      index <- x@index
      M <- length(index)
      N <- dim(Xtrue)[2] - 1
      delta <- T/N
      Xpred <- matrix(0, M, N + 1)
      times <- seq(0, T, by = delta)
      
      psipred <- rep(0, M)
      psipred <- 1/sqrt(rgamma(M, shape = x@a, rate = 1/x@lambda))
      
      if (sum(x@drift.random) == 1) {
        
        paramfixed <- x@mu[2]
        
        phipred <- rep(0, M)
        
        if (x@bic == 0) {
          p <- x@estimf/sum(x@estimf)
          for (i in 1:M) {
            phipred[i] <- discr(x@gridf, p)
          }
        }
        
        if (x@bic != 0) {
          for (j in 1:M) {
            phipred[j] <- rnorm(1, mean = x@mu[1], sd = x@omega[1] * psipred[j])
          }
        }
        
        if (x@model == "OU") {
          indexpred <- 1:M
          for (j in 1:M) {
            if (invariant == FALSE) {
              suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[j, 1], N = N,
                                                     delta = T/N, method = "EA", theta = c(phipred[j], paramfixed, psipred[j]),
                                                     model = "OU"))
            }
            if (invariant == TRUE) {
              X0 <- phipred[j]/paramfixed + (psipred[j]/(sqrt(2 * paramfixed))) *
                rnorm(1)
              suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N,
                                                     method = "EA", theta = c(phipred[j], paramfixed, psipred[j]), model = "OU"))
            }
          }
        }
        if (x@model == "CIR") {
          indexpred <- which(phipred > 0)
          phipred <- phipred[indexpred]
          psipred <- psipred[indexpred]
          Mpred <- length(phipred)
          Xpred <- matrix(0, Mpred, N + 1)
          
          for (j in 1:Mpred) {
            if (invariant == 0) {
              suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[indexpred[j],
                                                                       1], N = N, delta = T/N, method = "milstein", theta = c(phipred[j],
                                                                                                                              paramfixed, psipred[j]), model = "CIR", sigma.x = expression(psipred[j]/(2 *
                                                                                                                                                                                                         sqrt(x))), sigma = expression(psipred[j] * sqrt(x))))
            }
            if (invariant == 1) {
              X0 <- rgamma(1, 2 * phipred[j]/psipred[j]^2, scale = psipred[j]^2/(2 *
                                                                                   paramfixed))
              suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N,
                                                     method = "milstein", theta = c(phipred[j], paramfixed, psipred[j]),
                                                     model = "CIR", sigma.x = expression(psipred[j]/(2 * sqrt(x))),
                                                     sigma = expression(psipred[j] * sqrt(x))))
            }
          }
        }
      }
      
      if (sum(x@drift.random) == 2) {
        paramfixed <- x@mu[1]
        
        phipred <- rep(0, M)
        
        if (x@bic == 0) {
          p <- x@estimf/sum(x@estimf)
          for (i in 1:M) {
            phipred[i] <- discr(x@gridf, p)
          }
        }
        
        if (x@bic != 0) {
          for (j in 1:M) {
            phipred <- rnorm(M, x@mu[2], x@omega[2] * psipred[j])
          }
        }
        
        if (x@model == "OU") {
          indexpred <- which(phipred > 0)
          phipred <- phipred[indexpred]
          psipred <- psipred[indexpred]
          Mpred <- length(indexpred)
          Xpred <- matrix(0, Mpred, N + 1)
          
          for (j in 1:Mpred) {
            if (invariant == 0) {
              suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[j, 1], N = N,
                                                     delta = T/N, method = "EA", theta = c(paramfixed, phipred[j], psipred[j]),
                                                     model = "OU"))
            }
            if (invariant == 1) {
              X0 <- (sig/(sqrt(2 * phipred[j]))) * rnorm(1)
              suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N,
                                                     method = "EA", theta = c(paramfixed, phipred[j], psipred[j]), model = "OU"))
            }
          }
        }
        
        if (x@model == "CIR") {
          indexpred <- which(phipred > 0)
          phipred <- phipred[indexpred]
          psipred <- psipred[indexpred]
          Mpred <- length(phipred)
          Xpred <- matrix(0, Mpred, N + 1)
          
          for (j in 1:Mpred) {
            if (invariant == 1) {
              if (paramfixed == 0) {
                print("no invariant distribution, please fix X0")
              }
              if (paramfixed != 0) {
                X0 <- rgamma(1, 2 * paramfixed/psipred[j]^2, scale = psipred[j]^2/(2 *
                                                                                     phipred[j]))
                suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N,
                                                       method = "milstein", theta = c(paramfixed, phipred[j], psipred[j]),
                                                       model = "CIR", sigma.x = expression(psipred[j]/(2 * sqrt(x))),
                                                       sigma = expression(psipred[j] * sqrt(x))))
              }
            }
            if (invariant == 0) {
              suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[indexpred[j],
                                                                       1], N = N, delta = T/N, method = "milstein", theta = c(paramfixed,
                                                                                                                              phipred[j], psipred[j]), model = "CIR", sigma.x = expression(psipred[j]/(2 *
                                                                                                                                                                                                         sqrt(x))), sigma = expression(psipred[j] * sqrt(x))))
            }
          }
        }
      }
      
      if (plot.pred == TRUE) {
        op <- par(mfrow = c(2, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0),
                  oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
        
        plot(sort(x@estimphi[indexpred]), sort(phipred), pch = 18, xlim = c(min(x@estimphi[indexpred],
                                                                                phipred) * 0.8, max(x@estimphi[indexpred], phipred) * 1.2), ylim = c(min(x@estimphi[indexpred],
                                                                                                                                                         phipred) * 0.8, max(x@estimphi[indexpred], phipred) * 1.2), ylab = "",
             xlab = "", main = "Random effect in the drift")
        abline(0, 1)
        
        plot(sort(sqrt(x@estimpsi2[indexpred])), sort(psipred), pch = 18, xlim = c(min(sqrt(x@estimpsi2[indexpred]),
                                                                                       psipred) * 0.8, max(sqrt(x@estimpsi2[indexpred]), psipred) * 1.2), ylim = c(min(sqrt(x@estimpsi2[indexpred]),
                                                                                                                                                                       psipred) * 0.8, max(sqrt(x@estimpsi2[indexpred]), psipred) * 1.2), ylab = "",
             xlab = "", main = "Random effect in the diffusion")
        abline(0, 1)
        
        plot(timestrue, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "",
             ylim = c(min(Xtrue, Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
        for (j in indexpred) {
          lines(timestrue, Xtrue[j, ], col = j)
        }
        
        plot(times, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue,
                                                                               Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), main = "Predictive trajectories")
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
      
      return(list(phipred = phipred, psipred = psipred, Xpred = Xpred, indexpred = indexpred))
      
    }
    ## FIN NOUVEAU }
    
    
  }
  
  if (dim(x@gridf)[1] == 2) {
    
    
    index <- x@index
    M <- length(index)
    N <- dim(Xtrue)[2] - 1
    delta <- T/N
    Xpred <- matrix(0, M, N + 1)
    times <- seq(0, T, by = delta)
    
    phipred <- matrix(0, 2, M)
    
    if (x@diffusion.random == 0) {
      
      sig <- sqrt(x@sigma2)
      
      if (x@bic != 0) {
        mu <- x@mu
        omega <- x@omega
        mu1 <- mu[1]
        mu2 <- mu[2]
        omega1 <- omega[1]
        omega2 <- omega[2]
        phipred[1, ] <- rnorm(M, mu1, omega1)
        phipred[2, ] <- rnorm(M, mu2, omega2)
      }
      
      if (x@bic == 0) {
        gridf1 <- x@gridf[1, ]
        gridf2 <- x@gridf[2, ]
        
        
        marg1 <- ((max(gridf2) - min(gridf2))/length(gridf2)) * apply(x@estimf,
                                                                      1, sum)
        marg2 <- ((max(gridf1) - min(gridf1))/length(gridf1)) * apply(x@estimf,
                                                                      2, sum)
        p1 <- marg1/sum(marg1)
        p2 <- marg2/sum(marg2)
        phipred <- matrix(0, 2, sum(M))
        for (i in 1:M) {
          phipred[1, i] <- discr(gridf1, p1)
          phipred[2, i] <- discr(gridf2, p2)
        }
        
      }
      
      if (x@model == "OU") {
        indexpred <- which(phipred[2, ] > 0)
        phipred <- phipred[, indexpred]
        Mpred <- length(indexpred)
        Xpred <- matrix(0, Mpred, N + 1)
        for (j in 1:Mpred) {
          if (invariant == 0) {
            suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[j, 1], N = N,
                                                   delta = T/N, method = "EA", theta = c(phipred[, j], sig), model = "OU"))
          }
          if (invariant == 1) {
            X0 <- phipred[1, j]/phipred[2, j] + (sig/(sqrt(2 * phipred[2, j]))) *
              rnorm(1)
            suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N,
                                                   method = "EA", theta = c(phipred[, j], sig), model = "OU"))
          }
        }
      }
      if (x@model == "CIR") {
        indexpred <- which((phipred[1, ] > 0) & (phipred[2, ] > 0))
        phipred <- phipred[, indexpred]
        Mpred <- length(indexpred)
        Xpred <- matrix(0, Mpred, N + 1)
        
        for (j in 1:Mpred) {
          if (invariant == 0) {
            suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[indexpred[j],
                                                                     1], N = N, delta = T/N, method = "milstein", theta = c(phipred[,
                                                                                                                                    j], sig), model = "CIR", sigma.x = expression(sig/(2 * sqrt(x))),
                                                   sigma = expression(sig * sqrt(x))))
          }
          if (invariant == 1) {
            X0 <- rgamma(1, 2 * phipred[1, j]/sig^2, scale = sig^2/(2 * phipred[2,
                                                                                j]))
            suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N,
                                                   method = "milstein", theta = c(phipred[, j], sig), model = "CIR",
                                                   sigma.x = expression(sig/(2 * sqrt(x))), sigma = expression(sig *
                                                                                                                 sqrt(x))))
          }
        }
      }
      
      if (plot.pred == TRUE) {
        
        op <- par(mfrow = c(2, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0),
                  oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
        
        
        plot(sort(x@estimphi[1, indexpred]), sort(phipred[1, ]), pch = 18, xlim = c(min(x@estimphi[1,
                                                                                                   ], phipred[1, ]) * 0.8, max(x@estimphi[1, ], phipred[1, ]) * 1.2), ylim = c(min(x@estimphi[1,
                                                                                                                                                                                              ], phipred[1, ]) * 0.8, max(x@estimphi[1, ], phipred[1, ]) * 1.2), ylab = "",
             xlab = "", main = "First random effect in the drift")
        abline(0, 1)
        plot(sort(x@estimphi[2, indexpred]), sort(phipred[2, ]), pch = 18, xlim = c(min(x@estimphi[2,
                                                                                                   ], phipred[2, ]) * 0.8, max(x@estimphi[2, ], phipred[2, ]) * 1.2), ylim = c(min(x@estimphi[2,
                                                                                                                                                                                              ], phipred[2, ]) * 0.8, max(x@estimphi[2, ], phipred[2, ]) * 1.2), ylab = "",
             xlab = "", main = "Second random effect in the drift")
        abline(0, 1)
        
        
        plot(times, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue,
                                                                                          Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
        for (j in indexpred) {
          lines(times, Xtrue[j, ], col = j)
        }
        
        plot(timestrue, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue,
                                                                                   Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), main = "Predictive trajectories")
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
      
      return(list(phipred = phipred, Xpred = Xpred, indexpred = indexpred))
    }
    
    if (x@diffusion.random == 1) {
      
      if (x@bic != 0) {
        mu <- x@mu
        omega <- x@omega
        mu1 <- mu[1]
        mu2 <- mu[2]
        omega1 <- omega[1]
        omega2 <- omega[2]
        psipred <- rep(0, M)
        psipred <- 1/sqrt(rgamma(M, shape = x@a, rate = 1/x@lambda))
        for (j in 1:M) {
          phipred[1, j] <- rnorm(1, mu1, omega1 * psipred[j])
          phipred[2, j] <- rnorm(1, mu2, omega2 * psipred[j])
        }
      }
      
      if (x@bic == 0) {
        gridf1 <- x@gridf[1, ]
        gridf2 <- x@gridf[2, ]
        
        
        marg1 <- ((max(gridf2) - min(gridf2))/length(gridf2)) * apply(x@estimf,
                                                                      1, sum)
        marg2 <- ((max(gridf1) - min(gridf1))/length(gridf1)) * apply(x@estimf,
                                                                      2, sum)
        p1 <- marg1/sum(marg1)
        p2 <- marg2/sum(marg2)
        phipred <- matrix(0, 2, sum(M))
        for (i in 1:M) {
          phipred[1, i] <- discr(gridf1, p1)
          phipred[2, i] <- discr(gridf2, p2)
        }
        
      }
      
      if (x@model == "OU") {
        indexpred <- which(phipred[2, ] > 0)
        phipred <- phipred[, indexpred]
        Mpred <- length(indexpred)
        Xpred <- matrix(0, Mpred, N + 1)
        for (j in 1:Mpred) {
          if (invariant == 0) {
            suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[j, 1], N = N,
                                                   delta = T/N, method = "EA", theta = c(phipred[, j], psipred[j]),
                                                   model = "OU"))
          }
          if (invariant == 1) {
            X0 <- phipred[1, j]/phipred[2, j] + (psipred[j]/(sqrt(2 * phipred[2,
                                                                              j]))) * rnorm(1)
            suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N,
                                                   method = "EA", theta = c(phipred[, j], psipred[j]), model = "OU"))
          }
        }
      }
      if (x@model == "CIR") {
        indexpred <- which((phipred[1, ] > 0) & (phipred[2, ] > 0))
        phipred <- phipred[, indexpred]
        Mpred <- length(indexpred)
        Xpred <- matrix(0, Mpred, N + 1)
        
        for (j in 1:Mpred) {
          if (invariant == 0) {
            suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = Xtrue[indexpred[j],
                                                                     1], N = N, delta = T/N, method = "milstein", theta = c(phipred[,
                                                                                                                                    j], psipred[j]), model = "CIR", sigma.x = expression(psipred[j]/(2 *
                                                                                                                                                                                                       sqrt(x))), sigma = expression(psipred[j] * sqrt(x))))
          }
          if (invariant == 1) {
            X0 <- rgamma(1, 2 * phipred[1, j]/psipred[j]^2, scale = psipred[j]^2/(2 *
                                                                                    phipred[2, j]))
            suppressMessages(Xpred[j, ] <- sde.sim(T = T, X0 = X0, N = N, delta = T/N,
                                                   method = "milstein", theta = c(phipred[, j], psipred[j]), model = "CIR",
                                                   sigma.x = expression(psipred[j]/(2 * sqrt(x))), sigma = expression(psipred[j] *
                                                                                                                        sqrt(x))))
          }
        }
      }
      
      if (plot.pred == TRUE) {
        
        op <- par(mfrow = c(3, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0),
                  oma = c(0, 0, 0, 0), cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
        
        
        plot(sort(x@estimphi[1, indexpred]), sort(phipred[1, ]/psipred), pch = 18,
             xlim = c(min(x@estimphi[1, ], phipred[1, ]) * 0.8, max(x@estimphi[1,
                                                                               ], phipred[1, ]) * 1.2), ylim = c(min(x@estimphi[1, ], phipred[1, ]) *
                                                                                                                   0.8, max(x@estimphi[1, ], phipred[1, ]) * 1.2), ylab = "", xlab = "",
             main = "First random effect in the drift")
        abline(0, 1)
        plot(sort(x@estimphi[2, indexpred]), sort(phipred[2, ]/psipred), pch = 18,
             xlim = c(min(x@estimphi[2, ], phipred[2, ]) * 0.8, max(x@estimphi[2,
                                                                               ], phipred[2, ]) * 1.2), ylim = c(min(x@estimphi[2, ], phipred[2, ]) *
                                                                                                                   0.8, max(x@estimphi[2, ], phipred[2, ]) * 1.2), ylab = "", xlab = "",
             main = "Second random effect in the drift")
        abline(0, 1)
        
        plot(sort(sqrt(x@estimpsi2[indexpred])), sort(psipred), pch = 18, xlim = c(min(sqrt(x@estimpsi2[indexpred]),
                                                                                       psipred) * 0.8, max(sqrt(x@estimpsi2[indexpred]), psipred) * 1.2), ylim = c(min(sqrt(x@estimpsi2[indexpred]),
                                                                                                                                                                       psipred) * 0.8, max(sqrt(x@estimpsi2[indexpred]), psipred) * 1.2), ylab = "",
             xlab = "", main = "Random effect in the diffusion")
        abline(0, 1)
        
        plot(times, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue,
                                                                                          Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
        for (j in indexpred) {
          lines(times, Xtrue[j, ], col = j)
        }
        
        plot(timestrue, Xpred[1, ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue,
                                                                                   Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), main = "Predictive trajectories")
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
      
      return(list(phipred = phipred, psipred = psipred, Xpred = Xpred, indexpred = indexpred))
    }
  }
})

