# MsdeParEst R package ; file msde.fit.r (last modified: 2017-08-28)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017

#' Estimation Of The Random Effects In Mixed Stochastic Differential Equations
#' 
#' 
#' @description Parametric estimation of the joint density of the random effects in the mixed SDE
#' 
#'  \deqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j \ a(X_j(t)) dW_j(t),}
#'  \eqn{j=1,\ldots,M}, where the \eqn{(W_j(t))} are independant Wiener processes and the \eqn{(X_j(t))} are observed without noise. 
#'  There can be random effects either in the drift \eqn{(\alpha_j,\beta_j)} or in the diffusion coefficient \eqn{\sigma_j} or both
#'  \eqn{(\alpha_j,\beta_j,\sigma_j)}.
#'  
#' @param times vector of observation times
#' @param X matrix of the M trajectories (each row is a trajectory with as much columns as observations)
#' @param model name of the SDE: 'OU' (Ornstein-Uhlenbeck) or 'CIR' (Cox-Ingersoll-Ross)
#' @param drift.random random effects in the drift: 0 if only fixed effects, 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects. Default to c(1,2)
#' @param drift.fixed NULL if the fixed effect(s) in the drift is (are) estimated, value of the fixed effect(s) otherwise. Default to NULL
#' @param diffusion.random 1 if \eqn{\sigma} is random, 0 otherwise. Default to 0
#' @param diffusion.fixed NULL if \eqn{\sigma} is estimated (if fixed), value of \eqn{\sigma} otherwise. Default to NULL
#' @param nb.mixt number of mixture components for the distribution of the random effects in the drift. Default to 1 (no mixture) 
#' @param Niter number of iterations for the EM algorithm if the random effects in the drift follow a mixture distribution. Default to 10
#' @param discrete 1 for using a contrast based on discrete observations, 0 otherwise. Default to 1 
#' @param valid 1 if test validation, 0 otherwise. Default to 0
#' @param level alpha for the predicion intervals. Default 0.05
#' @param newwindow logical(1), if TRUE, a new window is opened for the plot. Default to FALSE
#'
#' @return 
#' \item{index}{is the vector of subscript in 1,...,M used for the estimation. Most of the time index=1:M, except for the CIR that requires positive trajectories.}
#' \item{estimphi}{matrix of estimators of the drift random effects \eqn{\hat{\alpha}_j}, or \eqn{\hat{\beta}_j} or \eqn{(\hat{\alpha}_j,\hat{\beta}_j)}}
#' \item{estimpsi2}{vector of estimators of the squared diffusion random effects \eqn{\hat{\sigma}_j^2}}
#' \item{gridf}{grid of values for the plots of the random effects distribution in the drift, matrix form}
#' \item{gridg}{grid of values for the plots of the random effects distribution in the diffusion, matrix form}
#' \item{estimf}{estimator of the density of \eqn{\alpha_j}, \eqn{\beta_j} or \eqn{(\alpha_j,\beta_j)}. Matrix form.}
#' \item{estimg}{estimator of the density of \eqn{\sigma_j^2}. Matrix form.}
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
#' \item{estim.drift.fix}{1 if the fixed effects in the drift are estimated, 0 otherwise.}
#' \item{estim.diffusion.fixed}{1 if the fixed effect in the diffusion is estimated, 0 otherwise.}
#' \item{discrete}{initial choice}
#' \item{times}{initial choice}
#' \item{X}{initial choice}
#' 
#' For mixture distributions in the drift:
#' \item{mu}{estimated value of the mean at each iteration of the algorithm. Niter x N x 2 array. }
#' \item{omega}{estimated value of the standard deviation at each iteration of the algorithm. Niter x N x 2 array.}
#' \item{mixt.prop}{estimated value of the mixture proportions at each iteration of the algorithm. Niter x N matrix.}
#' \item{probindi}{posterior component probabilites. M x N matrix.}
#' 
#' @importFrom stats dnorm
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @importFrom stats quantile
#' @importFrom stats density
#' @importFrom stats optim
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom methods new
#' @importFrom methods slot
#' @importFrom methods slotNames
#' @importFrom stats kmeans
#' @importFrom stats sd
#' @importFrom graphics lines
#' 
#' @details
#' Parametric estimation of the random effects density from M independent trajectories of the SDE:
#' \deqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j \ a(X_j(t)) dW_j(t),}
#' \eqn{j=1,\ldots,M}, where the \eqn{(W_j(t))} are independant Wiener processes and the \eqn{(X_j(t))} are observed without noise. 
#' 
#' \bold{Specification of the random effects:}
#' 
#' The drift includes no, one or two random effects: 
#' \enumerate{
#' \item if drift.random = 0: \eqn{\alpha_j \equiv \alpha} and \eqn{\beta_j \equiv \beta} are fixed
#' \item if drift.random = 1: \eqn{\beta_j \equiv \beta} is fixed and \eqn{\alpha_j} is random
#' \item if drift.random = 2: \eqn{\alpha_j \equiv \alpha} is fixed and \eqn{\beta_j} is random  
#' \item if drift.random = c(1,2): \eqn{\alpha_j} and \eqn{\beta_j} are random
#' }
#' 
#' The diffusion includes either a fixed effect or a random effect:
#' \enumerate{
#' \item if diffusion.random = 0: \eqn{\sigma_j \equiv \sigma} is fixed
#' \item if diffusion.random = 1: \eqn{\sigma_j} is random
#' }
#' 
#' \bold{Distribution of the random effects}
#' 
#' If there is no random effect in the diffusion (diffusion.random = 0), there is at least on random effect in the drift that follows
#' \enumerate{
#' \item a Gaussian distribution (nb.mixt=1): 
#' \eqn{\alpha_j \sim N(\mu,\Omega)} or \eqn{\beta_j \sim N(\mu,\Omega)} or \eqn{(\alpha_j,\beta_j) \sim N(\mu,\Omega)},
#' \item or a mixture of Gaussian distributions (nb.mixt=K, K>1):
#' \eqn{\alpha_j \sim \sum_{k=1}^{K} p_k N(\mu_k,\Omega_k)} or \eqn{\beta_j \sim \sum_{k=1}^{K} p_k N(\mu_k,\Omega_k)} or \eqn{(\alpha_j,\beta_j) \sim \sum_{k=1}^{K} p_k N(\mu_k,\Omega_k)},
#' where \eqn{\sum_{k=1}^{K} p_k=1.}
#' } 
#' 
#' If there is one random effect in the diffusion (diffusion.random = 1), \eqn{1/\sigma_j^2 \sim \Gamma(a,\lambda)}, and the coefficients
#' in the drift are conditionally Gaussian: \eqn{\alpha_j|\sigma_j \sim N(\mu,\sigma_j^2 \Omega)} or \eqn{\beta_j|\sigma_j \sim N(\mu,\sigma_j^2 \Omega)} 
#' or \eqn{(\alpha_j,\beta_j)|\sigma_j \sim N(\mu,\sigma_j^2 \Omega)}, or they are fixed \eqn{\alpha_j \equiv \alpha, \beta_j \equiv \beta}.
#' 
#' \bold{SDEs}
#' 
#' Two diffusions are implemented: 
#' \enumerate{
#' \item the Ornstein-Uhlenbeck model (OU) \eqn{a(X_j(t))=1}
#' \item the Cox-Ingersoll-Ross model (CIR) \eqn{a(X_j(t))=\sqrt{X_j(t)}}
#' }
#' 
#' 
#' \bold{Estimation}
#' 
#' \itemize{
#' \item If discrete = 0, the estimation is based on the exact likelihood associated with continuous observations ([1],[3]). This is only possible if diffusion.random = 0 and
#' \eqn{\sigma} is not estimated by maximum likelihood but empirically by means of the quadratic variations. 
#' \item If discrete = 1, the likelihood of the Euler scheme of the mixed SDE is computed and maximized for estimating all the parameters.
#' \item If nb.mixt > 1, an EM algorithm is implemented and the number of iterations of the algorithm must be specified with Niter.
#' \item If valid = 1, two-thirds of the sample trajectories are used for estimation, while the rest is used for validation. A plot is then provided for
#' visual comparison between the true trajectories of the test sample and some predicted trajectories simulated under the estimated model.
#' }
#'  
#'  
#'  
#' @examples
#'
#' 
#' # Example 1 : One random effect in the drift and one random effect in the diffusion
#'
#' sim <- msde.sim(M = 25, T = 1, N = 1000, model = 'OU', 
#'                 drift.random = 2, drift.param = c(0,0.5,0.5), 
#'                 diffusion.random = 1, diffusion.param = c(8,1/2))
#'
#' res <- msde.fit(times = sim$times, X = sim$X, model = 'OU', drift.random = 2, 
#' diffusion.random = 1)
#' 
#' summary(res)
#' plot(res)
#' 
#' \dontrun{
#'
#' # Example 2 : one mixture of two random effects in the drift, and one fixed effect in
#' # the diffusion coefficient
#'
#' sim <- msde.sim(M = 100, T = 5, N = 5000, model = 'OU', drift.random = c(1,2),
#'                 diffusion.random = 0, 
#'                 drift.param = matrix(c(0.5,1.8,0.25,0.25,1,2,0.25,0.25),nrow=2,byrow=FALSE), 
#'                 diffusion.param = 0.1, nb.mixt = 2, mixt.prop = c(0.5,0.5))
#'
#' # -- Estimation without validation
#' res <- msde.fit(times = sim$times, X = sim$X, model = 'OU', drift.random = c(1,2),
#'                 nb.mixt=2, Niter = 10)
#'
#' summary(res)
#' plot(res)
#' 
#' # -- Estimation with prediction
#' res.valid <- msde.fit(times = sim$times, X = sim$X, model = 'OU', drift.random = c(1,2),
#'                       nb.mixt=2, Niter = 10, valid = 1)
#'
#' summary(res.valid)
#' plot(res.valid)
#' 
#' # Example 3 : CIR with one random effect in the drift and one random effect in the diffusion 
#' # coefficient
#'
#' sim <- msde.sim(M = 100, T = 5, N = 5000, model = 'CIR', drift.random = 2,
#'                 diffusion.random = 1, drift.param = c(4,1,0.1), diffusion.param = c(8,0.5),
#'                 X0 = 1)
#'
#' res <- msde.fit(times = sim$times, X = sim$X, model = 'CIR', drift.random = 2,
#'                 diffusion.random = 1)
#'
#' summary(res)
#'   }
#' 
#' @references See  
#' 
#' \bold{[1]} Maximum Likelihood Estimation for Stochastic Differential Equations with Random Effects, Delattre, M., Genon-Catalot, V. and Samson, A. \emph{Scandinavian Journal of Statistics 40(2) 2012} \bold{322-343} 
#' 
#' \bold{[2]} Estimation of population parameters in stochastic differential equations with random effects in the diffusion coefficient, Delattre, M., Genon-Catalot, V. and Samson, A. \emph{ESAIM:PS 19 2015} \bold{671-688}
#' 
#' \bold{[3]} Mixtures of stochastic differential equations with random effects: application to data clustering, Delattre, M., Genon-Catalot, V. and Samson, A. \emph{Journal of Statistical Planning and Inference 173 2016} \bold{109-124} 
#' 
#' \bold{[4]} Parametric inference for discrete observations of diffusion processes with mixed effects, Delattre, M., Genon-Catalot, V. and Laredo, C. \emph{hal-01332630 2016}
#' 
#' \bold{[5]} Estimation of the joint distribution of random effects for a discretely observed diffusion with random effects, Delattre, M., Genon-Catalot, V. and Laredo, C. \emph{hal-01446063 2017}
#' 
#' @author Maud Delattre and Charlotte Dion
#' @export

msde.fit <- function(times, X, model = c("OU", "CIR"), drift.random = c(1,2), drift.fixed = NULL, 
                     diffusion.random = 0, diffusion.fixed = NULL, nb.mixt = 1,  
                     Niter = 10, discrete = 1, valid = 0, level = 0.05, newwindow = FALSE) {
  
  #model <- match.arg(model)
  
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
  
  if (((nb.mixt - round(nb.mixt)) != 0) || (nb.mixt <= 0)){
    stop("The number of mixture components (nb.mixt) should be a positive integer")
  } else 
  {
    if ((nb.mixt > 1) && (((Niter - round(Niter)) != 0) || (Niter <= 0))){
      stop("Invalid value for Niter. The number of iterations of the EM algorithm should be a positive integer")
    }
    
  }

  if (!is.null(diffusion.fixed) && (diffusion.fixed <= 0)){
    stop("Invalid number for diffusion.fixed. The diffusion coefficient should be positive")
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
  
 
  
  if ((valid == 1) && ((level <= 0) || (level >=1))){
    stop("The value of alpha for the prediction intervals should belong to ]0,1[")
  } 
  
  if ((valid != 0) && (valid != 1)){
    stop("Invalid value for valid: should either 0 or 1")
  }
  
  if ((discrete != 0) && (discrete != 1)){
    stop("Invalid value for discrete: should either 0 or 1")
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
  
  aic <- 0
  bic <- 0
  
  
  ## Preparation of sample
  
  # -- Valid trajectories
  
  if (model == "OU") {
    Mindex <- M
    index <- 1:Mindex
    Xindex <- X[index, ] 
  }
  
  if (model == "CIR") {
    
    index <- which(rowSums(X <= 0) == 0)
    Mindex <- length(index)
    Xindex <- X[index, ] 
    if (Mindex == 0) {
      warning("All the trajectories have non positive values the model CIR cannot be used", call. = FALSE)
    }
    
  }
  
  # -- Learning sample vs test sample
  
  
  delta <- round(diff(times), 10)[1] 
  Tend <- times[length(times)] 
  
  if ((Mindex < 100) && (valid == 1)) {
    valid <- 0
    warning("The sample is too small. No test sample is taken.")
  }
  
  
  if (valid == 1){
    Mpred <- floor(Mindex * 2/3) + 1
    ipred <- sample(seq(1,Mindex,1),Mpred,replace=F)
    Xtrue <- Xindex[-ipred,]
    Xestim <- Xindex[ipred,]
    indexestim <- index[ipred]
    
    if ((Mpred < 30) && (valid == 1)) {
      valid <- 0
      warning("The sample is too small. No test sample is taken.")
    }
    
    Mindex <- Mpred
  }
  
  if (valid == 0){
    Xestim <- Xindex
    Xtrue <- as.matrix(0)
    indexestim <- index
  }
  
  
  
  
  
  ## Estimation
  
  if (nb.mixt == 1) {
    
    # -- Errors and warnings
    
    if (((sum(drift.random)) == 0) && (diffusion.random == 0)) {
      stop("There should be at least one random effect either in the drift or the diffusion coefficient.")
    }
    
    
    # -- Computation of the sufficient statistics
    
    U <- matrix(0, 2, Mindex)
    V <- as.list(1:Mindex)
    S <- rep(0, Mindex)
    SigDelta <- rep(0, Mindex)
    
    estimUV <- UVS(Xestim, model, times)
    
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
        
        gridf <- matrix(0, 2, 1000)
        gridf[1, ] <- seq(mu[1] - 3 * omega[1], mu[1] + 3 * omega[1], length = 1000)
        gridf[2, ] <- seq(mu[2] - 3 * omega[2], mu[2] + 3 * omega[2], length = 1000)
        
        
        estimf1 <- dnorm(gridf[1, ], mean = mu[1], sd = abs(omega[1]))
        estimf2 <- dnorm(gridf[2, ], mean = mu[2], sd = abs(omega[2]))
        
        estimf <- rbind(estimf1, estimf2)
        
      }
      
      if (sum(drift.random) == 2) {
        
        estimphi <- estimphi[2, ]
        
        
        gridf <- seq(mu[2] - 3 * omega[2], mu[2] + 3 * omega[2], length = 1000)
        
        
        estimf <- matrix(dnorm(gridf, mean = mu[2], sd = omega[2]), 1, length(gridf), byrow = TRUE)
        gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
        estimphi <- matrix(estimphi, 1, length(estimphi), byrow = TRUE)
        
      }
      
      if (sum(drift.random) == 1) {
        
        estimphi <- estimphi[1, ]
        
        gridf <- seq(mu[1] - 3 * omega[1], mu[1] + 3 * omega[1], length = 1000)
        
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
      
      gridg <- seq(min(estimpsi2) * 0.8, max(estimpsi2) * 1.2, length = 1000)
      
      simupsi2 <- 1/rgamma(1000, a, rate = 1/lambda)
      
      testpsi <- suppressMessages(density(simupsi2, from = min(simupsi2), to = max(simupsi2), bw = "ucv", 
                                          n = 1000))
      if (testpsi$bw < 0.1) {
        testpsi <- suppressMessages(density(simupsi2, from = min(simupsi2), to = max(simupsi2), n = 1000))
      }
      
      gridg <- seq(min(estimpsi2) * 0.8, max(estimpsi2) * 1.2, length = 1000)
      estimg <- matrix(testpsi$y, nrow = 1)
      gridg <- matrix(gridg, nrow = 1)
      
      estimpsi2 <- matrix(estimpsi2, 1, length(estimpsi2), byrow = TRUE)
    
      estimphi <- as.matrix(0)  
      
    }
    
    if (((sum(drift.random)) >= 1) && (diffusion.random == 1)) {
      
      res <- EstParamNormalGamma(U, V, S, SigDelta, K, drift.random, drift.fixed)
      bic <- res$BIChere
      aic <- res$AIChere
      mu <- res$mu
      omega <- res$omega
      a <- res$a
      lambda <- res$lambda
      
      # -- Estimator of the density
      
      estimpsi2 <- matrix(estimpsi2, 1, length(estimpsi2), byrow = TRUE)
      
      simupsi2 <- 1/rgamma(1000, a, rate = 1/lambda)
      testpsi <- suppressMessages(density(simupsi2, from = min(simupsi2), to = max(simupsi2), bw = "ucv", 
                                          n = 1000))
      if (testpsi$bw < 0.1) {
        testpsi <- suppressMessages(density(simupsi2, from = min(simupsi2), to = max(simupsi2), n = 1000))
      }
      
      gridg <- seq(min(estimpsi2) * 0.8, max(estimpsi2) * 1.2, length = 1000)
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
                                                                                                   ]), bw = "ucv", n = 1000))
        test2 <- suppressMessages(density(simuphi[2, ], from = min(simuphi[2, ]), to = max(simuphi[2, 
                                                                                                   ]), bw = "ucv", n = 1000))
        
        if (test1$bw < 0.1) {
          test1 <- suppressMessages(density(simuphi[1, ], from = min(simuphi[1, ]), to = max(simuphi[1, 
                                                                                                     ]), n = 1000))
        }
        estimf1 <- test1$y
        gridf1 <- test1$x
        
        if (test2$bw < 0.1) {
          test2 <- suppressMessages(density(simuphi[2, ], from = min(simuphi[2, ]), to = max(simuphi[2, 
                                                                                                     ]), n = 1000))
        }
        estimf2 <- test2$y
        gridf2 <- test2$x
        
        gridf <- matrix(0, 2, 1000)
        gridf[1, ] <- gridf1
        gridf[2, ] <- gridf2
        
        estimf <- matrix(0, 2, 1000)
        estimf[1, ] <- estimf1
        estimf[2, ] <- estimf2
        
      }
      
      if (sum(drift.random) == 2) {
        
        estimphi <- estimphi[2,]
        estimphi <- matrix(estimphi, 1, length(estimphi), byrow = TRUE)
        
        simuphi <- rep(NA, length(simupsi2))
        simuphi <- sapply(1:length(simupsi2), function(s) {
          rnorm(1, mean = mu[2], sd = abs(omega[2] * sqrt(simupsi2[s])))
        })
        
        
        test <- suppressMessages(density(simuphi, from = min(simuphi), to = max(simuphi), bw = "ucv", 
                                         n = 1000))
        if (test$bw < 0.1) {
          test <- suppressMessages(density(simuphi, from = min(simuphi), to = max(simuphi), n = 1000))
        }
        
        gridf <- test$x
        estimf <- test$y
        
        estimf <- matrix(estimf, nrow = 1)
        gridf <- matrix(gridf, nrow = 1)
        
      }
      
      if (sum(drift.random) == 1) {
        
        estimphi <- estimphi[1,]
        estimphi <- matrix(estimphi, 1, length(estimphi), byrow = TRUE)
        
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
    
    
    return(new(Class = "Fit.class", model = model, drift.random = drift.random, diffusion.random = diffusion.random, 
               gridf = gridf, mu = mu, omega = omega, a = a, lambda = lambda, sigma2 = sigma2, index = index, indexestim = indexestim,
               estimphi = estimphi, estimpsi2 = estimpsi2, estimf = estimf, estim.drift.fix = as.numeric(is.null(drift.fixed)), 
               estim.diffusion.fix = as.numeric(is.null(diffusion.fixed)), 
               discrete = discrete, bic = bic, aic = aic, times = times, X = X, gridg = gridg, estimg = estimg))
  }
  
  
  if (nb.mixt > 1) {
    
    gridf <- NULL
    estimf <- NULL
    
    
    
    if (diffusion.random == 1) {
      warning("For the considered mixtures of SDE, there should not be random effects in the diffusion coefficient. diffusion.random is set to 0.")
      diffusion.random <- 0
    }
    
    if (sum(drift.random) == 0) {
      stop("There should be at least one random effect in the drift.")
    }
    
    # -- Computation of the sufficient statistics
    
    U <- matrix(0, 2, Mindex)
    V <- as.list(1:Mindex)
    S <- rep(0, Mindex)
    SigDelta <- rep(0, Mindex)
    
    estimUV <- UVS(Xestim, model, times)
    
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
      
      gridf <- matrix(0, 2, 1000)
      bg1 <- max(abs(estimphi[1, ])) * 1.2
      bg2 <- max(abs(estimphi[2, ])) * 1.2
      gridf[1, ] <- seq(min(abs(estimphi[1, ])) * 0.8, max(abs(estimphi[1, ])) * 1.2, length = 1000)
      gridf[2, ] <- seq(min(abs(estimphi[2, ])) * 0.8, max(abs(estimphi[2, ])) * 1.2, length = 1000)
      
      
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
      
      gridf <- matrix(0, 2, 1000)
      gridf <- seq(min(abs(estimphi[1, ])) * 0.8, max(abs(estimphi[1, ])) * 1.2, length = 1000)
      
      
      estimf <- rep(0, length(gridf), byrow = TRUE)
      
      for (n in 1:nb.mixt) {
        estimf <- estimf + mixt.prop[n] * dnorm(gridf, mean = mu[Niter, n, 1], sd = abs(omega[Niter, n, 
                                                                                              1]))
      }
      
      gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
      estimf <- matrix(estimf, 1, length(estimf), byrow = TRUE)
      estimphi <- estimphi[1,]
      estimphi <- matrix(estimphi, 1, length(estimphi), byrow = TRUE)
    }
    
    if (sum(drift.random) == 2) {
      
      gridf <- matrix(0, 2, 1000)
      gridf <- seq(min(abs(estimphi[2, ])) * 0.8, max(abs(estimphi[2, ])) * 1.2, length = 1000)
      
      
      estimf <- rep(0, length(gridf), byrow = TRUE)
      
      for (n in 1:nb.mixt) {
        estimf <- estimf + mixt.prop[n] * dnorm(gridf, mean = mu[Niter, n, 2], sd = abs(omega[Niter, n, 
                                                                                              2]))
      }
      
      gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
      estimf <- matrix(estimf, 1, length(estimf), byrow = TRUE)
      estimphi <- estimphi[2,]
      estimphi <- matrix(estimphi, 1, length(estimphi), byrow = TRUE)
    }
    
    return(new(Class = "Mixture.fit.class", model = model, drift.random = drift.random, gridf = gridf, mu = mu, 
               omega = omega, mixt.prop = mixt.prop, sigma2 = sigma2, index = index, indexestim = indexestim, probindi = probindi, 
               estimf = estimf, estimphi = estimphi, bic = bic, aic = aic, estim.drift.fix = as.numeric(is.null(drift.fixed)), times = times, 
               X = X))
    
    
  }
  
  ## Validation plots, if a test sample is taken 
  
  if (valid == 1){
    
    
    
    if (nb.mixt > 1) {
      
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
        
        #if (plot.pred == TRUE) {
          op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                    cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
          
          plot(times, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 
                                                                                              0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
          for (j in indexpred) {
            lines(times, Xtrue[j, ], col = j)
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
        #}
        
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
        
        #if (plot.pred == TRUE) {
          
          op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                    cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
          
          plot(times, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 
                                                                                          0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
          for (j in indexpred) {
            lines(times, Xtrue[j, ], col = j)
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
        #}
      }
      
    }
    
    if (nb.mixt == 1) {
      
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
          
          #if (plot.pred == TRUE) {
            op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                      cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
            
            plot(times, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, 
                                                                                                  Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
            for (j in indexpred) {
              lines(times, Xtrue[j, ], col = j)
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
          #}
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
          
          #if (plot.pred == TRUE) {
            
            op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                      cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
            
            
            plot(times, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 
                                                                                            0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
            for (j in indexpred) {
              lines(times, Xtrue[j, ], col = j)
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
          #}
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
          
          #if (plot.pred == TRUE) {
            op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                      cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
            
            
            plot(times, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, 
                                                                                                  Xpred) * 0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
            for (j in indexpred) {
              lines(times, Xtrue[j, ], col = j)
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
          #}
        
        
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
          
          #if (plot.pred == TRUE) {
            
            op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.8, 1.8), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0), 
                      cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
            
            
            plot(times, Xtrue[indexpred[1], ], type = "l", xlab = "", ylab = "", ylim = c(min(Xtrue, Xpred) * 
                                                                                            0.8, max(Xtrue, Xpred) * 1.5), main = "True trajectories")
            for (j in indexpred) {
              lines(times, Xtrue[j, ], col = j)
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
          #}
        }
        
        
      }
    }
  }
  
  }
  
}


#' S4 class for the estimation results in the mixed SDE with random effects in the drift, in the diffusion or both 
#'  
#' @slot model 'OU' or 'CIR' (character)
#' @slot drift.random 0, 1, 2, or c(1,2) (numeric)
#' @slot diffusion.random 0 or 1 (numeric)
#' @slot gridf matrix of values on which the estimation of the density of the random effects in the drift is done (matrix)
#' @slot gridg matrix of values on which the estimation of the density of the random effects in the diffusion is done (matrix)
#' @slot mu estimator of the mean mu of the drift random effects (numeric)
#' @slot omega estimator of the variance of the drift random effects (numeric)
#' @slot a estimator of the shape of the Gamma distribution for the diffusion random effect (numeric)
#' @slot lambda estimator of the scale of the Gamma distribution for the diffusion random effect (numeric)
#' @slot sigma2 estimated value of \eqn{\sigma^2} if the diffusion coefficient is not random (numeric)
#' @slot index index of the valid trajectories for the considered model (numeric)
#' @slot indexestim index of the trajectories used for the estimation (numeric)
#' @slot estimphi matrix of the estimator of the drift random effects (matrix)
#' @slot estimpsi2 vector of the estimator of the diffusion random effects \eqn{\sigma_j^2} (numeric)
#' @slot estimf estimator of the (conditional) density of the drift random effects (numeric)
#' @slot estimg estimator of the density of \eqn{\sigma_j^2} (numeric)
#' @slot estim.drift.fix 1 if the user asked for the estimation of fixed parameter in the drift (numeric)
#' @slot estim.diffusion.fix 1 if the user asked for the estimation of fixed diffusion coefficient (numeric)
#' @slot discrete 1 if the estimation is based on the likelihood of discrete observations, 0 otherwise (numeric)
#' @slot bic bic (numeric)
#' @slot aic aic (numeric)
#' @slot times vector of observation times, storage of input variable (numeric)
#' @slot X matrix of observations, storage of input variable (matrix)


setClass(Class = "Fit.class", representation = representation(model = "character", drift.random = "numeric", diffusion.random = "numeric", 
                                                              gridf = "matrix", gridg = "matrix", mu = "numeric", omega = "numeric", a = "numeric", lambda = "numeric", 
                                                              sigma2 = "numeric", index = "numeric", indexestim = "numeric", estimphi = "matrix", estimpsi2 = "matrix", estimf = "matrix", 
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
#' @slot index index of the valid trajectories for the considered model (numeric)
#' @slot indexestim index of the trajectories used for the estimation (numeric)
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
                                                                      index = "numeric", indexestim = "numeric", estimphi = "matrix", probindi = "matrix", estimf = "matrix", estim.drift.fix = "numeric", 
                                                                      bic = "numeric", aic = "numeric", times = "numeric", X = "matrix"))


########################################################### OUTPUTS

out <- function(x) {
  sN <- slotNames(x)
  res <- lapply(sN, function(name) slot(x, name))
  names(res) <- sN
  res
}



# Summary for Fit.class -----------------------------------------------------------------



############################################################### SUMMARY

#' Short summary of the results of class object Fit.class
#' @exportMethod summary
#' @description Method for the S4 class Fit.class
#' @param object Fit.class class
#'

setMethod(f = "summary", signature = "Fit.class", definition = function(object) {
  
  cat(c("Number of trajectories used for estimation: ",length(object@indexestim),"\n"))
  
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
#' @exportMethod summary
#' @description Method for the S4 class Mixture.fit.class
#' @param object Mixture.fit.class class

setMethod("summary", "Mixture.fit.class", function(object) {
  
  cat(c("Number of trajectories used for estimation: ",length(object@indexestim),"\n"))
  
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
    
    cat("\nFixed effect in the diffusion coefficient:\n")
    
      print(matrix(c("sigma", round(sqrt(object@sigma2), 6)), 1, 2, byrow = TRUE), quote = FALSE, 
            right = TRUE)
    
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
#' @exportMethod plot
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
        plot(x@gridg, x@estimg, main = expression(sigma[j]^2), xlab = "Value", ylab = "Density", type = "l")
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
        plot(x@gridg, x@estimg, main = expression(sigma[j]^2), xlab = "Value", ylab = "Density", type = "l")
        title("Estimated densities", outer = TRUE)
      }
    }
    
    if (dim(x@gridf)[1] == 2) {
      
      op <- par(mfrow = c(2, 2), mar = c(2.8, 2.8, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(1, 1, 1, 1), omi = c(0.2, 
                                                                                                             0.2, 0.2, 0.2), cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
      plot(x@gridf[1, ], x@estimf[1, ], main = expression(alpha[j]), xlab = "Value", ylab = "Density", type = "l")
      plot(x@gridf[2, ], x@estimf[2, ], main = expression(beta[j]), xlab = "Value", ylab = "Density", type = "l")
      plot(x@gridg, x@estimg, main = expression(sigma[j]^2), xlab = "Value", ylab = "Density", type = "l")
      title("Estimated densities", outer = TRUE)
      
    }
    
  }
  
})


########################################################### PLOT
#' Plot method for the mixture estimation class object
#'
#' @exportMethod plot
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




