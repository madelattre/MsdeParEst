#' Estimation In Mixed Stochastic Differential Equations with fixed effects in the drift and one random effect in the diffusion coefficient
#' 
#' @description Parameter estimation of the mixed SDE with Gamma distribution of the diffusion random effect 
#' and fixed effects in the drift:
#' 
#'  \eqn{dXj(t)= (\alpha- \beta Xj(t))dt + \sigma_j a(Xj(t)) dWj(t), 1/\sigma_j^2 \sim \Gamma(a,lambda)}, 
#'  
#' done with \code{\link{likelihoodGamma}}.
#' @param U matrix of M sufficient statistics U (see \code{\link{UVS}}).
#' @param V list of the M sufficient statistics matrix V (see \code{\link{UVS}}).
#' @param S vector of the M sufficient statistics S (see \code{\link{UVS}}).
#' @param SigDelta vector of the M constant terms of the individual likelihood (see \code{\link{UVS}}).
#' @param K number of times of observations.
#' @param drift.param values of the fixed effects in the drift. Defaults to c(0,0).
#' @param drift.estim 1 if the fixed effects in the drift are estimated, 0 otherwise. Defaults to 1.
#' @return
#' \item{mu}{values of the fixed effects in the drift.}
#' \item{a}{estimated value of the shape of the Gamma distribution.}
#' \item{lambda}{estimated value of the scale of the Gamma distribution.}
#' \item{BIChere}{BIC indicator.}
#' \item{AIChere}{AIC indicator.}
#' @importFrom stats optim 
#' @importFrom stats var
#' 
#' @references
#' Estimation of population parameters in stochastic differential equations with random effects in the diffusion coefficient, M. Delattre, V. Genon-Catalot and A. Samson, \emph{ESAIM: Probability and Statistics 2015}, Vol 19, \bold{671 -- 688}


EstParamGamma <- function(U, V, S, SigDelta, K, drift.param = c(0, 0), drift.estim = 1) {
  
  M <- length(S)
  
  k <- 0.1
  
  estimGamma <- (K/S) * ((S/K) >= (k/sqrt(K)))
  init.a <- (mean(estimGamma))^2/var(estimGamma)
  init.lambda <- var(estimGamma)/mean(estimGamma)
  
  if (drift.estim == 0) {
    ln = function(param) {
      contrastGamma(exp(param[1]), exp(param[2]), U, V, S, K, drift.param)
    }
    
    res <- optim(c(log(init.a), log(init.lambda)), fn = ln, method = "Nelder-Mead")
    
    a <- exp(res$par[1])
    lambda <- exp(res$par[2])
    mu <- drift.param
    
    if ((a+K/2)<=171){
      BIChere <- likelihoodGamma(a, lambda, U, V, S, SigDelta, K, mu) + 2 * log(M)
      AIChere <- likelihoodGamma(a, lambda, U, V, S, SigDelta, K, mu) + 4
    }
    else{
      BIChere <- Inf
      AIChere <- Inf
    }
  }
  
  if (drift.estim == 1) {
    
    estimphi <- c(0, 0)
    
    U1 <- rep(0, M)
    U2 <- rep(0, M)
    V11 <- rep(0, M)
    V22 <- rep(0, M)
    V12 <- rep(0, M)
    detVU <- rep(0, M)
    estimphi <- matrix(NA, 2, M)
    
    for (j in 1:M) {
      U1[j] <- U[1, j]
      U2[j] <- U[2, j]
      V11[j] <- V[[j]][1, 1]
      V22[j] <- V[[j]][2, 2]
      V12[j] <- V[[j]][1, 2]
    }
    
    for (j in 1:M) {
      detVU[j] <- V11[j] * V22[j] - V12[j]^2
      estimphi[1, j] <- (V22[j] * U1[j] - V12[j] * U2[j])/detVU[j]
      estimphi[2, j] <- (V11[j] * U2[j] - V12[j] * U1[j])/detVU[j]
    }
    
    init.phi <- apply(estimphi, 1, mean)
    
    ln = function(param) {
      contrastGamma(exp(param[1]), exp(param[2]), U, V, S, K, c(param[3], param[4]))
    }
    
    
    res <- optim(c(log(init.a), log(init.lambda), init.phi), fn = ln, method = "Nelder-Mead")
    
    a <- exp(res$par[1])
    lambda <- exp(res$par[2])
    mu <- c(res$par[3], res$par[4])
    
    if ((a+K/2)<=171){
      BIChere <- likelihoodGamma(a, lambda, U, V, S, SigDelta, K, mu) + 4 * log(M)
      AIChere <- likelihoodGamma(a, lambda, U, V, S, SigDelta, K, mu) + 4 * 2
    }
    else{
      BIChere <- Inf
      AIChere <- Inf
    }
    
    
  }
  
  return(list(mu = mu, a = a, lambda = lambda, BIChere = BIChere, AIChere = AIChere))
}
