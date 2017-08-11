# MsdeParEst R package ; file EstParamGamma.r (last modified: 2017-08-11)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017
# UMR 518 AgroParisTech/INRA, 16 rue Claude Bernard, 75 231 Paris Cedex 05


EstParamGamma <- function(U, V, S, SigDelta, K, drift.param = NULL) {
  
  M <- length(S)
  
  k <- 0.1
  
  estimGamma <- (K/S) * ((S/K) >= (k/sqrt(K)))
  init.a <- (mean(estimGamma))^2/var(estimGamma)
  init.lambda <- var(estimGamma)/mean(estimGamma)
  
  if (!is.null(drift.param)) {
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
  
  if (is.null(drift.param)) {
    
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
