#' EM algorithm for mixtures of stochastic differential equations with random effects in the drift
#' 
#' @description EM algorithm for parameter estimation in the mixed SDE 
#'  \eqn{dXj(t)= (\alpha_j - \beta_j Xj(t))dt + \sigma a(Xj(t)) dWj(t)} with random effects
#'  in the drift following a mixture of Gaussian distributions.
#' @param U matrix of M sufficient statistics U (see \code{\link{UVS}}).
#' @param V list of the M sufficient statistics matrix V (see \code{\link{UVS}}).
#' @param S vector of the M sufficient statistics S (see \code{\link{UVS}}).
#' @param K number of times of observations.
#' @param drift.random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random 
#'  effect or c(1,2) if 2 random effects.
#' @param start list of starting values: mu, omega, mixt.prop. mixt.prop is a vector of length N, the number of 
#' mixture components. mu is a N x 2 matrix, first (resp. second) column is the mean of \eqn{\alpha_j} (resp. 
#' \eqn{\beta_j}) if \eqn{\alpha_j} (resp. \eqn{\beta_j}) is random, the fixed effect value otherwise. omega is a 
#' N x 2 matrix, the components corresponding to a fixed effect should be set to 0. 
#' @param Niter number of iterations. Defaults to 10.
#' @param drift.fixed value for the fixed effects in the drift if known (not estimated). Vector of length N. Defaults to 0.
#' @param drift.estim.fix 1 if the fixed effects in the drift are estimated, 0 otherwise. Defaults to 1.
#' @param drift.fixed.mixt 1 if the value of the fixed effect in the drift is different from one mixture component to
#' another, 0 otherwise. Defaults to 1. 
#' @param sigma value for the diffusion parameter if known (not estimated). Defaults to 1.
#' @param sigma.estim 1 if the diffusion parameter is estimated, 0 otherwise. Defaults to 1.  
#' @return
#' \item{mu}{estimated value of the mean at each iteration of the algorithm. Niter x N x 2 array. }
#' \item{omega}{estimated value of the standard deviation at each iteration of the algorithm. Niter x N x 2 array.}
#' \item{mixt.prop}{estimated value of the mixture proportions at each iteration of the algorithm. Niter x N matrix.}
#' \item{sigma}{value of the diffusion parameter.}
#' \item{probindi}{posterior component probabilites. M x N matrix.}
#' \item{BIChere}{BIC indicator}
#' \item{AIChere}{AIC indicator}
#' @references
#' Mixtures of stochastic differential equations with random effects: application to data clustering, M. Delattre, V. Genon-Catalot and A. Samson, \emph{Journal of Statistical Planning and Inference 2016}, Vol 173, \bold{109--124}



EM <- function(U, V, S, K, drift.random, start, Niter = 10, drift.fixed = 0, drift.estim.fixed = 1, 
    drift.fixed.mixt = 1, sigma = 1, sigma.estim = 1) {
    
    M <- dim(U)[2]
    
    mu <- start$mu
    omega <- start$omega
    mixt.prop <- start$mixt.prop
    
    N <- length(mixt.prop)
    
    estimphi <- matrix(NA, 2, M)
    
    for (j in 1:M) {
        estimphi[, j] <- solve(V[[j]]) %*% U[, j]
    }
    
    if (sigma.estim == 1) {
        sigma <- sqrt(mean(S/(K - 1)))
    }
    
    omegahat <- array(0, c(Niter, N, 2))
    omegahat[1, , ] <- omega/sigma
    
    muhat <- array(0, c(Niter, N, 2))
    muhat[1, , ] <- mu
    
    probhat <- matrix(0, nrow = Niter, ncol = N)
    probhat[1, ] <- mixt.prop
    
    for (iter in 2:Niter) {
        probindi <- probind(mu, omega, mixt.prop, sigma, U, V, S, K, estimphi, drift.random)
        mixt.prop <- colMeans(probindi, na.rm = TRUE)
        probhat[iter, ] <- mixt.prop
        mu <- muhat[iter - 1, , ]
        omega <- omegahat[iter - 1, , ]
        
        if (sum(drift.random) > 2) {
            ln <- function(param) (Q_EM(matrix(param[1:(2 * N)], nrow = N, ncol = 2), matrix(param[(2 * 
                N + 1):(4 * N)], nrow = N, ncol = 2), sigma, probindi, U, V, S, K, estimphi, 
                drift.random))
            paraminit <- c(as.vector(mu), as.vector(omega))
            res <- optim(paraminit, f = ln, method = "Nelder-Mead")
            muhat[iter, , ] <- matrix(res$par[1:(2 * N)], nrow = N, ncol = 2)
            omegahat[iter, , ] <- matrix(abs(res$par[(2 * N + 1):(4 * N)]), nrow = N, ncol = 2)
        }
        
        if (drift.estim.fixed == 1) {
            if (drift.fixed.mixt == 0) {
                if (sum(drift.random) == 2) {
                  ln <- function(param) (Q_EM(matrix(c(rep(param[1], N), param[2:(N + 1)]), 
                    nrow = N, ncol = 2), matrix(c(rep(0, N), param[(N + 2):(2 * N + 1)]), 
                    nrow = N, ncol = 2), sigma, probindi, U, V, S, K, estimphi, drift.random))
                  paraminit <- c(mu[1, 1], as.vector(mu[, 2]), as.vector(omega[, 2]))
                  res <- optim(paraminit, f = ln, method = "Nelder-Mead")
                  muhat[iter, , ] <- matrix(c(rep(res$par[1], N), res$par[2:(N + 1)]), 
                    nrow = N, ncol = 2)
                  omegahat[iter, , ] <- matrix(c(rep(0, N), abs(res$par[(N + 2):(2 * N + 
                    1)])), nrow = N, ncol = 2)
                  
                }
                
                if (sum(drift.random) == 1) {
                  ln <- function(param) (Q_EM(matrix(c(param[1:N], rep(param[N + 1], N)), 
                    nrow = N, ncol = 2), matrix(c(param[(N + 2):(2 * N + 1)], rep(0, N)), 
                    nrow = N, ncol = 2), sigma, probindi, U, V, S, K, estimphi, drift.random))
                  paraminit <- c(as.vector(mu[, 1]), mu[2, 1], as.vector(omega[, 1]))
                  res <- optim(paraminit, f = ln, method = "Nelder-Mead")
                  muhat[iter, , ] <- matrix(c(res$par[1:N], rep(res$par[N + 1], N)), nrow = N, 
                    ncol = 2)
                  omegahat[iter, , ] <- matrix(c(abs(res$par[(N + 2):(2 * N + 1)]), rep(0, 
                    N)), nrow = N, ncol = 2)
                }
            }
            if (drift.fixed.mixt == 1) {
                if (sum(drift.random) == 2) {
                  ln <- function(param) (Q_EM(matrix(param[1:(2 * N)], nrow = N, ncol = 2), 
                    matrix(c(rep(0, N), param[(2 * N + 1):(3 * N)]), nrow = N, ncol = 2), 
                    sigma, probindi, U, V, S, K, estimphi, drift.random))
                  paraminit <- c(as.vector(mu), as.vector(omega[, 2]))
                  res <- optim(paraminit, f = ln, method = "Nelder-Mead")
                  muhat[iter, , ] <- matrix(c(res$par[1:N], res$par[(N + 1):(2 * N)]), 
                    nrow = N, ncol = 2)
                  omegahat[iter, , ] <- matrix(c(rep(0, N), abs(res$par[(2 * N + 1):(3 * 
                    N)])), nrow = N, ncol = 2)
                  
                }
                
                if (sum(drift.random) == 1) {
                  ln <- function(param) (Q_EM(matrix(param[1:(2 * N)], nrow = N, ncol = 2), 
                    matrix(c(param[(2 * N + 1):(3 * N)], rep(0, N)), nrow = N, ncol = 2), 
                    sigma, probindi, U, V, S, K, estimphi, drift.random))
                  paraminit <- c(as.vector(mu), as.vector(omega[, 1]))
                  res <- optim(paraminit, f = ln, method = "Nelder-Mead")
                  muhat[iter, , ] <- matrix(res$par[1:(2 * N)], nrow = N, ncol = 2)
                  omegahat[iter, , ] <- matrix(c(abs(res$par[(2 * N + 1):(3 * N)]), rep(0, 
                    N)), nrow = N, ncol = 2)
                }
            }
        }
        
        if (drift.estim.fixed == 0) {
            if (sum(drift.random) == 2) {
                ln <- function(param) (Q_EM(matrix(c(drift.fixed, param[1:N]), nrow = N, 
                  ncol = 2), matrix(c(rep(0, N), param[(N + 1):(2 * N)]), nrow = N, ncol = 2), 
                  sigma, probindi, U, V, S, K, estimphi, drift.random))
                paraminit <- c(as.vector(mu[, 2]), as.vector(omega[, 2]))
                res <- optim(paraminit, f = ln, method = "Nelder-Mead")
                muhat[iter, , ] <- matrix(c(drift.fixed, res$par[1:N]), nrow = N, ncol = 2)
                omegahat[iter, , ] <- matrix(c(rep(0, N), abs(res$par[(N + 1):(2 * N)])), 
                  nrow = N, ncol = 2)
                
            }
            
            if (sum(drift.random) == 1) {
                ln <- function(param) (Q_EM(matrix(c(param[1:N], drift.fixed), nrow = N, 
                  ncol = 2), matrix(c(param[(N + 1):(2 * N)], rep(0, N)), nrow = N, ncol = 2), 
                  sigma, probindi, U, V, S, K, estimphi, drift.random))
                paraminit <- c(as.vector(mu[, 1]), as.vector(omega[, 1]))
                res <- optim(paraminit, f = ln, method = "Nelder-Mead")
                muhat[iter, , ] <- matrix(c(res$par[1:N], drift.fixed), nrow = N, ncol = 2)
                omegahat[iter, , ] <- matrix(c(abs(res$par[(N + 1):(2 * N)]), rep(0, N)), 
                  nrow = N, ncol = 2)
            }
        }
        
        
        
    }
    
    if (sigma.estim == 1) {
        
        if (sum(drift.random) > 2) {
            nbparam <- 4 * N + 1
        }
        if (drift.estim.fixed == 1) {
            if (sum(drift.random) == 2) {
                nbparam <- 3 * N + 1
            }
            if (sum(drift.random) == 1) {
                nbparam <- 3 * N + 1
            }
        }
        if (drift.estim.fixed == 0) {
            if (sum(drift.random) == 2) {
                nbparam <- 2 * N + 1
            }
            if (sum(drift.random) == 1) {
                nbparam <- 2 * N + 1
            }
        }
    }
    if (sigma.estim == 0) {
        
        if (sum(drift.random) > 2) {
            nbparam <- 4 * N
        }
        if (drift.estim.fixed == 1) {
            if (sum(drift.random) == 2) {
                nbparam <- 3 * N
            }
            if (sum(drift.random) == 1) {
                nbparam <- 3 * N
            }
        }
        if (drift.estim.fixed == 0) {
            if (sum(drift.random) == 2) {
                nbparam <- 2 * N
            }
            if (sum(drift.random) == 1) {
                nbparam <- 2 * N
            }
        }
    }
    
    BIChere <- likelihoodMixtureNormal(muhat[Niter, , ], omegahat[Niter, , ], sigma, mixt.prop, 
        U, V, S, K, estimphi, drift.random) + nbparam * log(M)
    AIChere <- likelihoodMixtureNormal(muhat[Niter, , ], omegahat[Niter, , ], sigma, mixt.prop, 
        U, V, S, K, estimphi, drift.random) + nbparam * 2
    
    
    output = list(mu = muhat, omega = omegahat * sigma, mixt.prop = probhat, sigma = sigma, 
        probindi = probindi, BIChere = BIChere, AIChere = AIChere)
}


#' Computation of the E-step of the EM algorithm for mixtures of stochastic differential equations with random effects
#' 
#' @description Computation of the E-step of the EM algorithm for parameter estimation in the mixed SDE 
#'  \eqn{dXj(t)= (\alpha_j - \beta_j Xj(t))dt + \sigma a(Xj(t)) dWj(t)} with random effects
#'  in the drift following a mixture of Gaussian distributions.
#' @param mu mean of the random effects. N x 2 matrix, first (resp. second) column is the mean of \eqn{\alpha_j} (resp. 
#' \eqn{\beta_j}) in each mixture component if \eqn{\alpha_j} (resp. \eqn{\beta_j}) is random, the fixed effect value otherwise. 
#' @param omega standard deviation of the random effects. N x 2 matrix, the components corresponding to a fixed effect should be 
#' set to 0. 
#' @param sigma value of the diffusion parameter.
#' @param probindi M x N matrix of individual component probabilites.
#' @param U matrix of M sufficient statistics U (see \code{\link{UVS}}).
#' @param V list of the M sufficient statistics matrix V (see \code{\link{UVS}}).
#' @param S vector of the M sufficient statistics S (see \code{\link{UVS}}).
#' @param K number of times of observations.
#' @param estimphi matrix of estimators of the fixed/random effects. 2 x M matrix. 
#' @param drift.random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random 
#'  effect or c(1,2) if 2 random effects.
#' @return
#' \item{Q}{value of the E-step.}
#' Mixtures of stochastic differential equations with random effects: application to data clustering, M. Delattre, V. Genon-Catalot and A. Samson, \emph{Journal of Statistical Planning and Inference 2016}, Vol 173, \bold{109--124}

Q_EM <- function(mu, omega, sigma, probindi, U, V, S, K, estimphi, drift.random) {
    
    M <- dim(probindi)[1]
    N <- dim(probindi)[2]
    
    Lindicomponent <- matrix(NA, M, N)
    
    for (i in 1:M) {
        Lindicomponent[i, ] <- sapply(1:N, function(n) {
            likelihoodNormalindi(mu[n, ], omega[n, ], sigma, U[, i], V[[i]], S[i], K, estimphi[, 
                i], drift.random)
        })
    }
    
    Q <- sum(probindi * Lindicomponent)
    
    return(Q = Q)
}

#' Computation of the component probabilities
#' 
#' @description Computation of the individual component probabilities in the mixed SDE 
#'  \eqn{dXj(t)= (\alpha_j - \beta_j Xj(t))dt + \sigma a(Xj(t)) dWj(t)} with random effects
#'  in the drift following a mixture of Gaussian distributions.
#' @param mu mean of the random effects. N x 2 matrix, first (resp. second) column is the mean of \eqn{\alpha_j} (resp. 
#' \eqn{\beta_j}) in each mixture component if \eqn{\alpha_j} (resp. \eqn{\beta_j}) is random, the fixed effect value otherwise. 
#' @param omega standard deviation of the random effects. N x 2 matrix, the components corresponding to a fixed effect should be 
#' set to 0. 
#' @param mixt.prop vector of mixture proportions.
#' @param sigma value of the diffusion parameter.
#' @param U matrix of M sufficient statistics U (see \code{\link{UVS}}).
#' @param V list of the M sufficient statistics matrix V (see \code{\link{UVS}}).
#' @param S vector of the M sufficient statistics S (see \code{\link{UVS}}).
#' @param K number of times of observations.
#' @param estimphi matrix of estimators of the fixed/random effects. 2 x M matrix. 
#' @param drift.random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random 
#'  effect or c(1,2) if 2 random effects.
#' @return
#' \item{probindi}{M x N matrix of individual component probabilities.}
#' Mixtures of stochastic differential equations with random effects: application to data clustering, M. Delattre, V. Genon-Catalot and A. Samson, \emph{Journal of Statistical Planning and Inference 2016}, Vol 173, \bold{109--124}


probind <- function(mu, omega, mixt.prop, sigma, U, V, S, K, estimphi, drift.random) {
    
    M <- dim(U)[2]
    N <- dim(mu)[1]
    
    probRep <- matrix(mixt.prop, M, N, byrow = TRUE)
    
    Lindicomponent <- matrix(NA, M, N)
    
    for (i in 1:M) {
        Lindicomponent[i, ] <- sapply(1:N, function(n) {
            exp(-1/2 * likelihoodNormalindi(mu[n, ], omega[n, ], sigma, U[, i], V[[i]], 
                S[i], K, estimphi[, i], drift.random))
        })
    }
    
    probindi <- exp(log(probRep) + log(Lindicomponent))
    
    probindi <- probindi/rowSums(probindi)
    
    return(probindi)
}

#' Computation of the individual Log Likelihood In Mixed Stochastic Differential Equations
#' 
#' @description Computation of -2 loglikelihood of individual j in the mixed SDE with Normal distribution of the random effects
#'  \eqn{dXj(t)= (\alpha_j- \beta_j Xj(t))dt + \sigma a(Xj(t)) dWj(t)}.
#' @param mu vector of mean of the random effects. First (resp. second) value is the mean of \eqn{\alpha_j} (resp. 
#' \eqn{\beta_j}) if \eqn{\alpha_j} (resp. \eqn{\beta_j}) is random, the fixed effect value otherwise. 
#' @param omega standard deviation of the random effects. The components corresponding to a fixed effect should be 
#' set to 0. 
#' @param sigma value of the diffusion parameter.
#' @param Uj vector of the sufficient statistics U for individual j (see \code{\link{UVS}}).
#' @param Vj matrix of the sufficient statistics V for individual j (see \code{\link{UVS}}).
#' @param Sj value of the sufficient statistic S for individual j (see \code{\link{UVS}}).
#' @param K number of times of observations.
#' @param estimphij vector  ofestimators of the random effects for individual j.
#' @param drift.random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random 
#'  effect or c(1,2) if 2 random effects.
#' @return
#' \item{L}{value of -2 x loglikelihood}
#' @references
#' Maximum likelihood estimation for stochastic differential equations with random effects, M. Delattre, V. Genon-Catalot and A. Samson, \emph{Scandinavian Journal of Statistics 2012}, Vol 40, \bold{322--343}


likelihoodNormalindi <- function(mu, omega, sigma, Uj, Vj, Sj, K, estimphij, drift.random) {
    
    if (length(drift.random) == 2) {
        Omega <- matrix(c(omega[1]^2, 0, 0, omega[2]^2), 2, 2, byrow = TRUE)
    }
    
    if (sum(drift.random) == 1) {
        Omega <- matrix(c(omega[1]^2, 0, 0, 0), 2, 2, byrow = TRUE)
    }
    
    if (sum(drift.random) == 2) {
        Omega <- matrix(c(0, 0, 0, omega[2]^2), 2, 2, byrow = TRUE)
    }
    
    
    I2 <- diag(c(1, 1))
    A <- (I2 + Vj %*% Omega)
    Rinv <- solve(A) %*% Vj
    b <- mu - estimphij
    
    L <- log(det(A)) + 1/sigma^2 * (t(b) %*% Rinv %*% b - t(Uj) %*% estimphij)
    
    return(L)
}

#' Computation of the Log Likelihood In Mixed Stochastic Differential Equations
#' 
#' @description Computation of -2 loglikelihood the mixed SDE 
#'  \eqn{dXj(t)= (\alpha_j - \beta_j Xj(t))dt + \sigma a(Xj(t)) dWj(t)} with random effects
#'  in the drift following a mixture of Gaussian distributions.
#' @param mu mean of the random effects. N x 2 matrix, first (resp. second) column is the mean of \eqn{\alpha_j} (resp. 
#' \eqn{\beta_j}) in each mixture component if \eqn{\alpha_j} (resp. \eqn{\beta_j}) is random, the fixed effect value otherwise. 
#' @param omega standard deviation of the random effects. N x 2 matrix, the components corresponding to a fixed effect should be 
#' set to 0. 
#' @param sigma value of the diffusion parameter.
#' @param mixt.prop vector of mixture proportions.
#' @param U matrix of M sufficient statistics U (see \code{\link{UVS}}).
#' @param V list of the M sufficient statistics matrix V (see \code{\link{UVS}}).
#' @param S vector of the M sufficient statistics S (see \code{\link{UVS}}).
#' @param K number of times of observations.
#' @param drift.random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random 
#'  effect or c(1,2) if 2 random effects.
#' @return
#' \item{L}{value of -2 x loglikelihood}
#' @references
#' Mixtures of stochastic differential equations with random effects: application to data clustering, M. Delattre, V. Genon-Catalot and A. Samson, \emph{Journal of Statistical Planning and Inference 2016}, Vol 173, \bold{109--124}


likelihoodMixtureNormal <- function(mu, omega, sigma, mixt.prop, U, V, S, K, estimphi, 
    drift.random) {
    
    M <- dim(U)[2]
    N <- length(mixt.prop)
    
    probRep <- matrix(mixt.prop, M, N, byrow = TRUE)
    
    Lindicomponent <- matrix(NA, M, N)
    
    for (i in 1:M) {
        Lindicomponent[i, ] <- sapply(1:N, function(n) {
            exp(-1/2 * likelihoodNormalindi(mu[n, ], omega[n, ], sigma, U[, i], V[[i]], 
                S[i], K, estimphi[, i], drift.random))
        })
    }
    
    probindi <- exp(log(probRep) + log(Lindicomponent))
    
    loglik <- -2 * sum(log(rowSums(probindi)))
    
    return(loglik)
}

