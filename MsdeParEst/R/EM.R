# MsdeParEst R package ; file em.r (last modified: 2017-08-28)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017

EM <- function(U, V, S, K, drift.random, start, Niter = 10, drift.fixed = NULL,  
               sigma = NULL) {
    
    M <- dim(U)[2]
    
    mu <- start$mu
    omega <- start$omega
    mixt.prop <- start$mixt.prop
    
    N <- length(mixt.prop)
    
    estimphi <- matrix(NA, 2, M)
    
    for (j in 1:M) {
        estimphi[, j] <- solve(V[[j]]) %*% U[, j]
    }
    
    if (is.null(sigma)) {
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
            ln <- function(param){Q_EM(matrix(param[1:(2 * N)], nrow = N, ncol = 2), matrix(param[(2 * 
                N + 1):(4 * N)], nrow = N, ncol = 2), sigma, probindi, U, V, S, K, estimphi, 
                drift.random)}
            paraminit <- c(as.vector(mu), as.vector(omega))
            res <- optim(paraminit, fn = ln, method = "Nelder-Mead")
            muhat[iter, , ] <- matrix(res$par[1:(2 * N)], nrow = N, ncol = 2)
            omegahat[iter, , ] <- matrix(abs(res$par[(2 * N + 1):(4 * N)]), nrow = N, ncol = 2)
        }
        
        if (is.null(drift.fixed)) {

                if (sum(drift.random) == 2) {
                  ln <- function(param){Q_EM(matrix(param[1:(2 * N)], nrow = N, ncol = 2),
                    matrix(c(rep(0, N), param[(2 * N + 1):(3 * N)]), nrow = N, ncol = 2),
                    sigma, probindi, U, V, S, K, estimphi, drift.random)}
                  paraminit <- c(as.vector(mu), as.vector(omega[, 2]))
                  res <- optim(paraminit, fn = ln, method = "Nelder-Mead")
                  muhat[iter, , ] <- matrix(c(res$par[1:N], res$par[(N + 1):(2 * N)]),
                    nrow = N, ncol = 2)
                  omegahat[iter, , ] <- matrix(c(rep(0, N), abs(res$par[(2 * N + 1):(3 *
                    N)])), nrow = N, ncol = 2)

                }

                if (sum(drift.random) == 1) {
                  ln <- function(param){Q_EM(matrix(param[1:(2 * N)], nrow = N, ncol = 2),
                    matrix(c(param[(2 * N + 1):(3 * N)], rep(0, N)), nrow = N, ncol = 2),
                    sigma, probindi, U, V, S, K, estimphi, drift.random)}
                  paraminit <- c(as.vector(mu), as.vector(omega[, 1]))
                  res <- optim(paraminit, fn = ln, method = "Nelder-Mead")
                  muhat[iter, , ] <- matrix(res$par[1:(2 * N)], nrow = N, ncol = 2)
                  omegahat[iter, , ] <- matrix(c(abs(res$par[(2 * N + 1):(3 * N)]), rep(0,
                    N)), nrow = N, ncol = 2)
                }
            
        }
        
        if (!is.null(drift.fixed)) {
            if (sum(drift.random) == 2) {
                ln <- function(param){Q_EM(matrix(c(drift.fixed, param[1:N]), nrow = N,
                  ncol = 2), matrix(c(rep(0, N), param[(N + 1):(2 * N)]), nrow = N, ncol = 2),
                  sigma, probindi, U, V, S, K, estimphi, drift.random)}
                paraminit <- c(as.vector(mu[, 2]), as.vector(omega[, 2]))
                res <- optim(paraminit, fn = ln, method = "Nelder-Mead")
                muhat[iter, , ] <- matrix(c(drift.fixed, res$par[1:N]), nrow = N, ncol = 2)
                omegahat[iter, , ] <- matrix(c(rep(0, N), abs(res$par[(N + 1):(2 * N)])),
                  nrow = N, ncol = 2)

            }

            if (sum(drift.random) == 1) {
                ln <- function(param){Q_EM(matrix(c(param[1:N], drift.fixed), nrow = N,
                  ncol = 2), matrix(c(param[(N + 1):(2 * N)], rep(0, N)), nrow = N, ncol = 2),
                  sigma, probindi, U, V, S, K, estimphi, drift.random)}
                paraminit <- c(as.vector(mu[, 1]), as.vector(omega[, 1]))
                res <- optim(paraminit, fn = ln, method = "Nelder-Mead")
                muhat[iter, , ] <- matrix(c(res$par[1:N], drift.fixed), nrow = N, ncol = 2)
                omegahat[iter, , ] <- matrix(c(abs(res$par[(N + 1):(2 * N)]), rep(0, N)),
                  nrow = N, ncol = 2)
            }
        }
        
        
        
    }
    
    if (is.null(sigma)) {
        
        if (sum(drift.random) > 2) {
            nbparam <- 4 * N + 1
        }
        if (is.null(drift.fixed)) {
            if (sum(drift.random) == 2) {
                nbparam <- 3 * N + 1
            }
            if (sum(drift.random) == 1) {
                nbparam <- 3 * N + 1
            }
        }
        if (!is.null(drift.fixed)) {
            if (sum(drift.random) == 2) {
                nbparam <- 2 * N + 1
            }
            if (sum(drift.random) == 1) {
                nbparam <- 2 * N + 1
            }
        }
    }
    if (!is.null(sigma)) {
        
        if (sum(drift.random) > 2) {
            nbparam <- 4 * N
        }
        if (is.null(drift.fixed)) {
            if (sum(drift.random) == 2) {
                nbparam <- 3 * N
            }
            if (sum(drift.random) == 1) {
                nbparam <- 3 * N
            }
        }
        if (!is.null(drift.fixed)) {
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


