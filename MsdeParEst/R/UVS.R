#' Computation Of The Sufficient Statistics 
#' 
#' @description Computation of the sufficient statistics of the (approximate) likelihood of the mixed SDE
#'  \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j a(X_j(t)) dW_j(t)}.
#' @param X matrix of the M trajectories.
#' @param model name of the SDE: 'OU' (Ornstein-Uhlenbeck) or 'CIR' (Cox-Ingersoll-Ross).
#' @param times times vector of observation times.
#' @return
#' \item{U}{vector of the M statistics U(Tend)}
#' \item{V}{list of the M matrices V(Tend)}
#' \item{S}{vector of the M quadratic variations S(X(t_1),...,X(t_n))}
#' \item{SigDelta}{vector of the M constant contributions to the Euler scheme approximation to the likelihood SigDelta(X(t_1),...,X(t_n))}
#' @references See  
#' Maximum Likelihood Estimation for Stochastic Differential Equations with Random Effects, Delattre, M., Genon-Catalot, V. and Samson, A. \emph{Scandinavian Journal of Statistics 40(2) 2012} \bold{322-343} 
#' Estimation of population parameters in stochastic differential equations with random effects in the diffusion coefficient, Delattre, M., Genon-Catalot, V. and Samson, A. \emph{ESAIM:PS 19 2015} \bold{671-688}
#' Mixtures of stochastic differential equations with random effects: application to data clustering, Delattre, M., Genon-Catalot, V. and Samson, A. \emph{Journal of Statistical Planning and Inference 173 2016} \bold{109-124} 
#' Parametric inference for discrete observations of diffusion processes with mixed effects, Delattre, M., Genon-Catalot, V. and Larédo, C. \emph{hal-01332630 2016}
#' Estimation of the joint distribution of random effects for a discretely observed diffusion with random effects, Delattre, M., Genon-Catalot, V. and Larédo, C. \emph{hal-01446063 2017}
#' @keywords sufficient statistics, quadratic variations
#' @details 
#' Computation of the sufficient statistics of the (approximate) likelihood of the mixed SDE
#' \eqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j a(X_j(t)) dW_j(t) = (\alpha_j, \beta_j)b(X_j(t))dt + \sigma_j a(X_j(t)) dW_j(t)} with \eqn{b(x)=(1,-x)^t}:
#' 
#' U : \eqn{U(Tend) = \int_0^{Tend} b(X(s))/a^2(X(s))dX(s) }
#'  
#' V : \eqn{V(Tend) = \int_0^{Tend} b(X(s))^2/a^2(X(s))ds }
#' 
#' S : \eqn{S(X(t_1),...,X(t_n)) = 1/delta \sum_{j=1}^{n} (X(t_j)-X(t_{j-1}))^2/a^2(X(t_{j-1}))}
#'  
#' SigDelta: \eqn{SigDelta(X(t_1),...,X(t_n)) = n log(delta) + \sum_{j=1}^{n} log(a(X(t_j)))}

UVS <- function(X, model, times) {
    
    M <- dim(X)[1]
    K <- dim(X)[2]
    Xm <- X[, -K]
    delta <- diff(times)
    Tend <- times[K]
    
    S <- rep(NA, M)
    
    if (model == "OU") {
        for (j in 1:M) {
            S[j] <- sum((X[j, 2:K] - X[j, 1:(K - 1)])^2/delta)
        }
    }
    
    if (model == "CIR") {
        for (j in 1:M) {
            S[j] <- sum((X[j, 2:K] - X[j, 1:(K - 1)])^2/(delta * X[j, 1:(K - 1)]))
        }
    }
    
    
    U <- matrix(0, 2, M)
    
    V <- as.list(1:M)
    b <- as.list(1:M)
    
    Int1 <- rowSums(Xm * matrix(delta, M, length(delta)))  #Int1 <- apply(Xm * delta, 1, sum)
    
    if (model == "OU") {
        
        SigDelta <- rep(sum(log(delta)), M)
        
        Int2 <- rowSums(Xm^2 * matrix(delta, M, length(delta)))
        
        
        for (j in 1:M) {
            b[[j]] <- matrix(apply(matrix(X[j, ], 1, K), 2, bx, 0, c(1, 2)), 2, K)  # 2xK  matrix
            
            U[, j] <- rowSums((b[[j]][, 1:(K - 1)] * matrix((X[j, 2:K] - X[j, 1:(K - 1)]), 
                2, K - 1, byrow = TRUE)))
            
            V[[j]] <- matrix(c(Tend, -Int1[j], -Int1[j], Int2[j]), 2, 2)
            
        }
        
    }
    
    if (model == "CIR") {
        
        SigDelta <- rep(0, M)
        
        Int3 <- rowSums(1/Xm * matrix(delta, M, length(delta)))
        
        b <- as.list(1:M)
        bsig <- as.list(1:M)
        
        for (j in 1:M) {
            b[[j]] <- matrix(apply(matrix(X[j, ], 1, K), 2, bx, 0, c(1, 2)), 2, K)
            
            bsig[[j]] <- matrix(-1, 2, K)
            bsig[[j]][1, ] <- 1/X[j, ]
            
            U[, j] = rowSums((bsig[[j]][, 1:(K - 1)] * matrix((X[j, 2:K] - X[j, 1:(K - 
                1)]), 2, K - 1, byrow = TRUE)))
            
            V[[j]] <- matrix(c(Int3[j], -Tend, -Tend, Int1[j]), 2, 2)
            
            SigDelta[j] <- sum(log(delta) + log(X[j, 2:K]))
            
        }
    }
    
    
    return(list(U = U, V = V, S = S, SigDelta = SigDelta))
}
