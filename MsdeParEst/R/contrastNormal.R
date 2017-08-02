#' Computation of the contrast used for the estimation of the normal conditional distribution
#' of the random effects in the drift in mixed SDE with random effects in the drift and in the
#' diffusion coefficient
#' 
#' @description Computation of the contrast used for the estimation of the parameters of the 
#' Gaussian conditional distribution of the random effects in the drift \eqn{\alpha_j,\beta_j} when
#' the SDE includes random effects in the drift and in the diffusion coefficient:
#' 
#'  \eqn{dXj(t)= (\alpha_j - \beta_j Xj(t))dt + \sigma_j a(Xj(t)) dWj(t)}.
#'  
#' @param mu value of the mean of the Gaussian distribution.
#' @param omega value of the standard deviation of the Gaussian distribution.   
#' @param U matrix of M sufficient statistics U (see \code{\link{UVS}}).
#' @param V list of the M sufficient statistics matrix V (see \code{\link{UVS}}).
#' @param S vector of the M sufficient statistics S (see \code{\link{UVS}}).
#' @param K number of times of observations.
#' @param estimphi matrix of the M x 2 estimated parameters (\eqn{\alpha_j,\beta_j}).
#' @param drift.random random effects in the drift: 1 if one additive random effect, 2 if one multiplicative random effect or c(1,2) if 2 random effects.
#' @return
#' \item{L}{value of the contrast}
#' @references
#' Estimaton of the joint distribution of random effects for a discretely observed diffusion with random effects, M. Delattre, V. Genon-Catalot and C. Laredo, \emph{Preprint, hal-01446063}.

contrastNormal <- function(mu, omega, U, V, S, K, estimphi, drift.random) {
    
    if (length(drift.random) == 2) {
        Omega <- matrix(c(omega[1]^2, 0, 0, omega[2]^2), 2, 2, byrow = TRUE)
    }
    
    if (sum(drift.random) == 1) {
        Omega <- matrix(c(omega[1]^2, 0, 0, 0), 2, 2, byrow = TRUE)
    }
    
    if (sum(drift.random) == 2) {
        Omega <- matrix(c(0, 0, 0, omega[2]^2), 2, 2, byrow = TRUE)
    }
    
    M <- dim(U)[2]
    
    loglik <- vector(length = M)
    
    I2 <- diag(c(1, 1))
    
    for (j in 1:M) {
        A <- (I2 + V[[j]] %*% Omega)
        Rinv <- solve(A) %*% V[[j]]
        b <- mu - estimphi[, j]
        loglik[j] <- log(det(A))/2 + K/(2 * S[j]) * (t(b) %*% Rinv %*% b - t(U[, j]) %*% 
            estimphi[, j])
    }
    
    L <- sum(loglik)
    
    return(L)
}
