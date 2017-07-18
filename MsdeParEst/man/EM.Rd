\name{EM}
\alias{EM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EM(U, V, S, K, drift.random, start, Niter = 10, drift.fixed = 0, drift.estim.fixed = 1, drift.fixed.mixt = 1, sigma = 1, sigma.estim = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{U}{
%%     ~~Describe \code{U} here~~
}
  \item{V}{
%%     ~~Describe \code{V} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{drift.random}{
%%     ~~Describe \code{drift.random} here~~
}
  \item{start}{
%%     ~~Describe \code{start} here~~
}
  \item{Niter}{
%%     ~~Describe \code{Niter} here~~
}
  \item{drift.fixed}{
%%     ~~Describe \code{drift.fixed} here~~
}
  \item{drift.estim.fixed}{
%%     ~~Describe \code{drift.estim.fixed} here~~
}
  \item{drift.fixed.mixt}{
%%     ~~Describe \code{drift.fixed.mixt} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{sigma.estim}{
%%     ~~Describe \code{sigma.estim} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (U, V, S, K, drift.random, start, Niter = 10, drift.fixed = 0, 
    drift.estim.fixed = 1, drift.fixed.mixt = 1, sigma = 1, sigma.estim = 1) 
{
    M <- dim(U)[2]
    mu <- start$mu
    omega <- start$omega
    mixt.prop <- start$mixt.prop
    N <- length(mixt.prop)
    estimphi <- matrix(NA, 2, M)
    for (j in 1:M) {
        estimphi[, j] <- solve(V[[j]]) \%*\% U[, j]
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
        probindi <- probind(mu, omega, mixt.prop, sigma, U, V, 
            S, K, estimphi, drift.random)
        mixt.prop <- colMeans(probindi, na.rm = TRUE)
        probhat[iter, ] <- mixt.prop
        mu <- muhat[iter - 1, , ]
        omega <- omegahat[iter - 1, , ]
        if (sum(drift.random) > 2) {
            ln <- function(param) (Q_EM(matrix(param[1:(2 * N)], 
                nrow = N, ncol = 2), matrix(param[(2 * N + 1):(4 * 
                N)], nrow = N, ncol = 2), sigma, probindi, U, 
                V, S, K, estimphi, drift.random))
            paraminit <- c(as.vector(mu), as.vector(omega))
            res <- optim(paraminit, f = ln, method = "Nelder-Mead")
            muhat[iter, , ] <- matrix(res$par[1:(2 * N)], nrow = N, 
                ncol = 2)
            omegahat[iter, , ] <- matrix(abs(res$par[(2 * N + 
                1):(4 * N)]), nrow = N, ncol = 2)
        }
        if (drift.estim.fixed == 1) {
            if (drift.fixed.mixt == 0) {
                if (sum(drift.random) == 2) {
                  ln <- function(param) (Q_EM(matrix(c(rep(param[1], 
                    N), param[2:(N + 1)]), nrow = N, ncol = 2), 
                    matrix(c(rep(0, N), param[(N + 2):(2 * N + 
                      1)]), nrow = N, ncol = 2), sigma, probindi, 
                    U, V, S, K, estimphi, drift.random))
                  paraminit <- c(mu[1, 1], as.vector(mu[, 2]), 
                    as.vector(omega[, 2]))
                  res <- optim(paraminit, f = ln, method = "Nelder-Mead")
                  muhat[iter, , ] <- matrix(c(rep(res$par[1], 
                    N), res$par[2:(N + 1)]), nrow = N, ncol = 2)
                  omegahat[iter, , ] <- matrix(c(rep(0, N), abs(res$par[(N + 
                    2):(2 * N + 1)])), nrow = N, ncol = 2)
                }
                if (sum(drift.random) == 1) {
                  ln <- function(param) (Q_EM(matrix(c(param[1:N], 
                    rep(param[N + 1], N)), nrow = N, ncol = 2), 
                    matrix(c(param[(N + 2):(2 * N + 1)], rep(0, 
                      N)), nrow = N, ncol = 2), sigma, probindi, 
                    U, V, S, K, estimphi, drift.random))
                  paraminit <- c(as.vector(mu[, 1]), mu[2, 1], 
                    as.vector(omega[, 1]))
                  res <- optim(paraminit, f = ln, method = "Nelder-Mead")
                  muhat[iter, , ] <- matrix(c(res$par[1:N], rep(res$par[N + 
                    1], N)), nrow = N, ncol = 2)
                  omegahat[iter, , ] <- matrix(c(abs(res$par[(N + 
                    2):(2 * N + 1)]), rep(0, N)), nrow = N, ncol = 2)
                }
            }
            if (drift.fixed.mixt == 1) {
                if (sum(drift.random) == 2) {
                  ln <- function(param) (Q_EM(matrix(param[1:(2 * 
                    N)], nrow = N, ncol = 2), matrix(c(rep(0, 
                    N), param[(2 * N + 1):(3 * N)]), nrow = N, 
                    ncol = 2), sigma, probindi, U, V, S, K, estimphi, 
                    drift.random))
                  paraminit <- c(as.vector(mu), as.vector(omega[, 
                    2]))
                  res <- optim(paraminit, f = ln, method = "Nelder-Mead")
                  muhat[iter, , ] <- matrix(c(res$par[1:N], res$par[(N + 
                    1):(2 * N)]), nrow = N, ncol = 2)
                  omegahat[iter, , ] <- matrix(c(rep(0, N), abs(res$par[(2 * 
                    N + 1):(3 * N)])), nrow = N, ncol = 2)
                }
                if (sum(drift.random) == 1) {
                  ln <- function(param) (Q_EM(matrix(param[1:(2 * 
                    N)], nrow = N, ncol = 2), matrix(c(param[(2 * 
                    N + 1):(3 * N)], rep(0, N)), nrow = N, ncol = 2), 
                    sigma, probindi, U, V, S, K, estimphi, drift.random))
                  paraminit <- c(as.vector(mu), as.vector(omega[, 
                    1]))
                  res <- optim(paraminit, f = ln, method = "Nelder-Mead")
                  muhat[iter, , ] <- matrix(res$par[1:(2 * N)], 
                    nrow = N, ncol = 2)
                  omegahat[iter, , ] <- matrix(c(abs(res$par[(2 * 
                    N + 1):(3 * N)]), rep(0, N)), nrow = N, ncol = 2)
                }
            }
        }
        if (drift.estim.fixed == 0) {
            if (sum(drift.random) == 2) {
                ln <- function(param) (Q_EM(matrix(c(drift.fixed, 
                  param[1:N]), nrow = N, ncol = 2), matrix(c(rep(0, 
                  N), param[(N + 1):(2 * N)]), nrow = N, ncol = 2), 
                  sigma, probindi, U, V, S, K, estimphi, drift.random))
                paraminit <- c(as.vector(mu[, 2]), as.vector(omega[, 
                  2]))
                res <- optim(paraminit, f = ln, method = "Nelder-Mead")
                muhat[iter, , ] <- matrix(c(drift.fixed, res$par[1:N]), 
                  nrow = N, ncol = 2)
                omegahat[iter, , ] <- matrix(c(rep(0, N), abs(res$par[(N + 
                  1):(2 * N)])), nrow = N, ncol = 2)
            }
            if (sum(drift.random) == 1) {
                ln <- function(param) (Q_EM(matrix(c(param[1:N], 
                  drift.fixed), nrow = N, ncol = 2), matrix(c(param[(N + 
                  1):(2 * N)], rep(0, N)), nrow = N, ncol = 2), 
                  sigma, probindi, U, V, S, K, estimphi, drift.random))
                paraminit <- c(as.vector(mu[, 1]), as.vector(omega[, 
                  1]))
                res <- optim(paraminit, f = ln, method = "Nelder-Mead")
                muhat[iter, , ] <- matrix(c(res$par[1:N], drift.fixed), 
                  nrow = N, ncol = 2)
                omegahat[iter, , ] <- matrix(c(abs(res$par[(N + 
                  1):(2 * N)]), rep(0, N)), nrow = N, ncol = 2)
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
    BIChere <- likelihoodMixtureNormal(muhat[Niter, , ], omegahat[Niter, 
        , ], sigma, mixt.prop, U, V, S, K, estimphi, drift.random) + 
        nbparam * log(M)
    AIChere <- likelihoodMixtureNormal(muhat[Niter, , ], omegahat[Niter, 
        , ], sigma, mixt.prop, U, V, S, K, estimphi, drift.random) + 
        nbparam * 2
    output = list(mu = muhat, omega = omegahat * sigma, mixt.prop = probhat, 
        sigma = sigma, probindi = probindi, BIChere = BIChere, 
        AIChere = AIChere)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
