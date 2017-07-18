\name{EstParamNormalGamma}
\alias{EstParamNormalGamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EstParamNormalGamma(U, V, S, SigDelta, K, drift.random, drift.fixed = 0, estim.drift.fix = 0)
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
  \item{SigDelta}{
%%     ~~Describe \code{SigDelta} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{drift.random}{
%%     ~~Describe \code{drift.random} here~~
}
  \item{drift.fixed}{
%%     ~~Describe \code{drift.fixed} here~~
}
  \item{estim.drift.fix}{
%%     ~~Describe \code{estim.drift.fix} here~~
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
function (U, V, S, SigDelta, K, drift.random, drift.fixed = 0, 
    estim.drift.fix = 0) 
{
    M <- length(S)
    k <- 0.1
    estimGamma <- (K/S) * ((S/K) >= (k/sqrt(K)))
    init.a <- mean(estimGamma)^2/var(estimGamma)
    init.lambda <- var(estimGamma)/mean(estimGamma)
    V1 <- function(param) {
        contrastGamma(param[1], param[2], U, V, S, K, c(0, 0))
    }
    res.gamma <- suppressWarnings(optim(c(init.a, init.lambda), 
        fn = V1, method = "Nelder-Mead"))
    a <- res.gamma$par[1]
    lambda <- res.gamma$par[2]
    estimphi <- matrix(NA, 2, M)
    for (j in 1:M) {
        estimphi[, j] <- solve(V[[j]]) \%*\% U[, j]
    }
    init.mu <- apply(estimphi, 1, mean)
    if (length(drift.random) == 2) {
        init.omega2 <- c(mean(estimGamma * (estimphi[1, ] - mean(estimphi[1, 
            ]))^2), mean(estimGamma * (estimphi[2, ] - mean(estimphi[2, 
            ]))^2))
        V2 <- function(param) {
            contrastNormal(c(param[1], param[2]), c(param[3], 
                param[4]), U, V, S, K, estimphi, drift.random)
        }
        res.normal <- suppressWarnings(optim(c(init.mu, sqrt(init.omega2)), 
            fn = V2, method = "Nelder-Mead"))
        mu <- c(res.normal$par[1], res.normal$par[2])
        omega <- c(abs(res.normal$par[3]), abs(res.normal$par[4]))
        BIChere <- likelihoodNormalGamma(a, lambda, mu, omega, 
            U, V, S, SigDelta, K, drift.random) + 6 * log(M)
        AIChere <- likelihoodNormalGamma(a, lambda, mu, omega, 
            U, V, S, SigDelta, K, drift.random) + 12
    }
    if (estim.drift.fix == 1) {
        if (sum(drift.random) == 1) {
            init.omega2 <- mean(estimGamma * (estimphi[1, ] - 
                mean(estimphi[1, ]))^2)
            V2 <- function(param) {
                contrastNormal(c(param[1], param[2]), c(param[3], 
                  0), U = U, V = V, S = S, K = K, estimphi = estimphi, 
                  drift.random = drift.random)
            }
            res.normal <- suppressWarnings(optim(c(init.mu, sqrt(init.omega2)), 
                fn = V2, method = "Nelder-Mead"))
            mu <- c(res.normal$par[1], res.normal$par[2])
            omega <- c(abs(res.normal$par[3]), 0)
            BIChere <- likelihoodNormalGamma(a, lambda, mu, omega, 
                U, V, S, SigDelta, K, drift.random) + 5 * log(M)
            AIChere <- likelihoodNormalGamma(a, lambda, mu, omega, 
                U, V, S, SigDelta, K, drift.random) + 10
        }
        if (sum(drift.random) == 2) {
            init.omega2 <- mean(estimGamma * (estimphi[2, ] - 
                mean(estimphi[2, ]))^2)
            V2 <- function(param) {
                contrastNormal(mu = c(param[1], param[2]), omega = c(0, 
                  param[3]), U = U, V = V, S = S, K = K, estimphi = estimphi, 
                  drift.random = drift.random)
            }
            res.normal <- suppressWarnings(optim(c(init.mu, sqrt(init.omega2)), 
                fn = V2, method = "Nelder-Mead"))
            mu <- c(res.normal$par[1], res.normal$par[2])
            omega <- c(0, abs(res.normal$par[3]))
            BIChere <- likelihoodNormalGamma(a, lambda, mu, omega, 
                U, V, S, SigDelta, K, drift.random) + 5 * log(M)
            AIChere <- likelihoodNormalGamma(a, lambda, mu, omega, 
                U, V, S, SigDelta, K, drift.random) + 10
        }
    }
    if (estim.drift.fix == 0) {
        if (sum(drift.random) == 1) {
            init.omega2 <- mean(estimGamma * (estimphi[1, ] - 
                mean(estimphi[1, ]))^2)
            V2 <- function(param) {
                contrastNormal(c(param[1], drift.fixed), c(param[2], 
                  0), U = U, V = V, S = S, K = K, estimphi = estimphi, 
                  drift.random = drift.random)
            }
            res.normal <- suppressWarnings(optim(c(init.mu[1], 
                sqrt(init.omega2)), fn = V2, method = "Nelder-Mead"))
            mu <- c(res.normal$par[1], drift.fixed)
            omega <- c(abs(res.normal$par[2]), 0)
            BIChere <- likelihoodNormalGamma(a, lambda, mu, omega, 
                U, V, S, SigDelta, K, drift.random) + 5 * log(M)
            AIChere <- likelihoodNormalGamma(a, lambda, mu, omega, 
                U, V, S, SigDelta, K, drift.random) + 10
        }
        if (sum(drift.random) == 2) {
            init.omega2 <- mean(estimGamma * (estimphi[2, ] - 
                mean(estimphi[2, ]))^2)
            V2 <- function(param) {
                contrastNormal(mu = c(drift.fixed, param[1]), 
                  omega = c(0, param[2]), U = U, V = V, S = S, 
                  K = K, estimphi = estimphi, drift.random = drift.random)
            }
            res.normal <- suppressWarnings(optim(c(init.mu, sqrt(init.omega2)), 
                fn = V2, method = "Nelder-Mead"))
            mu <- c(drift.fixed, res.normal$par[1])
            omega <- c(0, abs(res.normal$par[2]))
            BIChere <- likelihoodNormalGamma(a, lambda, mu, omega, 
                U, V, S, SigDelta, K, drift.random) + 5 * log(M)
            AIChere <- likelihoodNormalGamma(a, lambda, mu, omega, 
                U, V, S, SigDelta, K, drift.random) + 10
        }
    }
    return(list(mu = mu, omega = omega, a = a, lambda = lambda, 
        BIChere = BIChere, AIChere = AIChere))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
