\name{EstParamNormal}
\alias{EstParamNormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EstParamNormal(U, V, S, SigDelta = 0, K, drift.fixed = 0, estim.drift.fix = 1, sigma = 0, drift.random, diffusion.estim = 1, discrete = 1)
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
  \item{drift.fixed}{
%%     ~~Describe \code{drift.fixed} here~~
}
  \item{estim.drift.fix}{
%%     ~~Describe \code{estim.drift.fix} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{drift.random}{
%%     ~~Describe \code{drift.random} here~~
}
  \item{diffusion.estim}{
%%     ~~Describe \code{diffusion.estim} here~~
}
  \item{discrete}{
%%     ~~Describe \code{discrete} here~~
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
function (U, V, S, SigDelta = 0, K, drift.fixed = 0, estim.drift.fix = 1, 
    sigma = 0, drift.random, diffusion.estim = 1, discrete = 1) 
{
    M <- length(S)
    estimphi <- matrix(NA, 2, M)
    for (j in 1:M) {
        estimphi[, j] <- solve(V[[j]]) \%*\% U[, j]
    }
    init.mu <- apply(estimphi, 1, mean)
    if (discrete == 1) {
        if (diffusion.estim == 1) {
            k <- 0.1
            estimGamma <- (K/S) * ((S/K) >= (k/sqrt(K)))
            init.gamma <- mean(estimGamma)
            init.sigma <- 1/sqrt(init.gamma)
            if (length(drift.random) == 2) {
                init.omega2 <- apply(estimphi, 1, var) * init.gamma
                ln = function(param) {
                  likelihoodNormal(c(param[1], param[2]), c(param[3], 
                    param[4]), exp(param[5]), U, V, S, SigDelta, 
                    K, estimphi, drift.random, discrete)
                }
                res <- optim(c(init.mu, sqrt(init.omega2), init.sigma), 
                  fn = ln, method = "Nelder-Mead")
                mu <- c(res$par[1], res$par[2])
                sigma <- exp(res$par[5])
                omega <- c(abs(res$par[3]), abs(res$par[4])) * 
                  sigma
                nbparam <- 5
            }
            if (sum(drift.random) == 1) {
                init.omega2 <- var(estimphi[1, ]) * init.gamma
                if (estim.drift.fix == 1) {
                  ln = function(param) {
                    likelihoodNormal(c(param[1], param[2]), c(param[3], 
                      0), exp(param[4]), U, V, S, SigDelta, K, 
                      estimphi, drift.random, discrete)
                  }
                  res <- optim(c(init.mu, sqrt(init.omega2), 
                    init.sigma), fn = ln, method = "Nelder-Mead")
                  mu <- c(res$par[1], res$par[2])
                  sigma <- exp(res$par[4])
                  omega <- c(abs(res$par[3]), 0) * sigma
                  nbparam <- 4
                }
                if (estim.drift.fix == 0) {
                  ln = function(param) {
                    likelihoodNormal(c(param[1], drift.fixed), 
                      c(param[2], 0), exp(param[3]), U, V, S, 
                      SigDelta, K, estimphi, drift.random, discrete)
                  }
                  res <- optim(c(init.mu[1], sqrt(init.omega2), 
                    init.sigma), fn = ln, method = "Nelder-Mead")
                  mu <- c(res$par[1], drift.fixed)
                  sigma <- exp(res$par[3])
                  omega <- c(abs(res$par[2]), 0) * sigma
                  nbparam <- 3
                }
            }
            if (sum(drift.random) == 2) {
                init.omega2 <- var(estimphi[2, ]) * init.gamma
                if (estim.drift.fix == 1) {
                  ln = function(param) {
                    likelihoodNormal(c(param[1], param[2]), c(0, 
                      param[3]), exp(param[4]), U, V, S, SigDelta, 
                      K, estimphi, drift.random, discrete)
                  }
                  res <- optim(c(init.mu, sqrt(init.omega2), 
                    init.sigma), fn = ln, method = "Nelder-Mead")
                  mu <- c(res$par[1], res$par[2])
                  sigma <- exp(res$par[4])
                  omega <- c(0, abs(res$par[3])) * sigma
                  nbparam <- 4
                }
                if (estim.drift.fix == 0) {
                  ln = function(param) {
                    likelihoodNormal(c(drift.fixed, param[1]), 
                      c(0, param[2]), exp(param[3]), U, V, S, 
                      SigDelta, K, estimphi, drift.random, discrete)
                  }
                  res <- optim(c(init.mu[2], sqrt(init.omega2), 
                    init.sigma), fn = ln, method = "Nelder-Mead")
                  mu <- c(drift.fixed, res$par[1])
                  sigma <- exp(res$par[3])
                  omega <- c(0, abs(res$par[2])) * sigma
                  nbparam <- 3
                }
            }
        }
        if (diffusion.estim == 0) {
            init.mu <- apply(estimphi, 1, mean)
            if (length(drift.random) == 2) {
                init.omega2 <- apply(estimphi, 1, var)/sigma^2
                ln = function(param) {
                  likelihoodNormal(c(param[1], param[2]), c(param[3], 
                    param[4]), sigma, U, V, S, SigDelta, K, estimphi, 
                    drift.random, discrete)
                }
                res <- optim(c(init.mu, sqrt(init.omega2)), fn = ln, 
                  method = "Nelder-Mead")
                mu <- c(res$par[1], res$par[2])
                omega <- c(abs(res$par[3]), abs(res$par[4])) * 
                  sigma
                nbparam <- 4
            }
            if (sum(drift.random) == 1) {
                init.omega2 <- var(estimphi[1, ])/sigma^2
                if (estim.drift.fix == 1) {
                  ln = function(param) {
                    likelihoodNormal(c(param[1], param[2]), c(param[3], 
                      0), sigma, U, V, S, SigDelta, K, estimphi, 
                      drift.random, discrete)
                  }
                  res <- optim(c(init.mu, sqrt(init.omega2)), 
                    fn = ln, method = "Nelder-Mead")
                  mu <- c(res$par[1], res$par[2])
                  omega <- c(abs(res$par[3]), 0) * sigma
                  nbparam <- 3
                }
                if (estim.drift.fix == 0) {
                  ln = function(param) {
                    likelihoodNormal(c(param[1], drift.fixed), 
                      c(param[2], 0), sigma, U, V, S, SigDelta, 
                      K, estimphi, drift.random, discrete)
                  }
                  res <- optim(c(init.mu[1], sqrt(init.omega2)), 
                    fn = ln, method = "Nelder-Mead")
                  mu <- c(res$par[1], drift.fixed)
                  omega <- c(abs(res$par[2]), 0) * sigma
                  nbparam <- 2
                }
            }
            if (sum(drift.random) == 2) {
                init.omega2 <- var(estimphi[2, ])/sigma^2
                if (estim.drift.fix == 1) {
                  ln = function(param) {
                    likelihoodNormal(c(param[1], param[2]), c(0, 
                      param[3]), sigma, U, V, S, SigDelta, K, 
                      estimphi, drift.random, discrete)
                  }
                  res <- optim(c(init.mu, sqrt(init.omega2)), 
                    fn = ln, method = "Nelder-Mead")
                  mu <- c(res$par[1], res$par[2])
                  omega <- c(0, abs(res$par[3])) * sigma
                  nbparam <- 3
                }
                if (estim.drift.fix == 0) {
                  ln = function(param) {
                    likelihoodNormal(c(drift.fixed, param[1]), 
                      c(0, param[2]), sigma, U, V, S, SigDelta, 
                      K, estimphi, drift.random, discrete)
                  }
                  res <- optim(c(init.mu[2], sqrt(init.omega2)), 
                    fn = ln, method = "Nelder-Mead")
                  mu <- c(drift.fixed, res$par[1])
                  omega <- c(0, abs(res$par[2])) * sigma
                  nbparam <- 2
                }
            }
        }
    }
    if (discrete == 0) {
        if (diffusion.estim == 1) {
            sigma <- sqrt(mean(S/(K - 1)))
        }
        if (length(drift.random) == 2) {
            init.omega2 <- apply(estimphi, 1, var)/sigma^2
            ln = function(param) {
                likelihoodNormal(c(param[1], param[2]), c(param[3], 
                  param[4]), sigma, U, V, S, SigDelta, K, estimphi, 
                  drift.random, discrete)
            }
            res <- optim(c(init.mu, sqrt(init.omega2)), fn = ln, 
                method = "Nelder-Mead")
            mu <- c(res$par[1], res$par[2])
            omega <- c(abs(res$par[3]), abs(res$par[4])) * sigma
            nbparam <- 4 + diffusion.estim
        }
        if (sum(drift.random) == 1) {
            init.omega2 <- var(estimphi[1, ])/sigma^2
            if (estim.drift.fix == 1) {
                ln = function(param) {
                  likelihoodNormal(c(param[1], param[2]), c(param[3], 
                    0), sigma, U, V, S, SigDelta, K, estimphi, 
                    drift.random, discrete)
                }
                res <- optim(c(init.mu, sqrt(init.omega2)), fn = ln, 
                  method = "Nelder-Mead")
                mu <- c(res$par[1], res$par[2])
                omega <- c(abs(res$par[3]), 0) * sigma
                nbparam <- 3 + diffusion.estim
            }
            if (estim.drift.fix == 0) {
                ln = function(param) {
                  likelihoodNormal(c(param[1], drift.fixed), 
                    c(param[2], 0), sigma, U, V, S, SigDelta, 
                    K, estimphi, drift.random, discrete)
                }
                res <- optim(c(init.mu[1], sqrt(init.omega2)), 
                  fn = ln, method = "Nelder-Mead")
                mu <- c(res$par[1], drift.fixed)
                omega <- c(abs(res$par[2]), 0) * sigma
                nbparam <- 2 + diffusion.estim
            }
        }
        if (sum(drift.random) == 2) {
            init.omega2 <- var(estimphi[2, ])/sigma^2
            if (estim.drift.fix == 1) {
                ln = function(param) {
                  likelihoodNormal(c(param[1], param[2]), c(0, 
                    param[3]), sigma, U, V, S, SigDelta, K, estimphi, 
                    drift.random, discrete)
                }
                res <- optim(c(init.mu, sqrt(init.omega2)), fn = ln, 
                  method = "Nelder-Mead")
                mu <- c(res$par[1], res$par[2])
                omega <- c(0, abs(res$par[3])) * sigma
                nbparam <- 3 + diffusion.estim
            }
            if (estim.drift.fix == 0) {
                ln = function(param) {
                  likelihoodNormal(c(drift.fixed, param[1]), 
                    c(0, param[2]), sigma, U, V, S, SigDelta, 
                    K, estimphi, drift.random, discrete)
                }
                res <- optim(c(init.mu[2], sqrt(init.omega2)), 
                  fn = ln, method = "Nelder-Mead")
                mu <- c(drift.fixed, res$par[1])
                omega <- c(0, abs(res$par[2])) * sigma
                nbparam <- 2 + diffusion.estim
            }
        }
    }
    BIChere <- likelihoodNormal(mu, omega, sigma, U, V, S, SigDelta, 
        K, estimphi, drift.random, discrete) + nbparam * log(M)
    AIChere <- likelihoodNormal(mu, omega, sigma, U, V, S, SigDelta, 
        K, estimphi, drift.random, discrete) + nbparam * 2
    return(list(mu = mu, omega = omega, sigma = sigma, BIChere = BIChere, 
        AIChere = AIChere))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
