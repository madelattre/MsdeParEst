\name{msde.sim}
\alias{msde.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
msde.sim(M, T, N = 100, model, drift.random, diffusion.random, density.phi, drift.fixed = 0, drift.param, diffusion.param, nb.mixt = 1, mixt.prop = 1, t0 = 0, X0 = 0.01, invariant = 0, delta = T/N, op.plot = 0, add.plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{T}{
%%     ~~Describe \code{T} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{drift.random}{
%%     ~~Describe \code{drift.random} here~~
}
  \item{diffusion.random}{
%%     ~~Describe \code{diffusion.random} here~~
}
  \item{density.phi}{
%%     ~~Describe \code{density.phi} here~~
}
  \item{drift.fixed}{
%%     ~~Describe \code{drift.fixed} here~~
}
  \item{drift.param}{
%%     ~~Describe \code{drift.param} here~~
}
  \item{diffusion.param}{
%%     ~~Describe \code{diffusion.param} here~~
}
  \item{nb.mixt}{
%%     ~~Describe \code{nb.mixt} here~~
}
  \item{mixt.prop}{
%%     ~~Describe \code{mixt.prop} here~~
}
  \item{t0}{
%%     ~~Describe \code{t0} here~~
}
  \item{X0}{
%%     ~~Describe \code{X0} here~~
}
  \item{invariant}{
%%     ~~Describe \code{invariant} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{op.plot}{
%%     ~~Describe \code{op.plot} here~~
}
  \item{add.plot}{
%%     ~~Describe \code{add.plot} here~~
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
function (M, T, N = 100, model, drift.random, diffusion.random, 
    density.phi, drift.fixed = 0, drift.param, diffusion.param, 
    nb.mixt = 1, mixt.prop = 1, t0 = 0, X0 = 0.01, invariant = 0, 
    delta = T/N, op.plot = 0, add.plot = FALSE) 
{
    con <- file(tempfile(), open = "w")
    on.exit(close(con))
    sde.sim <- function(...) {
        sink(con)
        res <- sde::sde.sim(...)
        sink(NULL)
        res
    }
    if (missing(X0) && missing(invariant)) {
        message("be careful, X0 and invariant are missing thus the initial value X0=0.01 is used")
    }
    delta <- T/N
    times <- seq(t0, T, length = N + 1)
    X <- matrix(0, M, N + 1)
    index <- NA
    if (diffusion.random == 0) {
        sig <- diffusion.param[1]
        psi <- NA
        if ((sum(drift.random) > 2)) {
            phi <- matrix(0, 2, M)
            if (density.phi == "normalnormal") {
                phi[1, ] <- rnorm(M, drift.param[1], drift.param[2])
                phi[2, ] <- rnorm(M, drift.param[3], drift.param[4])
            }
            if (density.phi == "mixture.normal") {
                phi <- mixture.sim(M, "mixture.normal", drift.param, 
                  mixt.prop)
            }
            if (density.phi == "gammagamma") {
                phi[1, ] <- rgamma(M, drift.param[1], scale = drift.param[2])
                phi[2, ] <- rgamma(M, drift.param[3], scale = drift.param[4])
            }
            if (density.phi == "normalgamma") {
                phi[1, ] <- drift.param[1] + drift.param[2] * 
                  rnorm(M, mean = 0, sd = 1)
                phi[2, ] <- rgamma(M, drift.param[3], scale = drift.param[4])
            }
            if (density.phi == "normalinvgamma") {
                phi[1, ] <- drift.param[1] + drift.param[2] * 
                  rnorm(M, mean = 0, sd = 1)
                phi[2, ] <- 1/rgamma(M, drift.param[3], scale = drift.param[4])
            }
            if (density.phi == "gammainvgamma") {
                phi[1, ] <- rgamma(M, drift.param[1], scale = drift.param[2])
                phi[2, ] <- 1/rgamma(M, drift.param[3], scale = drift.param[4])
            }
            if (density.phi == "gammagamma2") {
                phi[1, ] <- 2 * sig^2 + rgamma(M, drift.param[1], 
                  scale = drift.param[2])
                phi[2, ] <- 1 + rgamma(M, drift.param[3], scale = drift.param[4])
            }
            if (density.phi == "gammainvgamma2") {
                phi[1, ] <- 2 * sig^2 + rgamma(M, drift.param[1], 
                  scale = drift.param[2])
                phi[2, ] <- 1/rgamma(M, drift.param[3], scale = drift.param[4])
            }
            if (model == "OU") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- phi[1, j]/phi[2, j] + (sig/(sqrt(2 * 
                      phi[2, j]))) * rnorm(1)
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "EA", 
                      theta = c(phi[, j], sig), model = "OU"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "EA", 
                      theta = c(phi[, j], sig), model = "OU"))
                  }
                }
            }
            if (model == "CIR") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- rgamma(1, 2 * phi[1, j]/sig^2, scale = sig^2/(2 * 
                      phi[2, j]))
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "milstein", 
                      theta = c(phi[, j], sig), model = "CIR", 
                      sigma.x = expression(sig/(2 * sqrt(x))), 
                      sigma = expression(sig * sqrt(x))))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "milstein", 
                      theta = c(phi[, j], sig), sigma.x = expression(sig/(2 * 
                        sqrt(x))), sigma = expression(sig * sqrt(x)), 
                      model = "CIR"))
                  }
                }
            }
        }
        if (sum(drift.random) == 1) {
            phi <- rep(0, M)
            if (density.phi == "mixture.normal") {
                phi <- mixture.sim(M, "mixture.normal", drift.param, 
                  mixt.prop)
            }
            if (density.phi == "normal") {
                phi <- drift.param[1] + drift.param[2] * rnorm(M, 
                  mean = 0, sd = 1)
            }
            if (density.phi == "lognormal") {
                phi <- drift.param[1] + drift.param[2] * rnorm(M, 
                  mean = 0, sd = 1)
                phi <- exp(phi)
            }
            if (density.phi == "gamma") {
                phi <- rgamma(M, shape = drift.param[1], scale = drift.param[2])
            }
            if (density.phi == "gamma2") {
                phi <- 2 * sig^2 + rgamma(M, shape = drift.param[1], 
                  scale = drift.param[2])
            }
            if (density.phi == "mixture.gamma") {
                phi <- mixture.sim(M, density.phi, drift.param, 
                  mixt.prop)
            }
            if (density.phi == "mixture.gamma2") {
                phi <- 2 * sig^2 + mixture.sim(M, density.phi, 
                  drift.param)
            }
            if (model == "OU") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- phi[j]/drift.fixed + (sig/(sqrt(2 * 
                      drift.fixed))) * rnorm(1)
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "EA", 
                      theta = c(phi[j], drift.fixed, sig), model = "OU"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "EA", 
                      theta = c(phi[j], drift.fixed, sig), model = "OU"))
                  }
                }
            }
            if (model == "CIR") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- rgamma(1, 2 * phi[j]/sig^2, scale = sig^2/(2 * 
                      drift.fixed))
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      N = N, X0 = X0, delta = delta, method = "milstein", 
                      theta = c(phi[j], drift.fixed, sig), sigma.x = expression(sig/(2 * 
                        sqrt(x))), sigma = expression(sig * sqrt(x)), 
                      model = "CIR"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      N = N, X0 = X0, delta = delta, method = "milstein", 
                      theta = c(phi[j], drift.fixed, sig), sigma.x = expression(sig/(2 * 
                        sqrt(x))), sigma = expression(sig * sqrt(x)), 
                      model = "CIR"))
                  }
                }
            }
        }
        if (sum(drift.random) == 2) {
            phi <- rep(0, M)
            if (density.phi == "mixture.normal") {
                phi <- mixture.sim(M, "mixture.normal", drift.param, 
                  mixt.prop)
            }
            if (density.phi == "normal") {
                phi <- drift.param[1] + drift.param[2] * rnorm(M, 
                  mean = 0, sd = 1)
            }
            if (density.phi == "gamma") {
                phi <- rgamma(M, drift.param[1], scale = drift.param[2])
            }
            if (density.phi == "mixture.gamma") {
                phi <- mixture.sim(M, density.phi, drift.param, 
                  mixt.prop)
            }
            if (model == "OU") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- (sig/(sqrt(2 * phi[j]))) * rnorm(1)
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "EA", 
                      theta = c(drift.fixed, phi[j], sig), model = "OU"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "EA", 
                      theta = c(drift.fixed, phi[j], sig), model = "OU"))
                  }
                }
            }
            if (model == "CIR") {
                for (j in 1:M) {
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(t0, T, 
                      X0, N, delta, method = "milstein", theta = c(drift.fixed, 
                        phi[j], sig), sigma.x = expression(sig/(2 * 
                        sqrt(x))), sigma = expression(sig * sqrt(x)), 
                      model = "CIR"))
                  }
                  if (invariant == 1) {
                    if (drift.fixed == 0) {
                      message("no invariant distribution, please fix X0")
                    }
                    if (drift.fixed != 0) {
                      X0 <- rgamma(1, 2 * drift.fixed/sig^2, 
                        scale = sig^2/(2 * phi[j]))
                      suppressMessages(X[j, ] <- sde.sim(t0 = t0, 
                        T = T, X0 = X0, N, delta = delta, method = "milstein", 
                        theta = c(drift.fixed, phi[j], sig), 
                        sigma.x = expression(sig/(2 * sqrt(x))), 
                        sigma = expression(sig * sqrt(x)), model = "CIR"))
                    }
                  }
                }
            }
        }
    }
    if (diffusion.random == 1) {
        psi <- 1/sqrt(rgamma(M, shape = diffusion.param[1], rate = 1/diffusion.param[2]))
        if ((sum(drift.random) > 2)) {
            phi <- matrix(0, 2, M)
            if (density.phi == "normalnormal") {
                for (j in 1:M) {
                  phi[1, j] <- rnorm(1, drift.param[1], sd = drift.param[2] * 
                    psi[j])
                  phi[2, j] <- rnorm(1, drift.param[3], sd = drift.param[4] * 
                    psi[j])
                }
            }
            if (model == "OU") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- phi[1, j]/phi[2, j] + (psi[j]/(sqrt(2 * 
                      phi[2, j]))) * rnorm(1)
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "EA", 
                      theta = c(phi[, j], psi[j]), model = "OU"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "EA", 
                      theta = c(phi[, j], psi[j]), model = "OU"))
                  }
                }
            }
            if (model == "CIR") {
                message("to consider a random effect in the diffusion coefficient, alphaj should not be random")
            }
        }
        if (sum(drift.random) == 1) {
            phi <- rep(0, M)
            if (density.phi == "normal") {
                for (j in 1:M) {
                  phi[j] <- drift.param[1] + drift.param[2] * 
                    psi[j] * rnorm(1, mean = 0, sd = 1)
                }
            }
            if (model == "OU") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- phi[j]/drift.fixed + (psi[j]/(sqrt(2 * 
                      drift.fixed))) * rnorm(1)
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "EA", 
                      theta = c(phi[j], drift.fixed, psi[j]), 
                      model = "OU"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "EA", 
                      theta = c(phi[j], drift.fixed, psi[j]), 
                      model = "OU"))
                  }
                }
            }
            if (model == "CIR") {
                message("to consider a random effect in the diffusion coefficient, alphaj should not be random")
            }
        }
        if (sum(drift.random) == 2) {
            phi <- rep(0, M)
            if (density.phi == "normal") {
                for (j in 1:M) {
                  phi[j] <- drift.param[1] + drift.param[2] * 
                    psi[j] * rnorm(1, mean = 0, sd = 1)
                }
            }
            if (model == "OU") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- (psi[j]/(sqrt(2 * phi[j]))) * rnorm(1)
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "EA", 
                      theta = c(drift.fixed, phi[j], psi[j]), 
                      model = "OU"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "EA", 
                      theta = c(drift.fixed, phi[j], psi[j]), 
                      model = "OU"))
                  }
                }
            }
            if (model == "CIR") {
                for (j in 1:M) {
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(t0, T, 
                      X0, N, delta, method = "milstein", theta = c(drift.fixed, 
                        phi[j], psi[j]), sigma.x = expression(psi[j]/(2 * 
                        sqrt(x))), sigma = expression(psi[j] * 
                        sqrt(x)), model = "CIR"))
                  }
                  if (invariant == 1) {
                    if (drift.fixed == 0) {
                      message("no invariant distribution, please fix X0")
                    }
                    if (drift.fixed != 0) {
                      X0 <- rgamma(1, 2 * drift.fixed/psi[j]^2, 
                        scale = psi[j]^2/(2 * phi[j]))
                      suppressMessages(X[j, ] <- sde.sim(t0 = t0, 
                        T = T, X0 = X0, N, delta = delta, method = "milstein", 
                        theta = c(drift.fixed, phi[j], psi[j]), 
                        sigma.x = expression(psi[j]/(2 * sqrt(x))), 
                        sigma = expression(psi[j] * sqrt(x)), 
                        model = "CIR"))
                    }
                  }
                }
            }
        }
        if (sum(drift.random) == 0) {
            phi <- NA
            if (model == "OU") {
                for (j in 1:M) {
                  if (invariant == 1) {
                    X0 <- (psi[j]/(sqrt(2 * drift.fixed[2]))) * 
                      rnorm(1)
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "EA", 
                      theta = c(drift.fixed[1], drift.fixed[2], 
                        psi[j]), model = "OU"))
                  }
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(T = T, 
                      X0 = X0, N = N, delta = delta, method = "EA", 
                      theta = c(drift.fixed[1], drift.fixed[2], 
                        psi[j]), model = "OU"))
                  }
                }
            }
            if (model == "CIR") {
                for (j in 1:M) {
                  if (invariant == 0) {
                    suppressMessages(X[j, ] <- sde.sim(t0, T, 
                      X0, N, delta, method = "milstein", theta = c(drift.fixed[1], 
                        drift.fixed[2], psi[j]), sigma.x = expression(psi[j]/(2 * 
                        sqrt(x))), sigma = expression(psi[j] * 
                        sqrt(x)), model = "CIR"))
                  }
                  if (invariant == 1) {
                    if (drift.fixed[2] == 0) {
                      message("no invariant distribution, please fix X0")
                    }
                    if (drift.fixed[2] != 0) {
                      X0 <- rgamma(1, 2 * drift.fixed[1]/psi[j]^2, 
                        scale = psi[j]^2/(2 * drift.fixed[2]))
                      suppressMessages(X[j, ] <- sde.sim(t0 = t0, 
                        T = T, X0 = X0, N, delta = delta, method = "milstein", 
                        theta = c(drift.fixed[1], drift.fixed[2], 
                          psi[j]), sigma.x = expression(psi[j]/(2 * 
                          sqrt(x))), sigma = expression(psi[j] * 
                          sqrt(x)), model = "CIR"))
                    }
                  }
                }
            }
        }
    }
    output <- X
    if (op.plot) {
        if (add.plot) {
            for (j in 1:M) {
                lines(delta * (0:N), X[j, ], col = j)
            }
        }
        else {
            plot(delta * (0:N), X[1, ], type = "l", ylim = c(min(X), 
                max(X)), xlab = "time", ylab = "", col = 1)
            for (j in 2:M) {
                lines(delta * (0:N), X[j, ], col = j)
            }
        }
    }
    return(list(index = index, phi = phi, psi = psi, X = X, times = times))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
