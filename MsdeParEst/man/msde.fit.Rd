\name{msde.fit}
\alias{msde.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
msde.fit(times, X, model = c("OU", "CIR"), drift.random, drift.fixed = NULL, estim.drift.fix = 0, diffusion.random = 0, diffusion.fixed = NULL, estim.diffusion.fix = 0, estim.method = c("paramML", "paramMLmixture"), nb.mixt = 1, drift.fixed.mixed = 0, Niter = 10, discrete = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
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
  \item{diffusion.random}{
%%     ~~Describe \code{diffusion.random} here~~
}
  \item{diffusion.fixed}{
%%     ~~Describe \code{diffusion.fixed} here~~
}
  \item{estim.diffusion.fix}{
%%     ~~Describe \code{estim.diffusion.fix} here~~
}
  \item{estim.method}{
%%     ~~Describe \code{estim.method} here~~
}
  \item{nb.mixt}{
%%     ~~Describe \code{nb.mixt} here~~
}
  \item{drift.fixed.mixed}{
%%     ~~Describe \code{drift.fixed.mixed} here~~
}
  \item{Niter}{
%%     ~~Describe \code{Niter} here~~
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
function (times, X, model = c("OU", "CIR"), drift.random, drift.fixed = NULL, 
    estim.drift.fix = 0, diffusion.random = 0, diffusion.fixed = NULL, 
    estim.diffusion.fix = 0, estim.method = c("paramML", "paramMLmixture"), 
    nb.mixt = 1, drift.fixed.mixed = 0, Niter = 10, discrete = 1) 
{
    model <- match.arg(model)
    estim.method <- match.arg(estim.method)
    if (is.matrix(X)) {
        if (nrow(X) == length(times)) {
            X <- t(X)
        }
        else {
            if (ncol(X) != length(times)) {
                stop("Length of times has to be equal to the columns of X")
            }
        }
    }
    M <- dim(X)[1]
    K <- dim(X)[2]
    delta <- round(diff(times), 10)
    Tend <- times[length(times)]
    if (estim.method == "paramML") {
        a <- 0
        lambda <- 0
        mu <- 0
        omega <- 0
        sigma2 <- 0
        gridf <- NULL
        estimg <- 0
        gridg <- 0
        estimf <- 0
        gridf <- 0
        if (model == "OU") {
            Mindex <- M
            index <- 1:Mindex
        }
        if (model == "CIR") {
            index <- which(rowSums(X <= 0) == 0)
            Mindex <- length(index)
            if (Mindex == 0) {
                warning("All the trajectories have non positive values the model CIR cannot be used", 
                  call. = FALSE)
                estimf <- 0
                estimphi <- 0
                estimpsi2 <- 0
                bic <- 0
                aic <- 0
                gridf <- 0
                mu <- 0
                omega <- 0
                cutoff <- 0
                sigma2 <- 0
                a <- 0
                lambda <- 0
                estimf.trunc <- 0
                estimphi.trunc <- 0
                estimpsi2.trunc <- 0
            }
        }
        if (Mindex > 0) {
            if (((sum(drift.random)) == 0) && (diffusion.random == 
                0)) {
                stop("There should be at least one random effect either in the drift or the diffusion coefficient.")
            }
            if (diffusion.random == 0) {
                if ((is.null(diffusion.fixed)) && (estim.diffusion.fix == 
                  0)) {
                  warning("estim.diffusion.fix == 0. No value is specified for the diffusion coefficient. It is therefore estimated.")
                  estim.diffusion.fix <- 1
                }
                if ((!is.null(diffusion.fixed)) && (estim.diffusion.fix == 
                  1)) {
                  warning("estim.diffusion.fix == 1. The diffusion coefficient is estimated. The specified value is therefore not used.")
                }
            }
            if ((is.null(drift.fixed)) && (estim.drift.fix == 
                0)) {
                warning("estim.drift.fix == 0. No value is specified for the drift parameters. They are therefore estimated.")
                estim.drift.fix <- 1
            }
            if ((!is.null(drift.fixed)) && (estim.drift.fix == 
                1)) {
                warning("estim.drift.fix == 1. The drift parameters are estimated. The specified values are therefore not used.")
            }
            U <- matrix(0, 2, Mindex)
            V <- as.list(1:Mindex)
            S <- rep(0, Mindex)
            SigDelta <- rep(0, Mindex)
            estimUV <- UVS(X[index, ], model, times)
            U <- estimUV$U
            V <- estimUV$V
            S <- estimUV$S
            SigDelta <- estimUV$SigDelta
            deter <- lapply(V, det)
            index2 <- which((deter != Inf) & (deter != 0))
            Mindex2 <- length(index2)
            V <- V[index2]
            U <- U[, index2]
            S <- S[index2]
            SigDelta <- SigDelta[index2]
            estimphi <- matrix(0, 2, Mindex2)
            for (j in 1:Mindex2) {
                estimphi[, j] <- solve(V[[j]]) \%*\% U[, j]
            }
            estimpsi2 <- rep(0, Mindex2)
            estimpsi2 <- S/K
            if (((sum(drift.random)) >= 1) && (diffusion.random == 
                0)) {
                res <- EstParamNormal(U, V, S, SigDelta, K, drift.fixed, 
                  estim.drift.fix, diffusion.fixed, drift.random, 
                  estim.diffusion.fix, discrete)
                bic <- res$BIChere
                aic <- res$AIChere
                mu <- res$mu
                omega <- res$omega
                sigma2 <- res$sigma^2
                gridg <- as.matrix(gridg)
                estimg <- as.matrix(estimg)
                estimpsi2 <- matrix(estimpsi2, 1, length(estimpsi2), 
                  byrow = TRUE)
                if (sum(drift.random) == 3) {
                  if (is.null(gridf) == 1) {
                    gridf <- matrix(0, 2, 500)
                    gridf[1, ] <- seq(mu[1] - 3 * omega[1], mu[1] + 
                      3 * omega[1], length = 500)
                    gridf[2, ] <- seq(mu[2] - 3 * omega[2], mu[2] + 
                      3 * omega[2], length = 500)
                  }
                  if (is.null(gridf) == 0) {
                    gridf <- gridf
                  }
                  estimf1 <- dnorm(gridf[1, ], mean = mu[1], 
                    sd = abs(omega[1]))
                  estimf2 <- dnorm(gridf[2, ], mean = mu[2], 
                    sd = abs(omega[2]))
                  estimf <- rbind(estimf1, estimf2)
                  estimf.trunc <- estimf
                  estimphi.trunc <- estimphi
                  cutoff <- FALSE
                }
                if (sum(drift.random) == 2) {
                  estimphi <- estimphi[2, ]
                  if (is.null(gridf) == 1) {
                    gridf <- seq(mu[2] - 3 * omega[2], mu[2] + 
                      3 * omega[2], length = 500)
                  }
                  if (is.null(gridf) == 0) {
                    gridf <- gridf
                  }
                  estimf <- matrix(dnorm(gridf, mean = mu[2], 
                    sd = omega[2]), 1, length(gridf), byrow = TRUE)
                  gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
                  estimphi <- matrix(estimphi, 1, length(estimphi), 
                    byrow = TRUE)
                  estimphi.trunc <- estimphi
                  estimf.trunc <- estimf
                  cutoff <- FALSE
                }
                if (sum(drift.random) == 1) {
                  estimphi <- estimphi[1, ]
                  if (is.null(gridf) == 1) {
                    gridf <- seq(mu[1] - 3 * omega[1], mu[1] + 
                      3 * omega[1], length = 500)
                  }
                  if (is.null(gridf) == 0) {
                    gridf <- gridf
                  }
                  estimf <- matrix(dnorm(gridf, mean = mu[1], 
                    sd = omega[1]), 1, length(gridf), byrow = TRUE)
                  gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
                  estimphi <- matrix(estimphi, 1, length(estimphi), 
                    byrow = TRUE)
                  estimphi.trunc <- estimphi
                  estimf.trunc <- estimf
                  cutoff <- FALSE
                }
            }
            if ((sum(drift.random) == 0) && (diffusion.random == 
                1)) {
                if (discrete == 0) {
                  warning("Estimation cannot be performed from the likelihood associated with a continuous observation of the trajectories.")
                }
                res <- EstParamGamma(U, V, S, SigDelta, K, drift.fixed, 
                  estim.drift.fix)
                bic <- res$BIChere
                aic <- res$AIChere
                mu <- res$mu
                a <- res$a
                lambda <- res$lambda
                estimphi <- 0
                if (is.null(gridg) == 1) {
                  gridg <- seq(min(estimpsi2) * 0.8, max(estimpsi2) * 
                    1.2, length = 500)
                }
                if (is.null(gridg) == 0) {
                  gridg <- gridg
                }
                estimg <- matrix(dinvgamma(gridg, a, scale = lambda), 
                  1, length(gridg), byrow = TRUE)
                gridg <- matrix(gridg, 1, length(gridg), byrow = TRUE)
                estimpsi2 <- matrix(estimpsi2, 1, length(estimpsi2), 
                  byrow = TRUE)
                cutoff <- FALSE
                gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
                estimf <- matrix(estimf, 1, length(estimf), byrow = TRUE)
                estimphi <- matrix(estimphi, 1, length(estimphi), 
                  byrow = TRUE)
            }
            if (((sum(drift.random)) >= 1) && (diffusion.random == 
                1)) {
                res <- EstParamNormalGamma(U, V, S, SigDelta, 
                  K, drift.random, drift.fixed, estim.drift.fix)
                bic <- res$BIChere
                aic <- res$AIChere
                mu <- res$mu
                omega <- res$omega
                a <- res$a
                lambda <- res$lambda
                estimpsi2 <- matrix(estimpsi2, 1, length(estimpsi2), 
                  byrow = TRUE)
                simupsi2 <- 1/rgamma(500, a, rate = 1/lambda)
                gridg <- seq(min(estimpsi2) * 0.8, max(estimpsi2) * 
                  1.2, length = 500)
                estimg <- matrix(dinvgamma(gridg, a, rate = 1/lambda), 
                  nrow = 1)
                gridg <- matrix(gridg, nrow = 1)
                if (sum(drift.random) == 3) {
                  simuphi <- matrix(NA, 2, length(simupsi2))
                  simuphi[1, ] <- sapply(1:length(simupsi2), 
                    function(s) {
                      rnorm(1, mean = mu[1], sd = abs(omega[1] * 
                        sqrt(simupsi2[s])))
                    })
                  simuphi[2, ] <- sapply(1:length(simupsi2), 
                    function(s) {
                      rnorm(1, mean = mu[2], sd = abs(omega[2] * 
                        sqrt(simupsi2[s])))
                    })
                  test1 <- density(simuphi[1, ], from = min(simuphi[1, 
                    ]), to = max(simuphi[1, ]), bw = "ucv", n = 500)
                  test2 <- density(simuphi[2, ], from = min(simuphi[2, 
                    ]), to = max(simuphi[2, ]), bw = "ucv", n = 500)
                  if (test1$bw < 0.1) {
                    test1 <- density(simuphi[1, ], from = min(simuphi[1, 
                      ]), to = max(simuphi[1, ]), n = 500)
                  }
                  estimf1 <- test1$y
                  gridf1 <- test1$x
                  if (test2$bw < 0.1) {
                    test2 <- density(simuphi[2, ], from = min(simuphi[2, 
                      ]), to = max(simuphi[2, ]), n = 500)
                  }
                  estimf2 <- test2$y
                  gridf2 <- test2$x
                  gridf <- matrix(0, 2, 500)
                  gridf[1, ] <- gridf1
                  gridf[2, ] <- gridf2
                  estimf <- matrix(0, 2, 500)
                  estimf[1, ] <- estimf1
                  estimf[2, ] <- estimf2
                  estimf.trunc <- estimf
                  estimphi.trunc <- estimphi
                  cutoff <- FALSE
                }
                if (sum(drift.random) == 2) {
                  simuphi <- rep(NA, length(simupsi2))
                  simuphi <- sapply(1:length(simupsi2), function(s) {
                    rnorm(1, mean = mu[2], sd = abs(omega[2] * 
                      sqrt(simupsi2[s])))
                  })
                  test <- density(simuphi, from = min(simuphi), 
                    to = max(simuphi), bw = "ucv", n = 500)
                  if (test$bw < 0.1) {
                    test <- density(simuphi, from = min(simuphi), 
                      to = max(simuphi), n = 500)
                  }
                  gridf <- test$x
                  estimf <- test$y
                  estimf <- matrix(estimf, nrow = 1)
                  gridf <- matrix(gridf, nrow = 1)
                  estimf.trunc <- estimf
                  estimphi.trunc <- estimphi
                  cutoff <- FALSE
                }
                if (sum(drift.random) == 1) {
                  simuphi <- rep(NA, length(simupsi2))
                  simuphi <- sapply(1:length(simupsi2), function(s) {
                    rnorm(1, mean = mu[1], sd = abs(omega[1] * 
                      sqrt(simupsi2[s])))
                  })
                  test <- density(simuphi, from = min(simuphi), 
                    to = max(simuphi), n = length(simuphi), bw = "ucv")
                  if (test$bw < 0.1) {
                    test <- density(simuphi, from = min(simuphi), 
                      to = max(simuphi), n = length(simuphi))
                  }
                  gridf <- test$x
                  estimf <- test$y
                  estimf <- matrix(estimf, nrow = 1)
                  gridf <- matrix(gridf, nrow = 1)
                  estimf.trunc <- estimf
                  estimphi.trunc <- estimphi
                  cutoff <- FALSE
                }
            }
        }
        return(new(Class = "Freq.fit", model = model, drift.random = drift.random, 
            diffusion.random = diffusion.random, gridf = gridf, 
            mu = mu, omega = omega, a = a, lambda = lambda, cutoff = cutoff, 
            sigma2 = sigma2, index = index, estimphi = estimphi, 
            estimpsi2 = estimpsi2, estimf = estimf, estim.drift.fix = estim.drift.fix, 
            estim.diffusion.fix = estim.diffusion.fix, discrete = discrete, 
            bic = bic, aic = aic, times = times, X = X, gridg = gridg, 
            estimg = estimg))
    }
    if (estim.method == "paramMLmixture") {
        cutoff <- FALSE
        gridf <- NULL
        estimf <- NULL
        if (diffusion.random == 1) {
            warning("For the considered mixtures of SDE, there should not be random effects in the diffusion coefficient. diffusion.random is set to 0.")
            diffusion.random <- 0
        }
        if (sum(drift.random) == 0) {
            stop("There should be at least one random effect in the drift.")
        }
        if (model == "OU") {
            Mindex <- M
            index <- 1:Mindex
        }
        if (model == "CIR") {
            index <- which(rowSums(X <= 0) == 0)
            Mindex <- length(index)
            if (Mindex == 0) {
                warning("All the trajectories have non positive values the model CIR cannot be used", 
                  call. = FALSE)
                estimf <- 0
                estimphi <- 0
                estimpsi2 <- 0
                bic <- 0
                aic <- 0
                gridf <- 0
                mu <- 0
                omega <- 0
                cutoff <- 0
                sigma2 <- 0
                a <- 0
                lambda <- 0
                estimf.trunc <- 0
                estimphi.trunc <- 0
                estimpsi2.trunc <- 0
            }
        }
        U <- matrix(0, 2, Mindex)
        V <- as.list(1:Mindex)
        S <- rep(0, Mindex)
        SigDelta <- rep(0, Mindex)
        estimUV <- UVS(X[index, ], model, times)
        U <- estimUV$U
        V <- estimUV$V
        S <- estimUV$S
        SigDelta <- estimUV$SigDelta
        deter <- lapply(V, det)
        index2 <- which((deter != Inf) & (deter != 0))
        Mindex2 <- length(index2)
        V <- V[index2]
        U <- U[, index2]
        S <- S[index2]
        SigDelta <- SigDelta[index2]
        estimphi <- matrix(0, 2, Mindex2)
        for (j in 1:Mindex2) {
            estimphi[, j] <- solve(V[[j]]) \%*\% U[, j]
        }
        muinit <- matrix(0, nb.mixt, 2)
        omegainit <- matrix(0, nb.mixt, 2)
        km <- kmeans(t(estimphi), 2)
        if (sum(drift.random) == 3) {
            for (n in 1:nb.mixt) {
                label <- which(km$cluster == n)
                muinit[n, ] <- apply(estimphi[, label], 1, mean)
                omegainit[n, ] <- apply(estimphi[, label], 1, 
                  sd)
            }
        }
        if (sum(drift.random) == 2) {
            for (n in 1:nb.mixt) {
                label <- which(km$cluster == n)
                muinit[n, 2] <- mean(estimphi[2, label])
                omegainit[n, 2] <- sd(estimphi[2, label])
            }
            muinit[, 1] <- mean(estimphi[1, ])
        }
        if (sum(drift.random) == 1) {
            for (n in 1:nb.mixt) {
                label <- which(km$cluster == n)
                muinit[n, 1] <- mean(estimphi[1, label])
                omegainit[n, 1] <- sd(estimphi[1, label])
            }
            muinit[, 2] <- mean(estimphi[2, ])
        }
        probinit <- table(km$cluster)/Mindex2
        start <- list(mu = muinit, omega = omegainit, mixt.prop = probinit)
        res <- EM(U, V, S, K, drift.random, start, Niter, drift.fixed, 
            estim.drift.fix, drift.fixed.mixed, diffusion.fixed, 
            estim.diffusion.fix)
        bic <- res$BIChere
        aic <- res$AIChere
        mu <- res$mu
        omega <- res$omega
        sigma2 <- res$sigma^2
        mixt.prop <- res$mixt.prop
        probindi <- res$probindi
        if (sum(drift.random) == 3) {
            if (is.null(gridf) == 1) {
                gridf <- matrix(0, 2, 500)
                bg1 <- max(abs(estimphi[1, ])) * 1.2
                bg2 <- max(abs(estimphi[2, ])) * 1.2
                gridf[1, ] <- seq(min(abs(estimphi[1, ])) * 0.8, 
                  max(abs(estimphi[1, ])) * 1.2, length = 500)
                gridf[2, ] <- seq(min(abs(estimphi[2, ])) * 0.8, 
                  max(abs(estimphi[2, ])) * 1.2, length = 500)
            }
            if (is.null(gridf) == 0) {
                gridf <- gridf
            }
            estimf1 <- rep(0, length(gridf[1, ]), byrow = TRUE)
            estimf2 <- rep(0, length(gridf[2, ]), byrow = TRUE)
            for (n in 1:nb.mixt) {
                estimf1 <- estimf1 + mixt.prop[n] * dnorm(gridf[1, 
                  ], mean = mu[Niter, n, 1], sd = abs(omega[Niter, 
                  n, 1]))
                estimf2 <- estimf2 + mixt.prop[n] * dnorm(gridf[2, 
                  ], mean = mu[Niter, n, 2], sd = abs(omega[Niter, 
                  n, 2]))
            }
            estimf <- rbind(estimf1, estimf2)
            cutoff <- FALSE
        }
        if (sum(drift.random) == 1) {
            if (is.null(gridf) == 1) {
                gridf <- matrix(0, 2, 500)
                gridf <- seq(min(abs(estimphi[1, ])) * 0.8, max(abs(estimphi[1, 
                  ])) * 1.2, length = 500)
            }
            if (is.null(gridf) == 0) {
                gridf <- gridf
            }
            estimf <- rep(0, length(gridf), byrow = TRUE)
            for (n in 1:nb.mixt) {
                estimf <- estimf + mixt.prop[n] * dnorm(gridf, 
                  mean = mu[Niter, n, 1], sd = abs(omega[Niter, 
                    n, 1]))
            }
            gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
            estimf <- matrix(estimf, 1, length(estimf), byrow = TRUE)
            cutoff <- FALSE
        }
        if (sum(drift.random) == 2) {
            if (is.null(gridf) == 1) {
                gridf <- matrix(0, 2, 500)
                gridf <- seq(min(abs(estimphi[2, ])) * 0.8, max(abs(estimphi[2, 
                  ])) * 1.2, length = 500)
            }
            if (is.null(gridf) == 0) {
                gridf <- gridf
            }
            estimf <- rep(0, length(gridf), byrow = TRUE)
            for (n in 1:nb.mixt) {
                estimf <- estimf + mixt.prop[n] * dnorm(gridf, 
                  mean = mu[Niter, n, 2], sd = abs(omega[Niter, 
                    n, 2]))
            }
            gridf <- matrix(gridf, 1, length(gridf), byrow = TRUE)
            estimf <- matrix(estimf, 1, length(estimf), byrow = TRUE)
            cutoff <- FALSE
        }
        return(new(Class = "Freq.mixture.fit", model = model, 
            drift.random = drift.random, gridf = gridf, mu = mu, 
            omega = omega, mixt.prop = mixt.prop, cutoff = cutoff, 
            sigma2 = sigma2, index = index, probindi = probindi, 
            estimf = estimf, estimphi = estimphi, bic = bic, 
            aic = aic, estim.drift.fix = estim.drift.fix, times = times, 
            X = X))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
