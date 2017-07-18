\name{EstParamGamma}
\alias{EstParamGamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EstParamGamma(U, V, S, SigDelta, K, drift.param = c(0, 0), drift.estim = 1)
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
  \item{drift.param}{
%%     ~~Describe \code{drift.param} here~~
}
  \item{drift.estim}{
%%     ~~Describe \code{drift.estim} here~~
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
function (U, V, S, SigDelta, K, drift.param = c(0, 0), drift.estim = 1) 
{
    M <- length(S)
    k <- 0.1
    estimGamma <- (K/S) * ((S/K) >= (k/sqrt(K)))
    init.a <- (mean(estimGamma))^2/var(estimGamma)
    init.lambda <- var(estimGamma)/mean(estimGamma)
    if (drift.estim == 0) {
        ln = function(param) {
            contrastGamma(exp(param[1]), exp(param[2]), U, V, 
                S, K, drift.param)
        }
        res <- optim(c(log(init.a), log(init.lambda)), fn = ln, 
            method = "Nelder-Mead")
        a <- exp(res$par[1])
        lambda <- exp(res$par[2])
        mu <- drift.param
        BIChere <- likelihoodGamma(a, lambda, U, V, S, SigDelta, 
            K, mu) + 2 * log(M)
        AIChere <- likelihoodGamma(a, lambda, U, V, S, SigDelta, 
            K, mu) + 4
    }
    if (drift.estim == 1) {
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
            contrastGamma(exp(param[1]), exp(param[2]), U, V, 
                S, K, c(param[3], param[4]))
        }
        res <- optim(c(log(init.a), log(init.lambda), init.phi), 
            fn = ln, method = "Nelder-Mead")
        a <- exp(res$par[1])
        lambda <- exp(res$par[2])
        mu <- c(res$par[3], res$par[4])
        BIChere <- likelihoodGamma(a, lambda, U, V, S, SigDelta, 
            K, mu) + 4 * log(M)
        AIChere <- likelihoodGamma(a, lambda, U, V, S, SigDelta, 
            K, mu) + 4 * 2
    }
    return(list(mu = mu, a = a, lambda = lambda, BIChere = BIChere, 
        AIChere = AIChere))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
