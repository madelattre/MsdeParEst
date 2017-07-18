\name{likelihoodNormalGamma}
\alias{likelihoodNormalGamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
likelihoodNormalGamma(a, lambda, mu, omega, U, V, S, SigDelta, K, drift.random)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{omega}{
%%     ~~Describe \code{omega} here~~
}
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
function (a, lambda, mu, omega, U, V, S, SigDelta, K, drift.random) 
{
    if (is.infinite(suppressWarnings(log(gamma(a + K/2))))) {
        L <- Inf
    }
    else {
        invlambda <- 1/lambda
        M <- length(S)
        loglik <- vector(length = M)
        V1 <- vector(length = M)
        V2 <- vector(length = M)
        for (j in 1:M) {
            A <- (I2 + V[[j]] \%*\% Omega)
            Rinv <- solve(A) \%*\% V[[j]]
            b <- mu - estimphi[, j]
            V1[j] <- a * log(invlambda) + log(gamma(a + K/2)) - 
                log(gamma(a)) - (a + K/2) * log(invlambda + ST[j]/2)
            V2[j] <- -log(det(A))/2 - K/(2 * S[j]) * (t(b) \%*\% 
                Rinv \%*\% b - t(U[, j]) \%*\% estimphi[, j])
            loglik[j] <- K * log(2 * pi) + sum(SigDelta[j]) - 
                2 * V1[j] - 2 * V2[j]
        }
        L <- sum(loglik)
    }
    return(L)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
