\name{likelihoodGamma}
\alias{likelihoodGamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
likelihoodGamma(a, lambda, U, V, S, SigDelta, K, drift.fixed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
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
  \item{drift.fixed}{
%%     ~~Describe \code{drift.fixed} here~~
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
function (a, lambda, U, V, S, SigDelta, K, drift.fixed) 
{
    if (is.infinite(log(gamma(a + K/2)))) {
        L <- Inf
    }
    else {
        invlambda <- 1/lambda
        M <- length(S)
        loglik <- vector(length = M)
        phi1 <- drift.fixed[1]
        phi2 <- drift.fixed[2]
        U1 <- rep(0, M)
        U2 <- rep(0, M)
        V11 <- rep(0, M)
        V22 <- rep(0, M)
        V12 <- rep(0, M)
        ST <- rep(0, M)
        for (j in 1:M) {
            U1[j] <- U[1, j]
            U2[j] <- U[2, j]
            V11[j] <- V[[j]][1, 1]
            V22[j] <- V[[j]][2, 2]
            V12[j] <- V[[j]][1, 2]
        }
        for (j in 1:M) {
            ST[j] <- (S[j] - 2 * (phi1 * U1[j] + phi2 * U2[j]) + 
                phi1^2 * V11[j] + phi2^2 * V22[j] + 2 * phi1 * 
                phi2 * V12[j])
        }
        for (j in 1:M) {
            loglik[j] <- K * log(2 * pi) + sum(SigDelta[j]) - 
                2 * (a * log(invlambda) + log(gamma(a + K/2)) - 
                  log(gamma(a)) + (a + K/2) * log(invlambda + 
                  ST[j]/2))
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
