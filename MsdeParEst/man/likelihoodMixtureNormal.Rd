\name{likelihoodMixtureNormal}
\alias{likelihoodMixtureNormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
likelihoodMixtureNormal(mu, omega, sigma, mixt.prop, U, V, S, K, estimphi, drift.random)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{omega}{
%%     ~~Describe \code{omega} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{mixt.prop}{
%%     ~~Describe \code{mixt.prop} here~~
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
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{estimphi}{
%%     ~~Describe \code{estimphi} here~~
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
function (mu, omega, sigma, mixt.prop, U, V, S, K, estimphi, 
    drift.random) 
{
    M <- dim(U)[2]
    N <- length(mixt.prop)
    probRep <- matrix(mixt.prop, M, N, byrow = TRUE)
    Lindicomponent <- matrix(NA, M, N)
    for (i in 1:M) {
        Lindicomponent[i, ] <- sapply(1:N, function(n) {
            exp(-1/2 * likelihoodNormalindi(mu[n, ], omega[n, 
                ], sigma, U[, i], V[[i]], S[i], K, estimphi[, 
                i], drift.random))
        })
    }
    probindi <- exp(log(probRep) + log(Lindicomponent))
    loglik <- -2 * sum(log(rowSums(probindi)))
    return(loglik)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
