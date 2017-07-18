\name{likelihoodNormalindi}
\alias{likelihoodNormalindi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
likelihoodNormalindi(mu, omega, sigma, Uj, Vj, Sj, K, estimphij, drift.random)
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
  \item{Uj}{
%%     ~~Describe \code{Uj} here~~
}
  \item{Vj}{
%%     ~~Describe \code{Vj} here~~
}
  \item{Sj}{
%%     ~~Describe \code{Sj} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{estimphij}{
%%     ~~Describe \code{estimphij} here~~
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
function (mu, omega, sigma, Uj, Vj, Sj, K, estimphij, drift.random) 
{
    if (length(drift.random) == 2) {
        Omega <- matrix(c(omega[1]^2, 0, 0, omega[2]^2), 2, 2, 
            byrow = TRUE)
    }
    if (sum(drift.random) == 1) {
        Omega <- matrix(c(omega[1]^2, 0, 0, 0), 2, 2, byrow = TRUE)
    }
    if (sum(drift.random) == 2) {
        Omega <- matrix(c(0, 0, 0, omega[2]^2), 2, 2, byrow = TRUE)
    }
    I2 <- diag(c(1, 1))
    A <- (I2 + Vj \%*\% Omega)
    Rinv <- solve(A) \%*\% Vj
    b <- mu - estimphij
    L <- log(det(A)) + 1/sigma^2 * (t(b) \%*\% Rinv \%*\% b - t(Uj) \%*\% 
        estimphij)
    return(L)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
