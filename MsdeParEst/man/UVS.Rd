\name{UVS}
\alias{UVS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
UVS(X, model, times)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
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
function (X, model, times) 
{
    M <- dim(X)[1]
    K <- dim(X)[2]
    Xm <- X[, -K]
    delta <- diff(times)
    Tend <- times[K]
    S <- rep(NA, M)
    if (model == "OU") {
        for (j in 1:M) {
            S[j] <- sum((X[j, 2:K] - X[j, 1:(K - 1)])^2/delta)
        }
    }
    if (model == "CIR") {
        for (j in 1:M) {
            S[j] <- sum((X[j, 2:K] - X[j, 1:(K - 1)])^2/(delta * 
                X[j, 1:(K - 1)]))
        }
    }
    U <- matrix(0, 2, M)
    V <- as.list(1:M)
    b <- as.list(1:M)
    Int1 <- rowSums(Xm * matrix(delta, M, length(delta)))
    if (model == "OU") {
        SigDelta <- rep(sum(log(delta)), M)
        Int2 <- rowSums(Xm^2 * matrix(delta, M, length(delta)))
        for (j in 1:M) {
            b[[j]] <- matrix(apply(matrix(X[j, ], 1, K), 2, bx, 
                0, c(1, 2)), 2, K)
            U[, j] <- rowSums((b[[j]][, 1:(K - 1)] * matrix((X[j, 
                2:K] - X[j, 1:(K - 1)]), 2, K - 1, byrow = TRUE)))
            V[[j]] <- matrix(c(Tend, -Int1[j], -Int1[j], Int2[j]), 
                2, 2)
        }
    }
    if (model == "CIR") {
        SigDelta <- rep(0, M)
        Int3 <- rowSums(1/Xm * matrix(delta, M, length(delta)))
        b <- as.list(1:M)
        bsig <- as.list(1:M)
        for (j in 1:M) {
            b[[j]] <- matrix(apply(matrix(X[j, ], 1, K), 2, bx, 
                0, c(1, 2)), 2, K)
            bsig[[j]] <- matrix(-1, 2, K)
            bsig[[j]][1, ] <- 1/X[j, ]
            U[, j] = rowSums((bsig[[j]][, 1:(K - 1)] * matrix((X[j, 
                2:K] - X[j, 1:(K - 1)]), 2, K - 1, byrow = TRUE)))
            V[[j]] <- matrix(c(Int3[j], -Tend, -Tend, Int1[j]), 
                2, 2)
            SigDelta[j] <- sum(log(delta) + log(X[j, 2:K]))
        }
    }
    return(list(U = U, V = V, S = S, SigDelta = SigDelta))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
