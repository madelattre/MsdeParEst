\name{Freq.mixture.fit-class}
\Rdversion{1.1}
\docType{class}
\alias{Freq.mixture.fit-class}
\alias{plot,Freq.mixture.fit-method}
\alias{print,Freq.mixture.fit-method}
\alias{summary,Freq.mixture.fit-method}

\title{Class \code{"Freq.mixture.fit"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Freq.mixture.fit", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{model}:}{Object of class \code{"character"} ~~ }
    \item{\code{drift.random}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{gridf}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{mu}:}{Object of class \code{"array"} ~~ }
    \item{\code{omega}:}{Object of class \code{"array"} ~~ }
    \item{\code{mixt.prop}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{cutoff}:}{Object of class \code{"logical"} ~~ }
    \item{\code{sigma2}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{index}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{estimphi}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{probindi}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{estimf}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{estim.drift.fix}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{bic}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{aic}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{times}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{X}:}{Object of class \code{"matrix"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "Freq.mixture.fit")}: ... }
    \item{print}{\code{signature(x = "Freq.mixture.fit")}: ... }
    \item{summary}{\code{signature(object = "Freq.mixture.fit")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("Freq.mixture.fit")
}
\keyword{classes}
