\name{Freq.fit-class}
\Rdversion{1.1}
\docType{class}
\alias{Freq.fit-class}
\alias{plot,Freq.fit-method}
\alias{pred,Freq.fit-method}
\alias{print,Freq.fit-method}
\alias{summary,Freq.fit-method}
\alias{valid,Freq.fit-method}

\title{Class \code{"Freq.fit"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Freq.fit", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{model}:}{Object of class \code{"character"} ~~ }
    \item{\code{drift.random}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{diffusion.random}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{gridf}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{gridg}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{mu}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{omega}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{a}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lambda}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cutoff}:}{Object of class \code{"logical"} ~~ }
    \item{\code{sigma2}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{index}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{estimphi}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{estimpsi2}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{estimf}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{estimg}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{estim.drift.fix}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{estim.diffusion.fix}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{discrete}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{bic}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{aic}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{times}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{X}:}{Object of class \code{"matrix"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "Freq.fit")}: ... }
    \item{pred}{\code{signature(x = "Freq.fit")}: ... }
    \item{print}{\code{signature(x = "Freq.fit")}: ... }
    \item{summary}{\code{signature(object = "Freq.fit")}: ... }
    \item{valid}{\code{signature(x = "Freq.fit")}: ... }
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
showClass("Freq.fit")
}
\keyword{classes}
