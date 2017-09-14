\name{neuronal.data}
\alias{neuronal.data}
\title{
  Trajectories Interspike Of A Single Neuron Of A Ginea Pig
}
\description{
  The \code{neuronal.data} data has 240 measurements of the membrane potential in volts for one single neuron of a pig between the spikes, along time, with 2000 points for each. The step time is \eqn{delta= 0.00015} s. 
}
\usage{
  neuronal.data
}
\format{
  This data frame has a list form of length 2. The first element in the matrix named \code{Xreal}. Each row is a trajectory, that one can model by a diffusion process with random effect. The realisation can be assumed independent. The second element is a vector of times of observations \code{times}
}

\examples{

model <- "OU"

M <- 240     # number of trajectories, number of rows of the matrix of the data
T <- 0.3     # width of the interval of observation 
delta <- 0.00015   # step time
N <- T/delta  # number of points in the time interval 2000

data(neuronal.data)
# reduction of data for example to save running times
ind <- seq(1, 2000, by = 20)
X <- neuronal.data[[1]][1:100, ind]
times <- neuronal.data[[2]][ind]

  # plot(times, X[10, ], type = 'l', xlab = 'time', ylab = '', col = 'blue', ylim = c(0,0.016))
  
  
# - 1) Ornstein-Uhlenbeck with two random effects in the drift and one fixed effect in the diffusion

estim<- msde.fit(times=times, X=X, model="OU")
# summary(estim)

# - 2) Cox-Ingersoll-Ross with one random effect in the drift and one random effect in the diffusion

estim<- msde.fit(times=times, X=X, model="CIR", drift.random=1, diffusion.random=1)
# summary(estim)

# - 3) Cox-Ingersoll-Ross with one random effect in the drift and one fixed effect in the diffusion  

estim<- msde.fit(times=times, X=X, model="CIR", drift.random=1)
# summary(estim)

\donttest{
# - 4) Ornstein-Uhlenbeck with a mixture distribution for the two random effects in the drift and one fixed effect in the diffusion

estim<- msde.fit(times=times, X=X, model="OU", nb.mixt=2)
summary(estim)
}

}

\source{
  The parameters of the stochastic leaky integrate-and-fire neuronal model.
  Lansky, P., Sanda, P. and He, J. (2006).
  \emph{Journal of Computational Neuroscience} Vol 21, \bold{211--223}
}
\references{
  The parameters of the stochastic leaky integrate-and-fire neuronal model.
  Lansky, P., Sanda, P. and He, J. (2006).
  \emph{Journal of Computational Neuroscience} Vol 21, \bold{211--223}
}





\keyword{data}