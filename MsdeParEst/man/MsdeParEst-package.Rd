\name{MsdeParEst-package}
\alias{MsdeParEst-package}
\docType{package}
\title{Parametric estimation in mixed-effects stochastic differential equations}
\description{
Parametric estimation in mixed-effects stochastic differential equations
}
\details{
This package is dedicated to parametric estimation in the following mixed-effects stochastic differential equations:
\deqn{dX_j(t)= (\alpha_j- \beta_j X_j(t))dt + \sigma_j \ a(X_j(t)) dW_j(t),}
\eqn{j=1,\ldots,M}, where the \eqn{(W_j(t))} are independant Wiener processes and the \eqn{(X_j(t))} are observed without noise. The volatility function \eqn{a(x)} is known and can be either \eqn{a(x)=1} (Ornstein-Uhlenbeck process) or \eqn{a(x)=\sqrt{x}} (Cox-Ingersoll-Ross process).

Different estimation methods are implemented depending on whether there are random effects in the drift and/or in the diffusion coefficient:
\enumerate{
\item The diffusion coefficient is fixed \eqn{\sigma_j \equiv \sigma} and the parameters in the drift are Gaussian random variables: 
\enumerate{
\item either \eqn{\alpha_j \equiv \alpha} and \eqn{\beta_j \sim N(\mu,\Omega), j=1,\ldots,M},
\item or \eqn{\beta_j \equiv \beta} and \eqn{\alpha_j \sim N(\mu,\Omega), j=1,\ldots,M},
\item or \eqn{(\alpha_j,\beta_j) \sim N(\mu,\Omega), j=1,\ldots,M}. 
}
\eqn{\mu}, \eqn{\Omega} and potentially the fixed effects \eqn{\sigma}, \eqn{\alpha}, \eqn{\beta} are estimated as proposed in [1] and [4]. The extension to mixtures of Gaussian distributions is also implemented by following [3].

\item The coefficients in the drift are fixed \eqn{\alpha_j \equiv \alpha} and \eqn{\beta_j \equiv \beta} and the diffusion coefficient \eqn{1/\sigma_j^2} follows a Gamma distribution \eqn{1/\sigma_j^2 \sim \Gamma(a,\lambda), j=1,\ldots,M}. \eqn{a}, \eqn{\lambda}, and potentially the fixed effects \eqn{\alpha} and \eqn{\beta} are estimated by the method published in [2].

\item There are random effects in the drift and in the diffusion, such that \eqn{1/\sigma_j^2 \sim \Gamma(a,\lambda)} and 
\enumerate{
\item either \eqn{\alpha_j \equiv \alpha} and \eqn{\beta_j|\sigma_j \sim N(\mu,\sigma_j^2 \Omega)},
\item or \eqn{\beta_j \equiv \beta} and \eqn{\alpha_j|\sigma_j \sim N(\mu,\sigma_j^2 \Omega)},
\item or \eqn{(\alpha_j,\beta_j)|\sigma_j \sim N(\mu,\sigma_j^2 \Omega)}.
}
\eqn{a}, \eqn{\lambda}, \eqn{\mu}, \eqn{\Omega} and potentially the fixed effects \eqn{\alpha} and \eqn{\beta} are estimated by following [5].

}

Some examples and more details on how to use the different functions of the package are available on the package's vignette by typing \bold{vignette("MsdeParEst")} in the R console. 
}

\references{
\bold{[1]} Maximum Likelihood Estimation for Stochastic Differential Equations with Random Effects, Delattre, M., Genon-Catalot, V. and Samson, A. \emph{Scandinavian Journal of Statistics 40(2) 2012} \bold{322-343}. 

\bold{[2]} Estimation of population parameters in stochastic differential equations with random effects in the diffusion coefficient, Delattre, M., Genon-Catalot, V. and Samson, A. \emph{ESAIM:PS 19 2015} \bold{671-688}.

\bold{[3]} Mixtures of stochastic differential equations with random effects: application to data clustering, Delattre, M., Genon-Catalot, V. and Samson, A. \emph{Journal of Statistical Planning and Inference 173 2016} \bold{109-124}. 

\bold{[4]} Parametric inference for discrete observations of diffusion processes with mixed effects, Delattre, M., Genon-Catalot, V. and Laredo, C. \emph{Stochastic Processes and their Applications}. To appear. (Available pre-publication hal-0133263, 2016). 

\bold{[5]} Estimation of the joint distribution of random effects for a discretely observed diffusion with random effects, Delattre, M., Genon-Catalot, V. and Laredo, C. \emph{hal-01446063 2017}.

}
