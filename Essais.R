source('bx.R')
source('contrastNormal.R')
source('EM.R')
source('EstParamGamma.R')
source('EstParamNormal.R')
source('EstParamNormalGamma.R')
source('likelihoodGamma.R')
source('likelihoodNormal.R')
source('likelihoodNormalGamma.R')
source('mixture.sim.R')
source('msde.fit.R')
source('msde.sim.R')
source('simu.randomvariable.R')
source('UVS.R')


# Example 1: one random effect in the drift and one random effect in the diffusion coefficient.
# -- Simulation 
M <- 100
Tmax <- 5
N <- 5000
model <- "OU"
drift.random <- 2
diffusion.random <- 1
drift.fixed <- 0
drift.param <- c(0.5,0.5)
diffusion.param <- c(8,1/2)

sim1 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random, 
                  diffusion.random = diffusion.random, drift.fixed = drift.fixed, 
                  mixture = 0, drift.param = drift.param, diffusion.param = diffusion.param)
                  
# -- Estimation

# -----Fixed effect in the drift estimated
res1 <- msde.fit(times = sim1$times, X = sim1$X, model = "OU", drift.random = 2,
                 diffusion.random = 1, estim.drift.fix = 1, mixture = 0)
summary(res1)
#' print(res1)
#' #pred1 <- pred(res1, invariant = 0, level = 0.05, newwindow = FALSE, plot.pred = TRUE)
#' valid(res1)
#' plot(res1)
#' 
# ----- Fixed effect in the drift known and not estimated
res1bis <- msde.fit(times = sim1$times, X = sim1$X, model = "OU", drift.random = 2,
                    diffusion.random = 1, drift.fixed=0, mixture = 0)
summary(res1bis)
#' 
# Example 2: one random effect in the drift and one fixed effect in the diffusion coefficient

# -- Simulation
M <- 100
Tmax <- 5
N <- 5000
model <- "OU"
diffusion.random <- 0
diffusion.param <- 0.5
drift.random <- 2
drift.fixed <- 10
drift.param <- c(1,sqrt(0.4/4))

sim2 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random,
                 diffusion.random = diffusion.random, drift.fixed = drift.fixed,
                 mixture=0, drift.param = drift.param,
                 diffusion.param = diffusion.param)

# -- Estimation
res2 <- msde.fit(times = sim2$times, X = sim2$X, model = "OU", drift.random = 2,
                 diffusion.random = 0, estim.drift.fix = 1, mixture = 0)

summary(res2)
plot(res2)

# Example 3: two random effects in the drift and one random effect in the diffusion coefficient

# -- Simulation
M <- 100
Tmax <- 5
N <- 5000
model <- "OU"
drift.random <- c(1,2)
diffusion.random <- 1
drift.param <- c(1,0.5,0.5,0.5)
diffusion.param <- c(8,1/2)

sim3 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random,
                 diffusion.random = diffusion.random, mixture = 0,
                 drift.param = drift.param, diffusion.param = diffusion.param)

# -- Estimation

res3 <- msde.fit(times = sim3$times, X = sim3$X, model = "OU", drift.random = c(1,2),
                 diffusion.random = 1, mixture = 0)
summary(res3)
plot(res3)

# Example 4: fixed effects in the drift and one random effect in the diffusion coefficient

# -- Simulation
M <- 100
Tmax <- 5
N <- 5000

model <- "OU"
drift.random <- 0
diffusion.random <- 1
drift.fixed <- c(0,1)
diffusion.param <- c(5,3)
sim4 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random,
                 diffusion.random = diffusion.random, drift.fixed = drift.fixed,
                 diffusion.param = diffusion.param)

# -- Estimation
res4 <- msde.fit(times = sim4$times, X = sim4$X, model = "OU", drift.random = 0,
                 diffusion.random = 1, mixture = 0, estim.drift.fix = 0,
                 drift.fixed = c(0,1), discrete = 1)

summary(res4)

# Example 5: one fixed effect and one mixture random effect in the drift, and one fixed effect in
# the diffusion coefficient

# -- Simulation
M <- 100
Tmax <- 5
N <- 5000
diffusion.random <- 0
diffusion.param <- 0.1
model <- "OU"
drift.random <- 1
drift.fixed <- 1
density.phi <- 'mixture.normal'
nb.mixt <- 2
mixt.prop <- c(0.5,0.5)
param.ea1 <- c(0.5, 0.25, 1.8, 0.25)
param.ea2 <- c(1, 0.25, 1, 0.25)
drift.param <- param.ea1

sim5 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random,
                 diffusion.random = diffusion.random, drift.fixed = drift.fixed,
                 mixture=1, drift.param = drift.param,
                 diffusion.param = diffusion.param, nb.mixt = nb.mixt, mixt.prop = mixt.prop)

# -- Estimation
res5 <- msde.fit(times = sim5$times, X = sim5$X, model = "OU", drift.random = 1,
                 estim.drift.fix = 1, diffusion.random = 0, estim.diffusion.fix = 1,
                 mixture = 1, nb.mixt=2, Niter = 25)

summary(res5)
print(res5)
plot(res5)