rm(list = ls())

devtools::load_all(".")

#nettoie les codes
library("formatR")
formatR::tidy_dir("R")



# packages
require("MASS")
require("plot3D")
require("roxygen2")
require("devtools")
require("sde")
require("moments")
require("mvtnorm")

# Example 1: one random effect in the drift and one random effect in the diffusion coefficient.

# source('msde.pred.R')
# source('msde.sim.R')
# source('msde.fit.R')
# source('UVS.R')
# source('bx.R')
# source('contrastNormal.R')
# source('EM.R')
# source('EstParamGamma.R')
# source('EstParamNormal.R')
# source('EstParamNormalGamma.R')
# source('likelihoodGamma.R')
# source('likelihoodNormal.R')
# source('likelihoodNormalGamma.R')
# source('mixture.sim.R')
# source('disc.R')

# Example 1: One random effect in the drift and one random effect in the diffusion

# -- Simulation

sim1 <- msde.sim(M = 100, T = 5, N = 5000, model = 'OU', 
                 drift.random = 2, drift.param = c(0,0.5,0.5), 
                 diffusion.random = 1, diffusion.param = c(8,1/2))



# -- Fixed effect in the drift estimated

# -- Estimation 
res1 <- msde.fit(times = sim1$times, X = sim1$X, model = 'OU', drift.random = 2, diffusion.random = 1)

# -- Estimation and prediction

res1.valid <- msde.fit(times = sim1$times, X = sim1$X, model = 'OU', drift.random = 2, diffusion.random = 1,valid=1)

summary(res1.valid) # PB ESTIMATION DE L'ECART TYPE: effet du tirage
plot(res1.valid)

# ----- Fixed effect in the drift known and not estimated

# -- Estimation
res1bis <- msde.fit(times = sim1$times, X = sim1$X, model = "OU", drift.random = 2, 
                    diffusion.random = 1, drift.fixed = 0)
summary(res1bis)

# -- Estimation and prediction
res1bis.valid <- msde.fit(times = sim1$times, X = sim1$X, model = "OU", drift.random = 2, 
                          diffusion.random = 1, drift.fixed = 0, valid = 1)
summary(res1bis.valid)

# Example 2: one random effect in the drift and one fixed effect in the diffusion coefficient

# -- Simulation


sim2 <- msde.sim(M = 100, T = 5, N = 5000, model = 'OU', drift.random = 2,
                 diffusion.random = 0, drift.param = c(0,1,sqrt(0.4/4)), diffusion.param = 0.5)

# -- Estimation without prediction
res2 <- msde.fit(times = sim2$times, X = sim2$X, model = 'OU', drift.random = 2,
                 diffusion.random = 0)

summary(res2)
plot(res2)

# -- Estimation with prediction
res2.valid <- msde.fit(times = sim2$times, X = sim2$X, model = 'OU', drift.random = 2,
                       diffusion.random = 0, valid = 1)

summary(res2.valid)
plot(res2.valid)


# Example 3: two random effects in the drift and one random effect in the diffusion coefficient

# -- Simulation

sim3 <- msde.sim(M = 100, T = 5, N = 5000, model = 'OU', drift.random = c(1,2),
                 diffusion.random = 1, drift.param = c(1,0.5,0.5,0.5), diffusion.param = c(8,1/2))

# -- Estimation without prediction

res3 <- msde.fit(times = sim3$times, X = sim3$X, model = 'OU', drift.random = c(1,2),
                 diffusion.random = 1)
summary(res3)
plot(res3)

# -- Estimation with prediction

res3.valid <- msde.fit(times = sim3$times, X = sim3$X, model = 'OU', drift.random = c(1,2),
                 diffusion.random = 1, valid = 1)
summary(res3.valid)
plot(res3.valid)

# Example 4: fixed effects in the drift and one random effect in the diffusion coefficient

# -- Simulation

sim4 <- msde.sim(M = 100, T = 5, N = 5000, model = 'OU', drift.random = 0,
                 diffusion.random = 1, drift.param = c(0,1), diffusion.param = c(5,3))

# -- Estimation
res4 <- msde.fit(times = sim4$times, X = sim4$X, model = 'OU', drift.random = 0,
                 diffusion.random = 1, drift.fixed = c(0,0))
#, drift.fixed = c(0,0)

summary(res4)
plot(res4)

res4.valid <- msde.fit(times = sim4$times, X = sim4$X, model = 'OU', drift.random = 0,
                 diffusion.random = 1, drift.fixed = c(0,0), valid = 1)
#, drift.fixed = c(0,0)

summary(res4.valid)
plot(res4.valid)

# Example 5: one fixed effect and one mixture random effect in the drift, and one fixed effect in
# the diffusion coefficient

# -- Simulation


sim5 <- msde.sim(M = 100, T = 5, N = 5000, model = 'OU', 
                 drift.random = 1, drift.param = matrix(c(0.5,1.8,0.25,0.25,1,1),nrow=2,byrow=F),
                 diffusion.random = 0, diffusion.param = 0.1, 
                 nb.mixt = 2, mixt.prop = c(0.5,0.5))

# -- Estimation
res5 <- msde.fit(times = sim5$times, X = sim5$X, model = 'OU', drift.random = 1,
                 nb.mixt=2, Niter = 10)

summary(res5)
plot(res5)


res5.valid <- msde.fit(times = sim5$times, X = sim5$X, model = 'OU', drift.random = 1,
                 nb.mixt=2, Niter = 10, valid = 1)

summary(res5.valid)

# Example 6: one mixture of two random effects in the drift, and one fixed effect in
# the diffusion coefficient

sim6 <- msde.sim(M = 100, T = 5, N = 5000, model = 'OU', drift.random = c(1,2),
                 diffusion.random = 0, drift.param = matrix(c(0.5,1.8,0.25,0.25,1,2,0.25,0.25),nrow=2,byrow=F), 
                 diffusion.param = 0.1, nb.mixt = 2, mixt.prop = c(0.5,0.5))

# -- Estimation without prediction
res6 <- msde.fit(times = sim6$times, X = sim6$X, model = 'OU', drift.random = c(1,2),
                 nb.mixt=2, Niter = 10)

summary(res6)
plot(res6)

# -- Estimation with prediction
res6.valid <- msde.fit(times = sim6$times, X = sim6$X, model = 'OU', drift.random = c(1,2),
                 nb.mixt=2, Niter = 10, valid = 1)

summary(res6.valid)
plot(res6.valid)

# Example 7: CIR with one random effect in the drift and one random effect in the diffusion 
# coefficient

sim7 <- msde.sim(M = 100, T = 5, N = 5000, model = 'CIR', drift.random = 2,
                         diffusion.random = 1, drift.param = c(4,1,0.1), diffusion.param = c(8,0.5),
                 X0 = 1)

res7 <- msde.fit(times = sim7$times, X = sim7$X, model = 'CIR', drift.random = 2,
                       diffusion.random = 1)

summary(res7)


############

 sim <- msde.sim(M = 30, T = 5, N = 5000, model = 'OU', drift.random = 2,
                diffusion.random = 0, drift.param = c(0,1,sqrt(0.4/4)), diffusion.param = 0.5)

 # Example : two random effects in the drift and one random effect in the diffusion coefficient

 sim <- msde.sim(M = 30, T = 5, N = 5000, model = 'OU', drift.random = c(1,2),
                diffusion.random = 1, drift.param = c(1,0.5,0.5,0.5), diffusion.param = c(8,1/2))

 # Example 2 : one fixed effect and one mixture random effect in the drift, and one fixed effect in
 # the diffusion coefficient

 sim <- msde.sim(M = 30, T = 5, N = 5000, model = 'OU',
                 drift.random = 1, drift.param = matrix(c(0.5,1.8,0.25,0.25,1,1),nrow=2,byrow=F),
                 diffusion.random = 0, diffusion.param = 0.1,
                 nb.mixt = 2, mixt.prop = c(0.5,0.5))

# Example 3 : CIR with one random effect in the drift and one random effect in the diffusion
# coefficient

sim <- msde.sim(M = 30, T = 5, N = 5000, model = 'CIR', drift.random = 2,
                diffusion.random = 1, drift.param = c(4,1,0.1), diffusion.param = c(8,0.5),
                X0 = 1)


##############

sim <- msde.sim(M = 50, T = 3, N = 3000, model = 'OU',
                drift.random = 2, drift.param = c(0,0.5,0.5),
                diffusion.random = 1, diffusion.param = c(8,1/2))

# -- Fixed effect in the drift estimated
res <- msde.fit(times = sim$times, X = sim$X, model = 'OU', drift.random = 2,
diffusion.random = 1)

# ----- Fixed effect in the drift known and not estimated
resbis <- msde.fit(times = sim$times, X = sim$X, model = "OU", drift.random = 2,
                    diffusion.random = 1, drift.fixed = 0)
summary(resbis)

#######

sim <- msde.sim(M = 100, T = 5, N = 5000, model = 'OU', drift.random = c(1,2),
                                diffusion.random = 0,
                                drift.param = matrix(c(0.5,1.8,0.25,0.25,1,2,0.25,0.25),nrow=2,byrow=FALSE),
                                diffusion.param = 0.1, nb.mixt = 2, mixt.prop = c(0.5,0.5))

                # -- Estimation without validation
                res <- msde.fit(times = sim$times, X = sim$X, model = 'OU', drift.random = c(1,2),
                                nb.mixt=2, Niter = 10)