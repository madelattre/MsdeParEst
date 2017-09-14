rm(list = ls())

#setwd("~/Desktop/MsdeParEst/MsdeParEst")
setwd("~/Montages_reseau/home_maiage/Recherche/MsdeParEst/MsdeParEst/MsdeParEst/R")
devtools::load_all(".")


#nettoie les codes
library("formatR")
formatR::tidy_dir("R")



# packages
#library("invgamma")

require("MASS")
require("plot3D")
require("roxygen2")
require("devtools")
require("sde")
require("moments")
require("mvtnorm")

#pour générer le doc
document()

? mixture.sim
? msde.sim
? msde.fit
? msde.pred
? EM
? UVS
? EstParamGamma
? EstParamNormal
? EstParamNormalGamma
? likelihoodNormal
? likelihoodNormalGamma
? likelihoodGamma
? contrastNormal 
? bx
? discr # changer 

devtools::load_all(".")


#################################################
# EXAMPLE 1
################################################

#------- Simulation

M <- 100
Tmax <- 5
N <- 1000
model <- "OU"
drift.random <- 1
diffusion.random <- 0

drift.fixed <- 1
density.phi <- 'normal'
drift.param <- c(0.5,0.5)
diffusion.param <- 1    

sim1 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random,
                  diffusion.random = diffusion.random, diffusion.param = diffusion.param , drift.fixed = drift.fixed,
                 drift.param = drift.param)

# -- Estimation

# -----Fixed effect in the drift estimated
res1 <- msde.fit(times = sim1$times, X = sim1$X, model = "OU", drift.random = 1, diffusion.random = 0, diffusion.fixed = 1,
                  estim.drift.fix = 1)
# ne va pas package ‘invgamma’ is not available (for R version 3.0.2)

summary(res1)

print(res1)

test <- valid(res1, plot.valid = 1)
#plot(test)

plot(res1)


pred1 <- msde.pred( times = sim1$times, X = sim1$X, model = model, 
                    drift.random = drift.random,
               diffusion.random = diffusion.random, 
              level = 0.05)

# msde.pred <- function(times, X, model = c("OU", "CIR"), drift.random, drift.fixed = NULL, estim.drift.fix = 0, 
#                       diffusion.random = 0, diffusion.fixed = NULL, estim.diffusion.fix = 0, mixture = 0, nb.mixt = 1, drift.fixed.mixed = 0, 
#                       Niter = 10, discrete = 1, plot.pred = TRUE, level = 0.05, newwindow = FALSE)

##########################################################################
# Mixture 2
##########################################################################

M <- 100
Tmax <- 5
N <- 1000
diffusion.random <- 0
diffusion.fixed <- 0.1
model <- 'OU'
drift.random <- 1
drift.fixed <- 1
nb.mixt <- 2
mixture = 1
mixt.prop <- c(0.5,0.5)
param.ea1 <- c(0.5, 0.25, 1.8, 0.25)
#param.ea2 <- c(1, 0.25, 1, 0.25)
drift.param <- c(0.5, 0.25, 1.8, 0.25)
Niter = 10
drift.fixed.mixed = 0
discrete = 1
estim.diffusion.fix = 1
estim.drift.fix = 0

sim5 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random,
                 diffusion.random = diffusion.random, drift.fixed = drift.fixed,
                 mixture=1, drift.param = drift.param,
                 diffusion.param = diffusion.fixed, nb.mixt = nb.mixt, mixt.prop = mixt.prop)

# -- Estimation
res5 <- msde.fit(times = sim5$times, X = sim5$X, model = model, drift.random = drift.random,
                 estim.drift.fix = estim.drift.fix, diffusion.random = diffusion.random, estim.diffusion.fix = 1,
                 mixture = mixture, nb.mixt=nb.mixt, Niter = Niter)


#test = out(res5)
summary(res5)
print(res5)
plot(res5)

test <- valid(res1, plot.valid = 1)
#plot(test)
############################################# PREDICTION

res5 <- msde.pred(times = sim5$times, X = sim5$X, model = model, drift.random = drift.random,
                 estim.drift.fix = estim.drift.fix, diffusion.random = diffusion.random, estim.diffusion.fix = estim.diffusion.fix,
                 mixture = mixture, nb.mixt=nb.mixt, Niter = Niter)



##########################################################################
# Mixture 3
##########################################################################
M <- 100
Tmax <- 5
N <- 5000
model <- 'OU'
drift.random <- c(1,2)
diffusion.random <- 1
density.phi <- 'normalnormal'
drift.param <- c(1,0.5,0.5,0.5)
diffusion.param <- c(8,1/2)


nb.mixt <- 2
mixture = 0
mixt.prop <- c(0.5,0.5)
param.ea1 <- c(0.5, 0.25, 1.8, 0.25)
#param.ea2 <- c(1, 0.25, 1, 0.25)
drift.param <- c(0.5, 0.25, 1.8, 0.25)
Niter = 10
drift.fixed.mixed = 0
discrete = 1
estim.diffusion.fix = 1
estim.drift.fix = 0

sim6 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random,
                 diffusion.random = diffusion.random, drift.fixed = drift.fixed,
                 mixture=1, drift.param = drift.param,
                 diffusion.param = diffusion.fixed, nb.mixt = nb.mixt, mixt.prop = mixt.prop)

# -- Estimation
res6 <- msde.fit(times = sim5$times, X = sim5$X, model = model, drift.random = drift.random,
                 estim.drift.fix = estim.drift.fix, diffusion.random = diffusion.random, estim.diffusion.fix = 1,
                 mixture = mixture, nb.mixt=nb.mixt, Niter = Niter)


summary(res6)
print(res6)
plot(res6)

test <- valid(res6, plot.valid = 1)

### -- PREDICTION

res6 <- msde.pred(times = sim5$times, X = sim5$X, model = model, drift.random = drift.random,
                  estim.drift.fix = estim.drift.fix, diffusion.random = diffusion.random, estim.diffusion.fix = estim.diffusion.fix,
                  mixture = mixture, nb.mixt=nb.mixt, Niter = Niter)






# -- Simulation
M <- 100
Tmax <- 5
N <- 5000
model <- 'OU'
drift.random <- c(1,2)
diffusion.random <- 1
density.phi <- 'normalnormal'
drift.param <- c(1,0.5,0.5,0.5)
diffusion.param <- c(8,1/2)

sim3 <- msde.sim(M = M, T = Tmax, N = N, model = model, drift.random = drift.random,
                 diffusion.random = diffusion.random, mixture = 0,
                 drift.param = drift.param, diffusion.param = diffusion.param)

