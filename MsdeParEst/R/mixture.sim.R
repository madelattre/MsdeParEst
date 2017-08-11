# MsdeParEst R package ; file mixture.sim.r (last modified: 2017-08-11)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017
# UMR 518 AgroParisTech/INRA, 16 rue Claude Bernard, 75 231 Paris Cedex 05

mixture.sim <- function(M, param, prob) {
    
    N <- length(prob)
    
    index <- sample(1:N, M, replace = TRUE, prob)
    
  
        if (dim(param)[2] == 4) { #two series of random variables are simulated
            Y <- matrix(NA, dim(param)[2]/2, M)
            for (j in 1:M) {
                for (n in 1:2) {
                  Y[n, j] <- rnorm(1, param[index[j], 2 * n - 1], param[index[j], 2 * n])
                }
            }
        }
        if (dim(param)[2] == 2) { #one series of random variables is simulated
            Y <- rep(NA, M)
            for (j in 1:M) {
                Y[j] <- rnorm(1, param[index[j],1], param[index[j],2])
            }
        }
    
    res = list(Y=Y,index=index)
    
    return(res)
}


