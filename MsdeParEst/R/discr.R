# MsdeParEst R package ; file discr.r (last modified: 2017-09-15)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017

discr <- function(x, p) {
    cp <- cumsum(p)
    U <- runif(1, 0, max(cp))
    x[which(cp >= U)[1]]
}
