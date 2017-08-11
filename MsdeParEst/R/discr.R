# MsdeParEst R package ; file discr.r (last modified: 2017-08-11)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017
# UMR 518 AgroParisTech/INRA, 16 rue Claude Bernard, 75 231 Paris Cedex 05



discr <- function(x, p) {
    cp <- cumsum(p)
    U <- runif(1, 0, max(cp))
    x[which(cp >= U)[1]]
}
