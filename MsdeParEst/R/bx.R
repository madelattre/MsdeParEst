# MsdeParEst R package ; bx.r (last modified: 2017-08-28)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017

bx <- function(x) {
        return(matrix(c(1, -x), 2, 1))
}
