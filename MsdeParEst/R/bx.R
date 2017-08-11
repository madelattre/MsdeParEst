# MsdeParEst R package ; bx.r (last modified: 2017-08-11)
# Authors: M. Delattre, C. Dion
# Copyright INRA 2017
# UMR 518 AgroParisTech/INRA, 16 rue Claude Bernard, 75 231 Paris Cedex 05


bx <- function(x) {
        return(matrix(c(1, -x), 2, 1))
}
