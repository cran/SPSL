# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# License: GPL-3
# Package: Site Percolation on Square Lattice (SPSL)
# Author: Pavel V. Moskalev <moskalefff@gmail.com>
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function: 
# ssa20() and ssa30() functions provide a labeling of 
# anisotropic 2D & 3D clusters with von Neumann neighborhood.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# x - linear dimension of the percolation lattice; 
# p - vector of p[1:6] values, distributed by lattice 
#     directions: -x, +x, -y, +y, -z, +z;
# set - vector of linear indexes of initial sites subset;
# all - trigger "Mark all initial sites or accessible only?"
# Variables:
# e - linear indexes of sites from 2D & 3D von Neumann neighborhood.
# Value:
# acc - labeled accessibility matrix for the percolation lattice.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
ssa20 <- function(x=33, 
                  p=runif(4, max=0.9), 
                  set=(x^2+1)/2, all=TRUE) {
  p <- as.double(p)
  b <- as.integer(length(set))
  e <- as.integer(c(-1, 1,-x, x))
  acc <- array(runif(x^2), rep(x,2))
  acc[c(1,x),] <- acc[,c(1,x)] <- 1
  if (all) {
    cls <- rep(0L, b + max(p)*x^2) 
    acc[set] <- 2 
  } else {
    cls <- rep(0L, max(p)*x^2) 
    acc[set <- set[acc[set] < mean(p)]] <- 2
  }
  cls[seq_along(set)] <- as.integer(set - 1)
  .Call("ssTNd", p, acc, b, e, cls) 
  return(acc)
}
ssa30 <- function(x=33, 
                  p=runif(6, max=0.6), 
                  set=(x^3+1)/2, all=TRUE) {
  p <- as.double(p)
  b <- as.integer(length(set))
  e <- as.integer(c(-1, 1,-x, x,-x^2, x^2))
  acc <- array(runif(x^3), rep(x,3))
  acc[c(1,x),,] <- acc[,c(1,x),] <- acc[,,c(1,x)] <- 1
  if (all) {
    cls <- rep(0L, b + max(p)*x^3) 
    acc[set] <- 2 
  } else {
    cls <- rep(0L, max(p)*x^3) 
    acc[set <- set[acc[set] < mean(p)]] <- 2
  }
  cls[seq_along(set)] <- as.integer(set - 1)
  .Call("ssTNd", p, acc, b, e, cls) 
  return(acc)
}