# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function: 
# fssa20() and fssa30() functions calculates the relative 
# frequency distribution of anisotropic 2D & 3D clusters 
# with Moore (e,d)-neighborhood.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# n - sample size;
# x - linear dimension of the percolation lattice; 
# p0 - vector of p-values, distributed by lattice 
#      directions: (-x, +x, -y, +y, -z, +z);
# p1, p2 - double and triple combinations of p0-components, 
#          weighted by 2D & 3D (e,d)-neighborhood;
# set - vector of linear indexes of starting sites subset;
# all - trigger "Do we mark all starting sites or only accessible?";
# shape - vector of shape parameters of beta-distributed random variables, 
#         weighting the percolation lattice sites.
# Value:
# rfq - matrix of relative frequencies for sites of the percolation lattice.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
fssa2d <- function(n=1000, x=33, 
                   p0=runif(4, max=0.8), 
                   p1=colMeans(matrix(p0[c(
                     1,3, 2,3, 1,4, 2,4)], nrow=2))/2,
                   set=(x^2+1)/2, all=TRUE,
                   shape=c(1,1)) {
  rfq <- array(0, dim=rep(x, times=2))
  for (i in seq(n))
    rfq <- rfq + (ssa2d(x, p0, p1, set, all, shape) > 1)
  return(rfq/n)
}
fssa3d <- function(n=1000, x=33, 
                   p0=runif(6, max=0.4),
                   p1=colMeans(matrix(p0[c(
                     1,3, 2,3, 1,4, 2,4,
                     1,5, 2,5, 1,6, 2,6,
                     3,5, 4,5, 3,6, 4,6)], nrow=2))/2,
                   p2=colMeans(matrix(p0[c(
                     1,3,5, 2,3,5, 1,4,5, 2,4,5,
                     1,3,6, 2,3,6, 1,4,6, 2,4,6)], nrow=3))/3,
                   set=(x^3+1)/2, all=TRUE,
                   shape=c(1,1)) {
  rfq <- array(0, dim=rep(x, times=3))
  for (i in seq(n))
    rfq <- rfq + (ssa3d(x, p0, p1, p2, set, all, shape) > 1)
  return(rfq/n)
}