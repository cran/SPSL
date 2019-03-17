# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function: 
# fssa20() and fssa30() functions calculates the relative 
# frequency distribution of anisotropic 2D & 3D clusters 
# with von Neumann neighborhood.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# n - sample size;
# x - linear dimension of the percolation lattice; 
# p - vector of p[1:6] values, distributed by lattice 
#     directions: -x, +x, -y, +y, -z, +z;
# set - vector of linear indexes of starting sites subset;
# all - trigger "Do we mark all starting sites or only accessible?";
# shape - vector of shape parameters of beta-distributed random variables, 
#         weighting the percolation lattice sites.
# Value:
# rfq - matrix of relative frequencies for sites of the percolation lattice.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
fssa20 <- function(n=1000, 
                   x=33, p=runif(4, max=0.9),
                   set=(x^2+1)/2, all=TRUE,
                   shape=c(1,1)) {
  rfq <- array(0, dim=rep(x, times=2))
  for (i in seq(n))
    rfq <- rfq + (ssa20(x, p, set, all, shape) > 1)
  return(rfq/n)
}
fssa30 <- function(n=1000,
                   x=33, p=runif(6, max=0.6), 
                   set=(x^3+1)/2, all=TRUE,
                   shape=c(1,1)) {
  rfq <- array(0, dim=rep(x, times=3))
  for (i in seq(n))
    rfq <- rfq + (ssa30(x, p, set, all, shape) > 1)
  return(rfq/n)
}