# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function: 
# ssi2d() and ssi3d() functions provide a labeling of 
# isotropic 2D & 3D clusters with Moore d-neighborhood.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# x - linear dimension of the percolation lattice; 
# p0 - relative fraction of accessible sites
#      (occupation probability) for percolation lattice;
# p1, p2 - p1 value, weighted by 2D & 3D d-neighborhood;
# set - vector of linear indexes of starting sites subset;
# all - trigger "Do we mark all starting sites or only accessible?";
# shape - vector of shape parameters of beta-distributed random variables, 
#         weighting the percolation lattice sites.
# Variables:
# e0, e1, e2 - linear indexes of sites combinations from 
#              2D & 3D Moore neighborhood.
# Value:
# acc - labeled accessibility matrix for the percolation lattice.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
ssi2d <- function(x=33, 
                  p0=0.5, p1=p0/2, 
                  set=(x^2+1)/2, all=TRUE,
                  shape=c(1,1)) {
  e0 <- c(-1, 1,-x, x)
  e1 <- colSums(matrix(e0[c(
    1,3, 2,3, 1,4, 2,4)], nrow=2))
  acc <- array(rbeta(x^2,shape[1],shape[2]), rep(x,2))
  if (all) acc[set] <- 2
  else acc[set <- set[acc[set] < p0]] <- 2
  acc[c(1,x),] <- acc[,c(1,x)] <- 1
  repeat {
    acc[set <- unique(c(
      set[acc[set+e0[1]] < p0] +e0[1],
      set[acc[set+e0[2]] < p0] +e0[2],
      set[acc[set+e0[3]] < p0] +e0[3],
      set[acc[set+e0[4]] < p0] +e0[4],
      set[acc[set+e1[1]] < p1] +e1[1],
      set[acc[set+e1[2]] < p1] +e1[2],
      set[acc[set+e1[3]] < p1] +e1[3],
      set[acc[set+e1[4]] < p1] +e1[4] ))] <- 2
    if (length(set) < 1) break
  }
  return(acc)
}
ssi3d <- function(x=33, 
                  p0=0.2, p1=p0/2, p2=p0/3,
                  set=(x^3+1)/2, all=TRUE,
                  shape=c(1,1)) {
  e0 <- c(-1, 1,-x, x,-x^2, x^2)
  e1 <- colSums(matrix(e0[c(
    1,3, 2,3, 1,4, 2,4,
    1,5, 2,5, 1,6, 2,6,
    3,5, 4,5, 3,6, 4,6)], nrow=2))
  e2 <- colSums(matrix(e0[c(
    1,3,5, 2,3,5, 1,4,5, 2,4,5,
    1,3,6, 2,3,6, 1,4,6, 2,4,6)], nrow=3))
  acc <- array(rbeta(x^3,shape[1],shape[2]), rep(x,3))
  if (all) acc[set] <- 2
  else acc[set <- set[acc[set] < p0]] <- 2
  acc[c(1,x),,] <- acc[,c(1,x),] <- acc[,,c(1,x)] <- 1
  repeat {
    acc[set <- unique(c(
      set[acc[set+e0[1]] < p0] +e0[1],
      set[acc[set+e0[2]] < p0] +e0[2],
      set[acc[set+e0[3]] < p0] +e0[3],
      set[acc[set+e0[4]] < p0] +e0[4],
      set[acc[set+e0[5]] < p0] +e0[5],
      set[acc[set+e0[6]] < p0] +e0[6],
      set[acc[set+e1[1]] < p1] +e1[1],
      set[acc[set+e1[2]] < p1] +e1[2],
      set[acc[set+e1[3]] < p1] +e1[3],
      set[acc[set+e1[4]] < p1] +e1[4],
      set[acc[set+e1[5]] < p1] +e1[5],
      set[acc[set+e1[6]] < p1] +e1[6],
      set[acc[set+e1[7]] < p1] +e1[7],
      set[acc[set+e1[8]] < p1] +e1[8],
      set[acc[set+e1[9]] < p1] +e1[9],
      set[acc[set+e1[10]]< p1] +e1[10],
      set[acc[set+e1[11]]< p1] +e1[11],
      set[acc[set+e1[12]]< p1] +e1[12],
      set[acc[set+e2[1]] < p2] +e2[1],
      set[acc[set+e2[2]] < p2] +e2[2],
      set[acc[set+e2[3]] < p2] +e2[3],
      set[acc[set+e2[4]] < p2] +e2[4],
      set[acc[set+e2[5]] < p2] +e2[5],
      set[acc[set+e2[6]] < p2] +e2[6],
      set[acc[set+e2[7]] < p2] +e2[7],
      set[acc[set+e2[8]] < p2] +e2[8] ))] <- 2
    if (length(set) < 1) break
  }
  return(acc)
}