# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function: 
# ssa2d() and ssa3d() functions provide a labeling of 
# anisotropic 2D & 3D clusters with d-neighborhood.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# x - linear dimension of the percolation lattice; 
# p0 - vector of p-values, distributed by lattice 
#      directions: (-x, +x, -y, +y, -z, +z);
# p1, p2 - double and triple combinations of p0-components, 
#          weighted by 2D & 3D d-neighborhood;
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
ssa2d <- function(x=33, 
                  p0=runif(4, max=0.8), 
                  p1=colMeans(matrix(p0[c(
                    1,3, 2,3, 1,4, 2,4)], nrow=2))/2,
                  set=(x^2+1)/2, all=TRUE,
                  shape=c(1,1)) {
  e0 <- c(-1, 1,-x, x)
  e1 <- colSums(matrix(e0[c(
    1,3, 2,3, 1,4, 2,4)], nrow=2))
  acc <- array(rbeta(x^2,shape[1],shape[2]), rep(x,2))
  if (all) acc[set] <- 2
  else acc[set <- set[acc[set] < mean(p0)]] <- 2
  acc[c(1,x),] <- acc[,c(1,x)] <- 1
  repeat {
    acc[set <- unique(c(
        set[acc[set+e0[1]] < p0[1]] +e0[1],
        set[acc[set+e0[2]] < p0[2]] +e0[2],
        set[acc[set+e0[3]] < p0[3]] +e0[3],
        set[acc[set+e0[4]] < p0[4]] +e0[4],
        set[acc[set+e1[1]] < p1[1]] +e1[1],
        set[acc[set+e1[2]] < p1[2]] +e1[2],
        set[acc[set+e1[3]] < p1[3]] +e1[3],
        set[acc[set+e1[4]] < p1[4]] +e1[4]))] <- 2
    if (length(set) < 1) break
  }
  return(acc)
}
ssa3d <- function(x=33, 
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
  else acc[set <- set[acc[set] < mean(p0)]] <- 2
  acc[c(1,x),,] <- acc[,c(1,x),] <- acc[,,c(1,x)] <- 1
  repeat {
    acc[set <- unique(c(
        set[acc[set+e0[1]] < p0[1]] +e0[1],
        set[acc[set+e0[2]] < p0[2]] +e0[2],
        set[acc[set+e0[3]] < p0[3]] +e0[3],
        set[acc[set+e0[4]] < p0[4]] +e0[4],
        set[acc[set+e0[5]] < p0[5]] +e0[5],
        set[acc[set+e0[6]] < p0[6]] +e0[6],
        set[acc[set+e1[1]] < p1[1]] +e1[1],
        set[acc[set+e1[2]] < p1[2]] +e1[2],
        set[acc[set+e1[3]] < p1[3]] +e1[3],
        set[acc[set+e1[4]] < p1[4]] +e1[4],
        set[acc[set+e1[5]] < p1[5]] +e1[5],
        set[acc[set+e1[6]] < p1[6]] +e1[6],
        set[acc[set+e1[7]] < p1[7]] +e1[7],
        set[acc[set+e1[8]] < p1[8]] +e1[8],
        set[acc[set+e1[9]] < p1[9]] +e1[9],
        set[acc[set+e1[10]]< p1[10]]+e1[10],
        set[acc[set+e1[11]]< p1[11]]+e1[11],
        set[acc[set+e1[12]]< p1[12]]+e1[12],
        set[acc[set+e2[1]] < p2[1]] +e2[1],
        set[acc[set+e2[2]] < p2[2]] +e2[2],
        set[acc[set+e2[3]] < p2[3]] +e2[3],
        set[acc[set+e2[4]] < p2[4]] +e2[4],
        set[acc[set+e2[5]] < p2[5]] +e2[5],
        set[acc[set+e2[6]] < p2[6]] +e2[6],
        set[acc[set+e2[7]] < p2[7]] +e2[7],
        set[acc[set+e2[8]] < p2[8]] +e2[8]))] <- 2
    if (length(set) < 1) break
  }  
  return(acc)
}