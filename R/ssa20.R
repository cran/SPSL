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
  e <- c(-1, 1,-x, x)
  acc <- array(runif(x^2), rep(x,2))
  if (all) acc[set] <- 2
  else acc[set <- set[acc[set] < mean(p)]] <- 2
  acc[c(1,x),] <- acc[,c(1,x)] <- 1
  repeat {
    acc[set <- unique(c(
      set[acc[set+e[1]] < p[1]] +e[1],
      set[acc[set+e[2]] < p[2]] +e[2],
      set[acc[set+e[3]] < p[3]] +e[3],
      set[acc[set+e[4]] < p[4]] +e[4]))] <- 2
    if (length(set) < 1) break
  }
  return(acc)
}
ssa30 <- function(x=33, 
                  p=runif(6, max=0.6), 
                  set=(x^3+1)/2, all=TRUE) {
  e <- c(-1, 1,-x, x,-x^2, x^2)
  acc <- array(runif(x^3), rep(x,3))
  if (all) acc[set] <- 2
  else acc[set <- set[acc[set] < mean(p)]] <- 2
  acc[c(1,x),,] <- acc[,c(1,x),] <- acc[,,c(1,x)] <- 1
  repeat {
    acc[set <- unique(c(
      set[acc[set+e[1]] < p[1]] +e[1],
      set[acc[set+e[2]] < p[2]] +e[2],
      set[acc[set+e[3]] < p[3]] +e[3],
      set[acc[set+e[4]] < p[4]] +e[4],
      set[acc[set+e[5]] < p[5]] +e[5],      
      set[acc[set+e[6]] < p[6]] +e[6]))] <- 2
    if (length(set) < 1) break
  }
  return(acc)
}