# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function: 
# ssi20() and ssi30() functions provide a labeling of 
# isotropic 2D & 3D clusters with von Neumann neighborhood.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# x - linear dimension of the percolation lattice; 
# p - relative fraction of accessible sites 
#     (occupation probability) for percolation lattice;
# set - vector of linear indexes of initial sites subset;
# all - trigger "Mark all initial sites or accessible only?"
# Variables:
# e - linear indexes of sites from 2D & 3D von Neumann neighborhood.
# Value:
# acc - labeled accessibility matrix for the percolation lattice.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
ssi20 <- function(x=33, 
                  p=0.592746, 
                  set=(x^2+1)/2, all=TRUE) {
  e <- c(-1, 1,-x, x)
  acc <- array(runif(x^2), rep(x,2))
  if (all) acc[set] <- 2
  else acc[set <- set[acc[set] < p]] <- 2
  acc[c(1,x),] <- acc[,c(1,x)] <- 1
  repeat {
    acc[set <- unique(c(
      set[acc[set+e[1]] < p] +e[1],
      set[acc[set+e[2]] < p] +e[2],
      set[acc[set+e[3]] < p] +e[3],
      set[acc[set+e[4]] < p] +e[4] ))] <- 2
    if (length(set) < 1) break
  }
  return(acc)
}
ssi30 <- function(x=33, 
                  p=0.311608, 
                  set=(x^3+1)/2, all=TRUE) {
  e <- c(-1, 1,-x, x,-x^2, x^2)
  acc <- array(runif(x^3), rep(x,3))
  if (all) acc[set] <- 2
  else acc[set <- set[acc[set] < p]] <- 2
  acc[c(1,x),,] <- acc[,c(1,x),] <- acc[,,c(1,x)] <- 1
  repeat {
    acc[set <- unique(c(
      set[acc[set+e[1]] < p] +e[1],
      set[acc[set+e[2]] < p] +e[2],
      set[acc[set+e[3]] < p] +e[3],
      set[acc[set+e[4]] < p] +e[4],
      set[acc[set+e[5]] < p] +e[5],
      set[acc[set+e[6]] < p] +e[6] ))] <- 2
    if (length(set) < 1) break
  }
  return(acc)
}