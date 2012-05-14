\name{ssa3d}
\alias{ssa3d}
\title{Site cluster on Square Anisotropic 3D lattice with d-neighborhood}
\description{
\code{ssa3d()} function provide a sites labeling of the anisotropic cluster on 3D square lattice with Moore d-neighborhood.
}
\usage{
ssa3d(x=33, p0=runif(6, max=0.4),
            p1=colMeans(matrix(p0[c(
               1,3, 2,3, 1,4, 2,4,
               1,5, 2,5, 1,6, 2,6,
               3,5, 4,5, 3,6, 4,6)], nrow=2))/2,
            p2=colMeans(matrix(p0[c(
               1,3,5, 2,3,5, 1,4,5, 2,4,5,
               1,3,6, 2,3,6, 1,4,6, 2,4,6)], nrow=3))/3,                 
            set=(x^3+1)/2, all=TRUE)
}
\arguments{
\item{x}{ a linear dimension of 3D square percolation lattice. }
\item{p0}{ a vector of relative fractions \code{(0<p0)&(p0<1)} of accessible sites (occupation probability) for lattice directions: \code{(-x,+x,-y,+y,-z,+z)}. }
\item{p1}{ averaged double combinations of \code{p0}-components weighted by non-metrical Minkowski distance: \code{p1=colMeans(matrix(p0[c(1,3,...)], nrow=2))/rho.mink.e1}. }
\item{p2}{ averaged triple combinations of \code{p0}-components weighted by non-metrical Minkowski distance: \code{p2=colMeans(matrix(p0[c(1,3,5,...)], nrow=3))/rho.mink.e2}. }
\item{set}{ a vector of linear indexes of initial sites subset. }
\item{all}{ logical; if \code{all=TRUE}, mark all sites from initial subset; if \code{all=FALSE}, mark accessible sites from initial subset. }
}
\value{
\item{acc}{ 
an accessibility matrix for 3D square percolation lattice:\cr
if \code{acc[e]<pN[n]} then \code{acc[e]} is accessible site;\cr
if \code{acc[e]==1} then \code{acc[e]} is non-accessible site;\cr
if \code{acc[e]==2} then \code{acc[e]} belongs to a sites cluster. }
}
\details{
The percolation is simulated on 3D square lattice with uniformly weighted sites \code{acc} and the vectors \code{p0}, \code{p1}, and \code{p2}, distributed over the lattice directions, and their combinations.

The anisotropic cluster is formed from the accessible sites connected with the initial subset, and depends on the direction in 3D square lattice.

To form the cluster iteratively tested the condition \code{acc[set+eN[n]]<pN[n]} for sites of the Moore d-neighborhood \code{eN} for the current cluster perimeter \code{set}, where \code{eN} is equal to \code{e0}, \code{e1}, or \code{e2} vector; \code{pN} is equal to \code{p0}, \code{p1}, or \code{p2} vector; \code{n} is equal to direction in 3D square lattice.

Moore d-neighborhood on 3D square lattice consists of sites, at least one coordinate which is different from the current site by one: \code{e=c(e0,e1,e2)}, where 
\code{e0=c(-1,} \code{1,} \code{-x,} \code{x,} \code{-x^2,} \code{x^2)}; 
\code{e1=colSums(matrix(e0[c(1,3,} \code{2,3,} \code{1,4,} \code{2,4,} \code{1,5,} \code{2,5,} \code{1,6,} \code{2,6,} \code{3,5,} \code{4,5,} \code{3,6,} \code{4,6)], nrow=2))};
\code{e2=colMeans(matrix(p0[c(1,3,5,} \code{2,3,5,} \code{1,4,5,} \code{2,4,5,} \code{1,3,6,} \code{2,3,6,} \code{1,4,6,} \code{2,4,6)], nrow=3))}.

Minkowski non-metrical distance between sites \code{a} and \code{b} depend on the exponent \code{d}:\cr
\code{rho.mink <- function(a, b, d=1)} \cr
\code{if (is.infinite(d)) return(apply(abs(b-a), 2, max))} \cr
\code{else return(apply(abs(b-a)^d, 2, sum)^(1/d))}.

Minkowski non-metrical distance for a sites from \code{e1} and \code{e2} subsets with the exponent \code{d=1} is equal to \code{rho.mink.e1=2} and \code{rho.mink.e2=3}.

Forming cluster ends with the exhaustion of accessible sites in Moore d-neighborhood of the current cluster perimeter.
}
\seealso{
\link{ssa2d},
\link{ssa20}, \link{ssa30},
\link{ssi2d}, \link{ssi3d}
}
\concept{square lattice}
\concept{site percolation}
\concept{anisotropic cluster}
\concept{Moore neighborhood}
\concept{Moore d-neighborhood}
\concept{non-metrical Minkowski distance}
\author{Pavel V. Moskalev <moskalefff@gmail.com>}
\examples{
require(lattice)
set.seed(20120505)
x <- y <- z <- seq(33)
cls <- which(ssa3d()>1, arr.ind=TRUE)
cloud(cls[,3] ~ cls[,1]*cls[,2],
      xlim=range(x), ylim=range(y), zlim=range(z), 
      col=rgb(1,0,0,0.4), xlab="x", ylab="y", zlab="z",
      main="Anisotropic 3D cluster with 1-neighborhood")
}
