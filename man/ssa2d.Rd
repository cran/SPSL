\name{ssa2d}
\alias{ssa2d}
\title{Site cluster on Square Anisotropic 2D lattice with d-neighborhood}
\description{
\code{ssa2d()} function provide a sites labeling of the anisotropic cluster on 2D square lattice with Moore d-neighborhood.
}
\usage{
ssa2d(x=33, p0=runif(4, max=0.8), 
            p1=colMeans(matrix(p0[c(
               1,3, 2,3, 1,4, 2,4)], nrow=2))/2,
            set=(x^2+1)/2, all=TRUE)
}
\arguments{
\item{x}{ a linear dimension of 2D square percolation lattice. }
\item{p0}{ a vector of relative fractions \code{(0<p0)&(p0<1)} of accessible sites (occupation probability) for lattice directions: \code{(-x,+x,-y,+y)}. }
\item{p1}{ averaged double combinations of \code{p0}-components weighted by non-metrical Minkowski distance: \code{p1=colMeans(matrix(p0[c(1,3,...)], nrow=2))/rho.mink.e1}. }
\item{set}{ a vector of linear indexes of initial sites subset. }
\item{all}{ logical; if \code{all=TRUE}, mark all sites from initial subset; if \code{all=FALSE}, mark accessible sites from initial subset. }
}
\value{
\item{acc}{ 
an accessibility matrix for 2D square percolation lattice:\cr
if \code{acc[e]<pN[n]} then \code{acc[e]} is accessible site;\cr
if \code{acc[e]==1} then \code{acc[e]} is non-accessible site;\cr
if \code{acc[e]==2} then \code{acc[e]} belongs to a sites cluster. }
}
\details{
The percolation is simulated on 2D square lattice with uniformly weighted sites \code{acc} and the vectors \code{p0} and \code{p1}, distributed over the lattice directions, and their combinations.

The anisotropic cluster is formed from the accessible sites connected with the initial subset, and depends on the direction in 2D square lattice.

To form the cluster iteratively tested the condition \code{acc[set+eN[n]]<pN[n]} for sites of the Moore d-neighborhood \code{eN} for the current cluster perimeter \code{set}, where \code{eN} is equal to \code{e0} or \code{e1} vector; \code{pN} is equal to \code{p0} or \code{p1} vector; \code{n} is equal to direction in 2D square lattice.

Moore d-neighborhood on 2D square lattice consists of sites, at least one coordinate which is different from the current site by one: \code{e=c(e0,e1)}, where 
\code{e0=c(-1,} \code{1,} \code{-x,} \code{x)}; 
\code{e1=colSums(matrix(e0[c(1,3,} \code{2,3,} \code{1,4,} \code{2,4)], nrow=2))}.

Minkowski non-metrical distance between sites \code{a} and \code{b} depend on the exponent \code{d}:\cr
\code{rho.mink <- function(a, b, d=1)} \cr
\code{if (is.infinite(d)) return(apply(abs(b-a), 2, max))} \cr
\code{else return(apply(abs(b-a)^d, 2, sum)^(1/d))}.

Minkowski non-metrical distance for a sites from \code{e1} subset with the exponent \code{d=1} is equal to \code{rho.mink.e1=2}.

Forming cluster ends with the exhaustion of accessible sites in Moore d-neighborhood of the current cluster perimeter.
}
\seealso{
\link{ssa3d},
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
set.seed(20120507)
x <- y <- seq(33)
image(x, y, ssa2d(), asp=1, zlim=c(0,2), 
      main="Anisotropic 2D cluster with 1-neighborhood")
}
