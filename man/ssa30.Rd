\name{ssa30}
\alias{ssa30}
\title{Site cluster on Square Anisotropic 3D lattice with 0-neighborhood}
\description{
\code{ssa30()} function provide a sites labeling of the anisotropic cluster on 3D square lattice with von Neumann 0-neighborhood.
}
\usage{
ssa30(x=33, p=runif(6, max=0.6), set=(x^3+1)/2, all=TRUE)
}
\arguments{
\item{x}{ a linear dimension of 3D square percolation lattice. }
\item{p}{ a vector of relative fractions \code{(0<p)&(p<1)} of accessible sites (occupation probability) for lattice directions: \code{(-x,+x,-y,+y,-z,+z)}. }
\item{set}{ a vector of linear indexes of initial sites subset. }
\item{all}{ logical; if \code{all=TRUE}, mark all sites from initial subset; if \code{all=FALSE}, mark accessible sites from initial subset. }
}
\value{
\item{acc}{ 
an accessibility matrix for 3D square percolation lattice:\cr
if \code{acc[e]<p[n]} then \code{acc[e]} is accessible site;\cr
if \code{acc[e]==1} then \code{acc[e]} is non-accessible site;\cr
if \code{acc[e]==2} then \code{acc[e]} belongs to a sites cluster. }
}
\details{
The percolation is simulated on 3D square lattice with uniformly weighted sites \code{acc} and the vector \code{p}, distributed over the lattice directions.

The anisotropic cluster is formed from the accessible sites connected with the initial subset, and depends on the direction in 3D square lattice.

To form the cluster iteratively tested the condition \code{acc[set+e[n]]<p[n]} for sites of the von Neumann 0-neighborhood \code{e} for the current cluster perimeter \code{set}, where \code{n} is equal to direction in 3D square lattice.

Von Neumann 0-neighborhood on 3D square lattice consists of sites, only one coordinate which is different from the current site by one: \code{e=c(-1,} \code{1,} \code{-x,} \code{x,} \code{-x^2,} \code{x^2)}.

Forming cluster ends with the exhaustion of accessible sites in von Neumann 0-neighborhood of the current cluster perimeter.
}
\seealso{
\link{ssa20},
\link{ssi20}, \link{ssi30},
\link{ssa2d}, \link{ssa3d}
}
\concept{square lattice}
\concept{site percolation}
\concept{anisotropic cluster}
\concept{von Neumann neighborhood}
\concept{von Neumann 0-neighborhood}
\author{Pavel V. Moskalev <moskalefff@gmail.com>}
\examples{
require(lattice)
set.seed(20120506)
x <- y <- z <- seq(33)
cls <- which(ssa30()>1, arr.ind=TRUE)
cloud(cls[,3] ~ cls[,1]*cls[,2],
      xlim=range(x), ylim=range(y), zlim=range(z), 
      col=rgb(1,0,0,0.4), xlab="x", ylab="y", zlab="z",
      main="Anisotropic 3D cluster with 0-neighborhood")
}