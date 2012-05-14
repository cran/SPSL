\name{ssi30}
\alias{ssi30}
\title{Site cluster on Square Isotropic 3D lattice with 0-neighborhood}
\description{
\code{ssi30()} function provide a sites labeling of the isotropic cluster on 3D square lattice with von Neumann 0-neighborhood.
}
\usage{
ssi30(x=33, p=0.311608, set=(x^3+1)/2, all=TRUE)
}
\arguments{
\item{x}{ a linear dimension of 3D square percolation lattice. }
\item{p}{ the relative fractions \code{(0<p)&(p<1)} of accessible sites (occupation probability) for percolation lattice. }
\item{set}{ a vector containing the linear indexes of sites from initial subset. }
\item{all}{ logical; if \code{all=TRUE}, mark all sites from initial subset; if \code{all=FALSE}, mark accessible sites from initial subset. }
}
\value{
\item{acc}{ 
an accessiblity matrix for 3D square percolation lattice:\cr 
if \code{acc[e]<p} then \code{acc[e]} is accessible site;\cr 
if \code{acc[e]==1} then \code{acc[e]} is non-accessible site;\cr 
if \code{acc[e]==2} then \code{acc[e]} belongs to a sites cluster. }
}
\details{
The percolation is simulated on 3D square lattice with uniformly weighted sites \code{acc} and the constant parameter \code{p}.

The isotropic cluster is formed from the accessible sites connected with initial sites subset. 

To form the cluster iteratively tested the condition \code{acc[set+e]<p} for sites of the von Neumann 0-neighborhood \code{e} for the current cluster perimeter \code{set}.

Von Neumann 0-neighborhood on 3D square lattice consists of sites, only one coordinate which is different from the current site by one: \code{e=c(-1,} \code{1,} \code{-x,} \code{x,} \code{-x^2,} \code{x^2)}.

Forming cluster ends with the exhaustion of accessible sites in von Neumann 0-neighborhood of the current cluster perimeter.
}
\seealso{
\link{ssi20},
\link{ssa20}, \link{ssa30},
\link{ssi2d}, \link{ssi3d}
}
\concept{square lattice}
\concept{site percolation}
\concept{isotropic cluster}
\concept{von Neumann neighborhood}
\concept{von Neumann 0-neighborhood}
\author{Pavel V. Moskalev <moskalefff@gmail.com>}
\examples{
require(lattice)
set.seed(20120507)
x <- y <- z <- seq(33)
cls <- which(ssi30()>1, arr.ind=TRUE)
cloud(cls[,3] ~ cls[,1]*cls[,2],
      xlim=range(x), ylim=range(y), zlim=range(z), 
      col=rgb(1,0,0,0.4), xlab="x", ylab="y", zlab="z",
      main="Isotropic 3D cluster with 0-neighborhood")
}
