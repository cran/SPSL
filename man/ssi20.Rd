\name{ssi20}
\alias{ssi20}
\title{Site cluster on Square Isotropic 2D lattice with 0-neighborhood}
\description{
\code{ssi20()} function provide a sites labeling of the isotropic cluster on 2D square lattice with von Neumann 0-neighborhood.
}
\usage{
ssi20(x=33, p=0.592746, set=(x^2+1)/2, all=TRUE)
}
\arguments{
\item{x}{ a linear dimension of 2D square percolation lattice. }
\item{p}{ the relative fractions \code{(0<p)&(p<1)} of accessible sites (occupation probability) for percolation lattice. }
\item{set}{ a vector containing the linear indexes of sites from initial subset. }
\item{all}{ logical; if \code{all=TRUE}, mark all sites from initial subset; if \code{all=FALSE}, mark accessible sites from initial subset. }
}
\value{
\item{acc}{ 
an accessiblity matrix for 2D square percolation lattice:\cr 
if \code{acc[e]<p} then \code{acc[e]} is accessible site;\cr 
if \code{acc[e]==1} then \code{acc[e]} is non-accessible site;\cr 
if \code{acc[e]==2} then \code{acc[e]} belongs to a sites cluster. }
}
\details{
The percolation is simulated on 2D square lattice with uniformly weighted sites \code{acc} and the constant parameter \code{p}.

The isotropic cluster is formed from the accessible sites connected with initial sites subset. 

To form the cluster iteratively tested the condition \code{acc[set+e]<p} for sites of the von Neumann 0-neighborhood \code{e} for the current cluster perimeter \code{set}.

Von Neumann 0-neighborhood on 2D square lattice consists of sites, only one coordinate which is different from the current site by one: \code{e=c(-1,} \code{1,} \code{-x,} \code{x)}.

Forming cluster ends with the exhaustion of accessible sites in von Neumann 0-neighborhood of the current cluster perimeter.
}
\seealso{
\link{ssi30},
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
set.seed(20120507)
x <- y <- seq(33)
image(x, y, ssi20(), asp=1, zlim=c(0,2))
}
