\name{ssa20}
\alias{ssa20}
\title{Site cluster on Square Anisotropic 2D lattice with (1,0)-neighborhood}
\description{
\code{ssa20()} function provides sites labeling of the anisotropic cluster on 2D square lattice with von Neumann (1,0)-neighborhood.
}
\usage{
ssa20(x=33, p=runif(4, max=0.9), set=(x^2+1)/2, all=TRUE)
}
\arguments{
\item{x}{a linear dimension of 2D square percolation lattice.}
\item{p}{a vector of relative fractions \code{(0<p)&(p<1)} of accessible sites (occupation probability) for lattice directions: \code{(-x,+x,-y,+y)}.}
\item{set}{a vector of linear indexes of initial sites subset.}
\item{all}{logical; if \code{all=TRUE}, mark all sites from initial subset; if \code{all=FALSE}, mark accessible sites from initial subset.}
}
\value{
\item{acc}{ 
an accessibility matrix for 2D square percolation lattice:\cr
if \code{acc[e]<p[n]} then \code{acc[e]} is accessible site;\cr
if \code{acc[e]==1} then \code{acc[e]} is non-accessible site;\cr
if \code{acc[e]==2} then \code{acc[e]} belongs to a sites cluster. }
}
\details{
The percolation is simulated on 2D square lattice with uniformly weighted sites \code{acc} and the vector \code{p}, distributed over the lattice directions.

The anisotropic cluster is formed from the accessible sites connected with the initial subset, and depends on the direction in 2D square lattice.

To form the cluster the condition \code{acc[set+e[n]]<p[n]} is iteratively tested for sites of the von Neumann (1,0)-neighborhood \code{e} for the current cluster perimeter \code{set}, where \code{n} is equal to direction in 2D square lattice.

Von Neumann (1,0)-neighborhood on 2D square lattice consists of sites, only one coordinate of which is different from the current site by one: \code{e=c(-1,} \code{1,} \code{-x,} \code{x)}.

Forming cluster ends with the exhaustion of accessible sites in von Neumann (1,0)-neighborhood of the current cluster perimeter.
}
\seealso{
\link{fssa20}, \link{ssa30},
\link{ssi20}, \link{ssi30},
\link{ssa2d}, \link{ssa3d}
}
\concept{square lattice}
\concept{site percolation}
\concept{anisotropic cluster}
\concept{von Neumann neighborhood}
\author{Pavel V. Moskalev}
\references{
Moskalev, P.V. and Shitov, V.V. Mathematical modeling of porous structures. Moscow: Fizmatlit, 2007. 120 pp; in Russian.
}
\examples{
set.seed(20120507)
x <- y <- seq(33)
image(x, y, ssa20(), zlim=c(0,2), 
main="Anisotropic 2D cluster with (1,0)-neighborhood")
abline(h=17, lty=2); abline(v=17, lty=2)
}