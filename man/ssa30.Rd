\name{ssa30}
\alias{ssa30}
\title{Site cluster on Square Anisotropic 3D lattice with (1,0)-neighborhood}
\description{
\code{ssa30()} function provides sites labeling of the anisotropic cluster on 3D square lattice with von Neumann (1,0)-neighborhood.
}
\usage{
ssa30(x=33, p=runif(6, max=0.6), 
      set=(x^3+1)/2, all=TRUE, shape=c(1,1))
}
\arguments{
\item{x}{a linear dimension of 3D square percolation lattice.}
\item{p}{a vector of relative fractions \code{(0<p)&(p<1)} of accessible sites (occupation probability) for lattice directions: \code{(-x,+x,-y,+y,-z,+z)}.}
\item{set}{a vector of linear indexes of a starting sites subset.}
\item{all}{logical; if \code{all=TRUE}, mark all sites from a starting subset; if \code{all=FALSE}, mark only accessible sites from a starting subset.}
\item{shape}{a vector with two shape parameters of beta-distributed random variables, weighting the percolation lattice sites.}
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

The anisotropic cluster is formed from the accessible sites connected with the initial subset \code{set}, and depends on the direction in 3D square lattice.

To form the cluster the condition \code{acc[set+e[n]]<p[n]} is iteratively tested for sites of the von Neumann (1,0)-neighborhood \code{e} for the current cluster perimeter \code{set}, where \code{n} is equal to direction in 3D square lattice.

Von Neumann (1,0)-neighborhood on 3D square lattice consists of sites, only one coordinate of which is different from the current site by one: \code{e=c(-1,} \code{1,} \code{-x,} \code{x,} \code{-x^2,} \code{x^2)}.

Forming cluster ends with the exhaustion of accessible sites in von Neumann (1,0)-neighborhood of the current cluster perimeter.
}
\seealso{
\link{fssa30}, \link{ssa20},
\link{ssi20}, \link{ssi30},
\link{ssa2d}, \link{ssa3d}
}
\concept{square lattice}
\concept{site percolation}
\concept{anisotropic cluster}
\concept{von Neumann neighborhood}
\author{Pavel V. Moskalev}
\references{
[1] Moskalev, P.V. Percolation modeling of porous structures. Moscow: URSS, 2018. 240 pp; in Russian.\cr
[2] Moskalev, P.V. (2013) The structure of site percolation models on three-dimensional square lattices. \emph{Computer Research and Modeling}, Vol.5, No.4, pp.607–622; in Russian.}
\examples{
# Example No.1. Axonometric projection of 3D cluster
require(lattice)
set.seed(20120521)
x <- y <- z <- seq(33)
cls <- which(ssa30(p=.09*c(1,6,1,3,2,1))>1, arr.ind=TRUE)
cloud(cls[,3] ~ cls[,1]*cls[,2],
xlim=range(x), ylim=range(y), zlim=range(z), 
col=rgb(1,0,0,0.4), xlab="x", ylab="y", zlab="z", main.cex=1,
main="Anisotropic (1,0)-cluster")

# Example No.2. Z=17 slice of 3D cluster
set.seed(20120521)
x <- y <- z <- seq(33)
cls <- ssa30(p=.09*c(1,6,1,3,2,1))
image(x, y, cls[,,17], zlim=c(0,2), cex.main=1, 
main="Z=17 slice of an anisotropic (1,0)-cluster")
abline(h=17, lty=2); abline(v=17, lty=2)
}