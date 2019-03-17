\name{fssa30}
\alias{fssa30}
\title{Frequency of Sites on a Square Anisotropic 3D lattice with (1,0)-neighborhood}
\description{
\code{fssa30()} function calculates the relative frequency distribution of anisotropic clusters on 3D square lattice with von Neumann (1,0)-neighborhood.
}
\usage{
fssa30(n=1000, x=33, p=runif(6, max=0.6), 
       set=(x^3+1)/2, all=TRUE, shape=c(1,1))
}
\arguments{
\item{n}{a sample size.}
\item{x}{a linear dimension of 3D square percolation lattice.}
\item{p}{a vector of relative fractions \code{(0<p)&(p<1)} of accessible sites (occupation probability) for lattice directions: \code{(-x,+x,-y,+y,-z,+z)}.}
\item{set}{a vector of linear indexes of a starting sites subset.}
\item{all}{logical; if \code{all=TRUE}, mark all sites from a starting subset; if \code{all=FALSE}, mark only accessible sites from a starting subset.}
\item{shape}{a vector with two shape parameters of beta-distributed random variables, weighting the percolation lattice sites.}
}
\value{
\item{rfq}{ 
a 3D matrix of relative sampling frequencies for sites of the percolation lattice.}
}
\details{
The percolation is simulated on 3D square lattice with uniformly weighted sites and the vector \code{p}, distributed over the lattice directions.

The anisotropic cluster is formed from the accessible sites connected with the initial subset \code{set}, and depends on the direction in 3D square lattice.

Von Neumann (1,0)-neighborhood on 3D square lattice consists of sites, only one coordinate of which is different from the current site by one: \code{e=c(-1,} \code{1,} \code{-x,} \code{x,} \code{-x^2,} \code{x^2)}.

Each element of the 3D matrix \code{frq} is equal to the relative frequency with which the 3D square lattice site belongs to a cluster sample of size \code{n}.
}
\seealso{
\link{ssa30}, \link{fssa20},
\link{fssi20}, \link{fssi30},
\link{fssa2d}, \link{fssa3d}
}
\concept{square lattice}
\concept{site percolation}
\concept{isotropic cluster}
\concept{von Neumann neighborhood}
\author{Pavel V. Moskalev <moskalefff@gmail.com>}
\references{
[1] Moskalev, P.V. Percolation modeling of porous structures. Moscow: URSS, 2018. 240 pp; in Russian.}
\examples{
x <- y <- seq(33)
rfq <- fssa30(n=200, p=.17*c(.5,3,.5,1.5,1,.5))
image(x, y, rfq[,,17], cex.main=1,
main="Frequencies in z=17 slice of anisotropic (1,0)-clusters")
contour(x, y, rfq[,,17], levels=seq(.05,.3,.05), add=TRUE)
abline(h=17, lty=2); abline(v=17, lty=2)
}
