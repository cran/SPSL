\name{fssi30}
\alias{fssi30}
\title{Frequency of Sites on a Square Isotropic 3D lattice with (1,0)-neighborhood}
\description{
\code{fssi30()} function calculates the relative frequency distribution of isotropic clusters on 3D square lattice with von Neumann (1,0)-neighborhood.
}
\usage{
fssi30(n=1000, x=33, p=0.311608, set=(x^3+1)/2, all=TRUE)
}
\arguments{
\item{n}{a sample size.}
\item{x}{a linear dimension of 3D square percolation lattice.}
\item{p}{the relative fractions \code{(0<p)&(p<1)} of accessible sites (occupation probability) for percolation lattice.}
\item{set}{a vector containing the linear indexes of sites from initial subset.}
\item{all}{logical; if \code{all=TRUE}, mark all sites from initial subset; if \code{all=FALSE}, mark accessible sites from initial subset.}
}
\value{
\item{rfq}{ 
a 3D matrix of relative sampling frequencies for sites of the percolation lattice.}
}
\details{
The percolation is simulated on 3D square lattice with uniformly weighted sites and the constant parameter \code{p}.

The isotropic cluster is formed from the accessible sites connected with initial sites subset \code{set}.

Von Neumann (1,0)-neighborhood on 3D square lattice consists of sites, only one coordinate of which is different from the current site by one: \code{e=c(-1,} \code{1,} \code{-x,} \code{x,} \code{-x^2,} \code{x^2)}.

Each element of the matrix \code{frq} is equal to the relative frequency with which the 3D square lattice site belongs to a cluster sample of size \code{n}.
}
\seealso{
\link{ssi30}, \link{fssi20},
\link{fssa20}, \link{fssa30},
\link{fssi2d}, \link{fssi3d}
}
\concept{square lattice}
\concept{site percolation}
\concept{isotropic cluster}
\concept{von Neumann neighborhood}
\author{Pavel V. Moskalev <moskalefff@gmail.com>}
\examples{
x <- y <- seq(33)
rfq <- fssi30(p=0.37)
image(x, y, rfq[,,17], cex.main=1, 
main="Relative frequency distribution in the z=17 slice\n of isotropic 3D clusters with (1,0)-neighborhood")
contour(x, y, rfq[,,17], levels=c(0.2,0.25,0.3), add=TRUE)
abline(h=17, lty=2); abline(v=17, lty=2)
}
