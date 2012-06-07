\name{fssa3d}
\alias{fssa3d}
\title{Frequency of Sites on a Square Anisotropic 3D lattice with (1,d)-neighborhood}
\description{
\code{fssa3d()} function calculates the relative frequency distribution of anisotropic clusters on 3D square lattice with Moore (1,d)-neighborhood.
}
\usage{
fssa3d(n=1000, x=33, 
       p0=runif(6, max=0.4),
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
\item{n}{a sample size.}
\item{x}{a linear dimension of 2D square percolation lattice.}
\item{p0}{a vector of relative fractions \code{(0<p0)&(p0<1)} of accessible sites (occupation probability) for lattice directions: \code{(-x,+x,-y,+y,-z,+z)}.}
\item{p1}{averaged double combinations of \code{p0}-components weighted by non-metrical Minkowski distance: \code{p1=colMeans(matrix(p0[c(1,3,...)], nrow=2))/rhoMe1}.}
\item{p2}{averaged triple combinations of \code{p0}-components weighted by non-metrical Minkowski distance: \code{p2=colMeans(matrix(p0[c(1,3,5,...)], nrow=3))/rhoMe2}.}
\item{set}{a vector of linear indexes of initial sites subset.}
\item{all}{logical; if \code{all=TRUE}, mark all sites from initial subset; if \code{all=FALSE}, mark accessible sites from initial subset.}
}
\value{
\item{rfq}{ 
a 3D matrix of relative sampling frequencies for sites of the percolation lattice.}
}
\details{
The percolation is simulated on 3D square lattice with uniformly weighted sites \code{acc} and the vectors \code{p0}, \code{p1}, and \code{p2}, distributed over the lattice directions, and their combinations.

The anisotropic cluster is formed from the accessible sites connected with the initial subset \code{set}, and depends on the direction in 3D square lattice.

Moore (1,d)-neighborhood on 3D square lattice consists of sites, at least one coordinate of which is different from the current site by one: \code{e=c(e0,e1,e2)}, where 
\code{e0=c(-1,} \code{1,} \code{-x,} \code{x,} \code{-x^2,} \code{x^2)}; 
\code{e1=colSums(matrix(e0[c(1,3,} \code{2,3,} \code{1,4,} \code{2,4,} \code{1,5,} \code{2,5,} \code{1,6,} \code{2,6,} \code{3,5,} \code{4,5,} \code{3,6,} \code{4,6)], nrow=2))};
\code{e2=colMeans(matrix(p0[c(1,3,5,} \code{2,3,5,} \code{1,4,5,} \code{2,4,5,} \code{1,3,6,} \code{2,3,6,} \code{1,4,6,} \code{2,4,6)], nrow=3))}.

Minkowski non-metrical distance between sites \code{a} and \code{b} depends on the exponent \code{d}:\cr
\code{rhoM <- function(a, b, d=1)} \cr
\code{if (is.infinite(d)) return(apply(abs(b-a), 2, max))} \cr
\code{else return(apply(abs(b-a)^d, 2, sum)^(1/d))}.

Minkowski non-metrical distance for sites from \code{e1} and \code{e2} subsets with the exponent \code{d=1} is equal to \code{rhoMe1=2} and \code{rhoMe2=3}.

Each element of the matrix \code{frq} is equal to the relative frequency with which the 3D square lattice site belongs to a cluster sample of size \code{n}.
}
\seealso{
\link{ssa3d}, \link{fssa2d},
\link{fssa20}, \link{fssa30},
\link{fssi2d}, \link{fssi3d}
}
\concept{square lattice}
\concept{site percolation}
\concept{anisotropic cluster}
\concept{Moore neighborhood}
\concept{non-metrical Minkowski distance}
\author{Pavel V. Moskalev <moskalefff@gmail.com>}
\examples{
x <- y <- seq(33)
rfq <- fssa3d(p0=.17*c(.5,3,.5,1.5,1,.5))
image(x, y, rfq[,,17], cex.main=1,
main="Relative frequency distribution in the z=17 slice\n of anisotropic 3D clusters with (1,1)-neighborhood")
contour(x, y, rfq[,,17], levels=seq(.05,.3,.05), add=TRUE)
abline(h=17, lty=2); abline(v=17, lty=2)
}
