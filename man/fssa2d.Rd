\name{fssa2d}
\alias{fssa2d}
\title{Frequency of Sites on a Square Anisotropic 2D lattice with (1,d)-neighborhood}
\description{
\code{fssa2d()} function calculates the relative frequency distribution of anisotropic clusters on 2D square lattice with Moore (1,d)-neighborhood.
}
\usage{
fssa2d(n=1000, x=33, 
       p0=runif(4, max=0.8), 
       p1=colMeans(matrix(p0[c(1,3, 2,3, 1,4, 2,4)], nrow=2))/2,
       set=(x^2+1)/2, all=TRUE, shape=c(1,1))
}
\arguments{
\item{n}{a sample size.}
\item{x}{a linear dimension of 2D square percolation lattice.}
\item{p0}{a vector of relative fractions \code{(0<p0)&(p0<1)} of accessible sites (occupation probability) for lattice directions: \code{(-x,+x,-y,+y)}.}
\item{p1}{averaged double combinations of \code{p0}-components weighted by Minkowski distance: \code{p1=colMeans(matrix(p0[c(1,3,...)], nrow=2))/rhoMe1}.}
\item{set}{a vector of linear indexes of a starting sites subset.}
\item{all}{logical; if \code{all=TRUE}, mark all sites from a starting subset; if \code{all=FALSE}, mark only accessible sites from a starting subset.}
\item{shape}{a vector with two shape parameters of beta-distributed random variables, weighting the percolation lattice sites.}
}
\value{
\item{rfq}{ 
a 2D matrix of relative sampling frequencies for sites of the percolation lattice.}
}
\details{
The percolation is simulated on 2D square lattice with uniformly weighted sites and the vectors \code{p0} and \code{p1}, distributed over the lattice directions, and their combinations.

The anisotropic cluster is formed from the accessible sites connected with the initial subset \code{set}, and depends on the direction in 2D square lattice.

Moore (1,d)-neighborhood on 2D square lattice consists of sites, at least one coordinate of which is different from the current site by one: \code{e=c(e0,e1)}, where \cr
\code{e0=c(-1,} \code{1,} \code{-x,} \code{x)}; \cr
\code{e1=colSums(matrix(e0[c(1,3,} \code{2,3,} \code{1,4,} \code{2,4)], nrow=2))}.

Minkowski distance between sites \code{a} and \code{b} depends on the exponent \code{d}:\cr
\code{rhoM <- function(a, b, d=1)} \cr
\code{if (is.infinite(d)) return(apply(abs(b-a), 2, max))} \cr
\code{else return(apply(abs(b-a)^d, 2, sum)^(1/d))}.

Minkowski distance for sites from \code{e1} subset with the exponent \code{d=1} is equal to \code{rhoMe1=2}.

Each element of the matrix \code{frq} is equal to the relative frequency with which the 2D square lattice site belongs to a cluster sample of size \code{n}.
}
\seealso{
\link{ssa2d}, \link{fssa3d},
\link{fssa20}, \link{fssa30},
\link{fssi2d}, \link{fssi3d}
}
\concept{square lattice}
\concept{site percolation}
\concept{anisotropic cluster}
\concept{Moore neighborhood}
\concept{Minkowski distance}
\author{Pavel V. Moskalev <moskalefff@gmail.com>}
\references{
[1] Moskalev, P.V. Percolation modeling of porous structures. Moscow: URSS, 2018. 240 pp; in Russian.}
\examples{
x <- y <- seq(33) 
image(x, y, rfq <- fssa2d(n=200, p0=c(.3,.4,.75,.5)), cex.main=1,
main="Frequencies of anisotropic (1,1)-clusters")
contour(x, y, rfq, levels=seq(.2,.3,.05), add=TRUE)
abline(h=17, lty=2); abline(v=17, lty=2)
}
