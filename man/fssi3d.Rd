\name{fssi3d}
\alias{fssi3d}
\title{Frequency of Sites on a Square Isotropic 3D lattice with (1,d)-neighborhood}
\description{
\code{fssi3d()} function calculates the relative frequency distribution of isotropic clusters on 3D square lattice with Moore (1,d)-neighborhood.
}
\usage{
fssi3d(n=1000, x=33, p0=0.2, p1=p0/2, p2=p0/3, 
       set=(x^3+1)/2, all=TRUE, shape=c(1,1))
}
\arguments{
\item{n}{a sample size.}
\item{x}{a linear dimension of 3D square percolation lattice.}
\item{p0}{a relative fraction \code{(0<p0)&(p0<1)} of accessible sites (occupation probability) for percolation lattice.}
\item{p1}{\code{p0} value, weighted by Minkowski distance: \code{p1=p0/rhoMe1}.}
\item{p2}{\code{p0} value, weighted by Minkowski distance: \code{p2=p0/rhoMe2}.}
\item{set}{a vector of linear indexes of a starting sites subset.}
\item{all}{logical; if \code{all=TRUE}, mark all sites from a starting subset; if \code{all=FALSE}, mark only accessible sites from a starting subset.}
\item{shape}{a vector with two shape parameters of beta-distributed random variables, weighting the percolation lattice sites.}
}
\value{
\item{rfq}{ 
a 3D matrix of relative sampling frequencies for sites of the percolation lattice.}
}
\details{
The percolation is simulated on 3D square lattice with uniformly weighted sites and the constant parameters  \code{p0}, \code{p1}, and \code{p2}.

The isotropic cluster is formed from the accessible sites connected with initial sites subset \code{set}. 

Moore (1,d)-neighborhood on 3D square lattice consists of sites, at least one coordinate of which is different from the current site by one: \code{e=c(e0,e1,e2)}, where \cr
\code{e0=c(-1,} \code{1,} \code{-x,} \code{x,} \code{-x^2,} \code{x^2)}; \cr
\code{e1=colSums(matrix(e0[c(1,3,} \code{2,3,} \code{1,4,} \code{2,4,} \code{1,5,} \code{2,5,} \code{1,6,} \code{2,6,} \code{3,5,} \code{4,5,} \code{3,6,} \code{4,6)], nrow=2))}; \cr
\code{e2=colMeans(matrix(p0[c(1,3,5,} \code{2,3,5,} \code{1,4,5,} \code{2,4,5,} \code{1,3,6,} \code{2,3,6,} \code{1,4,6,} \code{2,4,6)], nrow=3))}.

Minkowski distance between sites \code{a} and \code{b} depends on the exponent \code{d}:\cr
\code{rhoM <- function(a, b, d=1)} \cr
\code{if (is.infinite(d)) return(apply(abs(b-a), 2, max))} \cr
\code{else return(apply(abs(b-a)^d, 2, sum)^(1/d))}.

Minkowski distance for sites from \code{e1} and \code{e2} subsets with the exponent \code{d=1} is equal to \code{rhoMe1=2} and \code{rhoMe2=3}.

Each element of the matrix \code{frq} is equal to the relative frequency with which the 3D square lattice site belongs to a cluster sample of size \code{n}.
}
\seealso{
\link{ssi3d}, \link{fssi2d},
\link{fssi20}, \link{fssi30},
\link{fssa2d}, \link{fssa3d}
}
\concept{square lattice}
\concept{site percolation}
\concept{isotropic cluster}
\concept{Moore neighborhood}
\concept{Minkowski distance}
\author{Pavel V. Moskalev <moskalefff@gmail.com>}
\references{
[1] Moskalev, P.V. Percolation modeling of porous structures. Moscow: URSS, 2018. 240 pp; in Russian.}
\examples{
x <- y <- seq(33)
rfq <- fssi3d(n=200, p0=.285)
image(x, y, rfq[,,17], cex.main=1, 
main="Frequencies in the z=17 slice of isotropic (1,1)-clusters")
contour(x, y, rfq[,,17], levels=c(0.2,0.25,0.3), add=TRUE)
abline(h=17, lty=2); abline(v=17, lty=2)
}