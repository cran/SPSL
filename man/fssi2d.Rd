\name{fssi2d}
\alias{fssi2d}
\title{Frequency of Sites on a Square Isotropic 2D lattice with (1,d)-neighborhood}
\description{
\code{fssi2d()} function calculates the relative frequency distribution of isotropic clusters on 2D square lattice with Moore (1,d)-neighborhood.
}
\usage{
fssi2d(n=1000, x=33, p0=0.5, p1=p0/2, set=(x^2+1)/2, all=TRUE)
}
\arguments{
\item{n}{a sample size.}
\item{x}{a linear dimension of 2D square percolation lattice.}
\item{p0}{a relative fraction \code{(0<p0)&(p0<1)} of accessible sites (occupation probability) for percolation lattice.}
\item{p1}{\code{p0} value, weighted by non-metrical Minkowski distance: \code{p1=p0/rhoMe1}.}
\item{set}{a vector of linear indexes of initial sites subset.}
\item{all}{logical; if \code{all=TRUE}, mark all sites from initial subset; if \code{all=FALSE}, mark accessible sites from initial subset.}
}
\value{
\item{rfq}{ 
a 2D matrix of relative sampling frequencies for sites of the percolation lattice.}
}
\details{
The percolation is simulated on 2D square lattice with uniformly weighted sites and the constant parameters \code{p0} and \code{p1}.

The isotropic cluster is formed from the accessible sites connected with initial sites subset \code{set}. 

Moore (1,d)-neighborhood on 2D square lattice consists of sites, at least one coordinate of which is different from the current site by one: \code{e=c(e0,e1)}, where 
\code{e0=c(-1,} \code{1,} \code{-x,} \code{x,} \code{-x^2,} \code{x^2)}; 
\code{e1=colSums(matrix(e0[c(1,3,} \code{2,3,} \code{1,4,} \code{2,4)], nrow=2))}.

Minkowski non-metrical distance between sites \code{a} and \code{b} depends on the exponent \code{d}:\cr
\code{rhoM <- function(a, b, d=1)} \cr
\code{if (is.infinite(d)) return(apply(abs(b-a), 2, max))} \cr
\code{else return(apply(abs(b-a)^d, 2, sum)^(1/d))}.

Minkowski non-metrical distance for sites from \code{e1} subset with the exponent \code{d=1} is equal to \code{rhoMe1=2}.

Each element of the matrix \code{frq} is equal to the relative frequency with which the 2D square lattice site belongs to a cluster sample of size \code{n}.
}
\seealso{
\link{ssi2d}, \link{fssi3d},
\link{fssi20}, \link{fssi30},
\link{fssa2d}, \link{fssa3d}
}
\concept{square lattice}
\concept{site percolation}
\concept{isotropic cluster}
\concept{Moore neighborhood}
\concept{non-metrical Minkowski distance}
\author{Pavel V. Moskalev <moskalefff@gmail.com>}
\examples{
x <- y <- seq(33)
image(x, y, rfq <- fssi2d(), cex.main=1,
main="Relative frequency distribution of\n isotropic 2D clusters with (1,1)-neighborhood")
contour(x, y, rfq, levels=seq(.2,.3,.05), add=TRUE)
abline(h=17, lty=2); abline(v=17, lty=2)
}
