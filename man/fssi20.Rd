\name{fssi20}
\alias{fssi20}
\title{Frequency of Sites on a Square Isotropic 2D lattice with (1,0)-neighborhood}
\description{
\code{fssi20()} function calculates the relative frequency distribution of isotropic clusters on 2D square lattice with von Neumann (1,0)-neighborhood.
}
\usage{
fssi20(n=1000, x=33, p=0.592746, 
       set=(x^2+1)/2, all=TRUE, shape=c(1,1))
}
\arguments{
\item{n}{a sample size.}
\item{x}{a linear dimension of 2D square percolation lattice.}
\item{p}{the relative fractions \code{(0<p)&(p<1)} of accessible sites (occupation probability) for percolation lattice.}
\item{set}{a vector of linear indexes of a starting sites subset.}
\item{all}{logical; if \code{all=TRUE}, mark all sites from a starting subset; if \code{all=FALSE}, mark only accessible sites from a starting subset.}
\item{shape}{a vector with two shape parameters of beta-distributed random variables, weighting the percolation lattice sites.}
}
\value{
\item{rfq}{ 
a 2D matrix of relative sampling frequencies for sites of the percolation lattice.}
}
\details{
The percolation is simulated on 2D square lattice with uniformly weighted sites and the constant parameter \code{p}.

The isotropic cluster is formed from the accessible sites connected with initial sites subset \code{set}.

Von Neumann (1,0)-neighborhood on 2D square lattice consists of sites, only one coordinate of which is different from the current site by one: \code{e=c(-1,} \code{1,} \code{-x,} \code{x)}.

Each element of the matrix \code{frq} is equal to the relative frequency with which the 2D square lattice site belongs to a cluster sample of size \code{n}.
}
\seealso{
\link{ssi20}, \link{fssi30},
\link{fssa20}, \link{fssa30},
\link{fssi2d}, \link{fssi3d}
}
\concept{square lattice}
\concept{site percolation}
\concept{isotropic cluster}
\concept{von Neumann neighborhood}
\author{Pavel V. Moskalev}
\references{
[1] Moskalev, P.V. Percolation modeling of porous structures. Moscow: URSS, 2018. 240 pp; in Russian.\cr
[2] Moskalev, P.V. (2014) Estimates of threshold and strength of percolation clusters on square lattices with (1,d)-neighborhood. \emph{Computer Research and Modeling}, Vol.6, No.3, pp.405–414; in Russian.}
\examples{
x <- y <- seq(33) 
image(x, y, rfq <- fssi20(n=200), cex.main=1,
main="Frequencies of isotropic (1,0)-clusters")
contour(x, y, rfq, levels=seq(.2,.3,.05), add=TRUE)
abline(h=17, lty=2); abline(v=17, lty=2)
}
