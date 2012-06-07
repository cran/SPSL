\name{fssi20}
\alias{fssi20}
\title{Frequency of Sites on a Square Isotropic 2D lattice with (1,0)-neighborhood}
\description{
\code{fssi20()} function calculates the relative frequency distribution of isotropic clusters on 2D square lattice with von Neumann (1,0)-neighborhood.
}
\usage{
fssi20(n=1000, x=33, p=0.592746, set=(x^2+1)/2, all=TRUE)
}
\arguments{
\item{n}{a sample size.}
\item{x}{a linear dimension of 2D square percolation lattice.}
\item{p}{the relative fractions \code{(0<p)&(p<1)} of accessible sites (occupation probability) for percolation lattice.}
\item{set}{a vector containing the linear indexes of sites from initial subset.}
\item{all}{logical; if \code{all=TRUE}, mark all sites from initial subset; if \code{all=FALSE}, mark accessible sites from initial subset.}
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
Moskalev, P.V., Grebennikov, K.V. and Shitov, V.V. (2011),
Statistical estimation of percolation cluster parameters.
\emph{Proceedings of Voronezh State University. Series: Systems Analysis and Information Technologies}, No.1 (January-June), pp.29-35; arXiv:1105.2334v1 [cond-mat.stat-mech]; in Russian.
}
\examples{
x <- y <- seq(33) 
image(x, y, rfq <- fssi20(), cex.main=1,
main="Relative frequency distribution of\n isotropic 2D clusters with (1,0)-neighborhood")
contour(x, y, rfq, levels=seq(.2,.3,.05), add=TRUE)
abline(h=17, lty=2); abline(v=17, lty=2)
}
