\name{ssi20}
\alias{ssi20}
\title{Site cluster on Square Isotropic 2D lattice with (1,0)-neighborhood}
\description{
\code{ssi20()} function provides sites labeling of the isotropic cluster on 2D square lattice with von Neumann (1,0)-neighborhood.
}
\usage{
ssi20(x=33, p=0.592746, 
      set=(x^2+1)/2, all=TRUE, shape=c(1,1))
}
\arguments{
\item{x}{a linear dimension of 2D square percolation lattice.}
\item{p}{the relative fractions \code{(0<p)&(p<1)} of accessible sites (occupation probability) for percolation lattice.}
\item{set}{a vector of linear indexes of a starting sites subset.}
\item{all}{logical; if \code{all=TRUE}, mark all sites from a starting subset; if \code{all=FALSE}, mark only accessible sites from a starting subset.}
\item{shape}{a vector with two shape parameters of beta-distributed random variables, weighting the percolation lattice sites.}
}
\value{
\item{acc}{ 
an accessiblity matrix for 2D square percolation lattice:\cr 
if \code{acc[e]<p} then \code{acc[e]} is accessible site;\cr 
if \code{acc[e]==1} then \code{acc[e]} is non-accessible site;\cr 
if \code{acc[e]==2} then \code{acc[e]} belongs to a sites cluster. }
}
\details{
The percolation is simulated on 2D square lattice with uniformly weighted sites \code{acc} and the constant parameter \code{p}.

The isotropic cluster is formed from the accessible sites connected with initial sites subset. 

To form the cluster the condition \code{acc[set+e]<p} is iteratively tested for sites of the von Neumann (1,0)-neighborhood \code{e} for the current cluster perimeter \code{set}.

Von Neumann (1,0)-neighborhood on 2D square lattice consists of sites, only one coordinate of which is different from the current site by one: \code{e=c(-1,} \code{1,} \code{-x,} \code{x)}.

Forming cluster ends with the exhaustion of accessible sites in von Neumann (1,0)-neighborhood of the current cluster perimeter.
}
\seealso{
\link{fssi20}, \link{ssi30},
\link{ssa20}, \link{ssa30},
\link{ssi2d}, \link{ssi3d}
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
set.seed(20120507)
x <- y <- seq(33)
image(x, y, ssi20(), zlim=c(0,2), 
main="Isotropic (1,0)-cluster")
abline(h=17, lty=2); abline(v=17, lty=2)
}
