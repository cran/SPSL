\name{SPSL-package}
\alias{SPSL-package}
\alias{SPSL}
\docType{package}
\title{Site Percolation on Square Lattice (SPSL)}
\description{
SPSL package provides functionality for labeling of percolation clusters on 2D & 3D square lattice with various lattice size, relative fraction of accessible sites (occupation probability), iso- & anisotropy, von Neumann & Moore (1,d)-neighborhood.
}
\details{
\tabular{ll}{
Package: \tab SPSL\cr
Type: \tab Package\cr
Version: \tab 0.1-6\cr
Date: \tab 2012-06-07\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}
\code{ssi20()} and \code{ssi30()} functions provide sites labeling of the isotropic cluster on 2D & 3D square lattice with von Neumann (1,0)-neighborhood.
\cr
\code{ssa20()} and \code{ssa30()} functions provide sites labeling of the anisotropic cluster on 2D & 3D square lattice with von Neumann (1,0)-neighborhood.
\cr
\code{ssi2d()} and \code{ssi3d()} functions provide sites labeling of the isotropic cluster on 2D & 3D square lattice with Moore (1,d)-neighborhood.
\cr
\code{ssa2d()} and \code{ssa3d()} functions provide sites labeling of the anisotropic cluster on 2D & 3D square lattice with Moore (1,d)-neighborhood.
\cr
\code{fssi20()} and \code{fssi30()} functions calculates the relative frequency distribution of isotropic clusters on 2D & 3D square lattice with von Neumann (1,0)-neighborhood.
\cr
\code{fssa20()} and \code{fssa30()} functions calculates the relative frequency distribution of anisotropic clusters on 2D & 3D square lattice with von Neumann (1,0)-neighborhood.
\cr
\code{fssi2d()} and \code{fssi3d()} functions calculates the relative frequency distribution of isotropic clusters on 2D & 3D square lattice with Moore (1,d)-neighborhood.
\cr
\code{fssa2d()} and \code{fssa3d()} functions calculates the relative frequency distribution of anisotropic clusters on 2D & 3D square lattice with Moore (1,d)-neighborhood.
}
\author{Pavel V. Moskalev <moskalefff@gmail.com>}
\references{
Moskalev, P.V. and Shitov, V.V. Mathematical modeling of porous structures. Moscow: Fizmatlit, 2007. 120 pp; in Russian.\cr
Moskalev, P.V. (2009), Analysis of the percolation cluster structure.
\emph{Technical Physics}, Vol.54, No.6, pp.763-769.\cr
Moskalev, P.V., Grebennikov, K.V. and Shitov, V.V. (2011),
Statistical estimation of percolation cluster parameters.
\emph{Proceedings of Voronezh State University. Series: Systems Analysis and Information Technologies}, No.1 (January-June), pp.29-35; arXiv:1105.2334v1 [cond-mat.stat-mech]; in Russian.
}