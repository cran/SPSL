\name{SPSL-package}
\alias{SPSL-package}
\alias{SPSL}
\docType{package}
\title{Site Percolation on Square Lattices (SPSL)}
\description{
Provides basic functionality for labeling iso- & anisotropic percolation clusters on 2D & 3D square lattices with various lattice sizes, occupation probabilities, von Neumann & Moore (1,d)-neighborhoods, and random variables weighting the percolation lattice sites.
}
\details{
\tabular{ll}{
Package: \tab SPSL\cr
Type: \tab Package\cr
Version: \tab 0.1-9\cr
Date: \tab 2019-03-17\cr
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
[1] Moskalev, P.V. Percolation modeling of porous structures. Moscow: URSS, 2018. 240 pp; in Russian.\cr
[2] Moskalev, P.V. (2014) Estimates of threshold and strength of percolation clusters on square lattices with (1,d)-neighborhood. \emph{Computer Research and Modeling}, Vol.6, No.3, pp.405–414; in Russian.\cr
[3] Moskalev, P.V. (2013) The structure of site percolation models on three-dimensional square lattices. \emph{Computer Research and Modeling}, Vol.5, No.4, pp.607–622; in Russian.}