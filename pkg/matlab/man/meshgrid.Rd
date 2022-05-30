% $Id$

\name{meshgrid}
\alias{meshgrid}
\title{MATLAB meshgrid functions}
\description{
  Generate X and Y matrices for three-dimensional plots.
}
\usage{
meshgrid(x, y, z, nargout=2)
}
\arguments{
  \item{x, y, z}{numeric vectors of values}
  \item{nargout}{numeric scalar that determines number of dimensions to return}
}
\details{
  In the first example below, the domain specified by vectors \code{x} and \code{y}
  are transformed into two arrays which can be used to evaluate functions of two
  variables and three-dimensional surface plots. The rows of the output array
  \code{x} are copies of the vector \code{x}; columns of the output array \code{y}
  are copies of the vector \code{y}.\cr

  The second example below is syntactic sugar for specifying
  \code{meshgrid(x, x)}.\cr

  The third example below produces three-dimensional arrays used to evaluate
  functions of three variables and three-dimensional volumetric plots. 
}
\value{
  Returns list containing eiher two or three matrices depending on the
  value of \code{nargout}.
    \item{x, y, z}{output matrices}
}
\note{
  Limited to two- or three-dimensional Cartesian space.
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\examples{
meshgrid(1:3, 10:14)            # example 1
meshgrid(1:3)                   # example 2
meshgrid(5:8, 10:14, 2:3, 3)    # example 3
}
\keyword{array}

