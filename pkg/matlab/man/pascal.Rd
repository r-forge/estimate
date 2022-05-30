% $Id$

\name{pascal}
\alias{pascal}
\title{MATLAB pascal function}
\description{
  Generate Pascal matrix.
}
\usage{
pascal(n, k=0)
}
\arguments{
  \item{n}{numeric scalar specifying order}
  \item{k}{numeric scalar specifying desired option. Valid values
are 0, 1, or 2}
}
\details{
  Specifying \eqn{k=0} returns symmetric positive definite matrix
  with integer entries taken from Pascal's triangle.

  Specifying \eqn{k=1} returns the lower triangular Cholesky factor
  (up to the signs of the columns) of the Pascal matrix.

  Specifying \eqn{k=2} returns a cube root of the identity matrix.
}
\value{
  Returns matrix of order \code{n} according to specified option \code{k}.
}
\examples{
pascal(4)
pascal(3, 2)
}
\author{
  P. Roebuck \email{proebuck@mdanderson.org}
}
\keyword{array}

