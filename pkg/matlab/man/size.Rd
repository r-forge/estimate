% $Id$

\name{size}
\alias{size}
\alias{size,array,integer-method}
\alias{size,array,missing-method}
\alias{size,array,numeric-method}
\alias{size,matrix,integer-method}
\alias{size,matrix,missing-method}
\alias{size,matrix,numeric-method}
\alias{size,vector,missing-method}
\alias{size,vector,numeric-method}
\alias{size,missing,ANY-method}
\title{MATLAB size function}
\description{
  Provides dimensions of \code{X}.
}
\usage{
size(X, dimen)
}
\arguments{
  \item{X}{vector, matrix, or array object}
  \item{dimen}{numeric scalar specifies particular dimension}
}
\details{
  This is an S4 generic function.
  Vector will be treated as a single row matrix. 
  Stored value is equivalent to \code{\link[base]{dim}}.
}
\value{
  Returns object of class \code{size_t} containing the dimensions of
  input argument \code{X} if invoked with a single argument. Returns
  integer value of specified dimension if invoked with two arguments.
  If \code{dimen} specifies a higher dimension than exists, returns 1
  representing the singleton dimension. 
}
\note{
  Handling of vectors is different than in initial implementation.
  Initial implementation returned \code{\link[base]{length}}.
}
\examples{
size(2:9)			# 1 8
size(matrix(1:8, 2, 4))		# 2 4
size(matrix(1:8, 2, 4), 2)	# 4
size(matrix(1:8, 2, 4), 3)	# 1
}
\author{
  P. Roebuck \email{proebuck@mdanderson.org}
}
\keyword{array}

