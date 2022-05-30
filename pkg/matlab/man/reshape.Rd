% $Id$

\name{reshape}
\alias{reshape}
\title{MATLAB reshape function}
\description{
  Reshape matrix or array.
}
\usage{
reshape(A, \dots)
}
\arguments{
  \item{A}{matrix or array containing the original data}
  \item{\dots}{numeric dimensions for the result}
}
\details{
  In the first example below, an \code{m}-by-\code{n} matrix is created whose
  elements are taken column-wise from \code{A}. An error occurs if \code{A}
  does not have \eqn{m*n} elements.\cr

  In the second example below, an \code{n}-dimensional array with the same
  elements as \code{A} but reshaped to have the size
  \code{m}-by-\code{n}-by-\code{p}. The product of the specified dimensions
  must be the same as \code{prod(size(A))}.\cr

  In the third example below, an \code{n}-dimensional array with the same
  elements as \code{A} but reshaped to \code{siz}, a vector representing the
  dimensions of the reshaped array. The quantity \code{prod(siz)} must be
  the same as \code{prod(size(A))}.
}   
\value{
  Returns matrix (or array) of requested dimensions containing the elements
  of \code{A}.
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\examples{
Xmat.2d <- matrix(1:12, nrow=4, ncol=3)
reshape(Xmat.2d, 6, 2)                   # example 1
reshape(Xmat.2d, c(6, 2))                    # same thing
Xarr.3d <- reshape(Xmat.2d, c(6, 2, 1))  # example 2
reshape(Xmat.2d, size(Xarr.3d))          # example 3
}
\keyword{array}

