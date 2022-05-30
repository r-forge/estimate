% $Id$

\name{repmat}
\alias{repmat}
\title{MATLAB repmat function}
\description{
  Replicate and tile a matrix.
}
\usage{
repmat(A, \dots)
}
\arguments{
  \item{A}{vector or matrix to be tiled. Must be numeric, logical,
    complex or character.}
  \item{\dots}{numeric dimensions for the result}
}
\value{
  Returns matrix with value \code{A} tiled to the number of dimensions specified.
  Defaults to square if dimension argument resolves to a single value.
}
\seealso{ 
  \code{\link{ones}},
  \code{\link{zeros}}
}   
\examples{
repmat(1, 3)                       # same as ones(3)
repmat(1, c(3, 3))                 # same thing 
repmat(1, 3, 3)                    # same thing
repmat(1, size(matrix(NA, 3, 3)))  # same thing
repmat(matrix(1:4, 2, 2), 3)
}
\author{
  P. Roebuck \email{proebuck@mdanderson.org}
}
\keyword{array}

