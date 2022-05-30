% $Id$

\name{fliplr}
\alias{fliplr}
\alias{fliplr,ANY-method}
\alias{fliplr,array-method}
\alias{fliplr,matrix-method}
\alias{fliplr,missing-method}
\alias{fliplr,vector-method}
\alias{flipud}
\alias{flipud,ANY-method}
\alias{flipud,array-method}
\alias{flipud,matrix-method}
\alias{flipud,missing-method}
\alias{flipud,vector-method}
\title{MATLAB matrix flip functions}
\description{
  Flips matrices either left-right or up-down.
}
\usage{
fliplr(object)
flipud(object)
}
\arguments{
  \item{object}{vector or matrix to be flipped}
}
\details{
  These are S4 generic functions.
}
\value{
  Return value is the same type as argument \code{object}.
}
\seealso{ 
  \code{\link{rot90}}
}   
\examples{
fliplr(1:9)
flipud(1:9)	# same as previous since vectors have no orientation in R
fliplr(matrix(1:9, 3, 3, byrow=TRUE))
flipud(matrix(1:9, 3, 3, byrow=TRUE))
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\keyword{array}

