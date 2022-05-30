% $Id$

\name{sum}
\alias{sum}
\alias{sum,array,logical-method}
\alias{sum,array,missing-method}
\alias{sum,matrix,logical-method}
\alias{sum,matrix,missing-method}
\alias{sum,vector,logical-method}
\alias{sum,vector,missing-method}
\alias{sum,ANY,ANY-method}
\alias{sum,missing,ANY-method}
\title{MATLAB sum function}
\description{
  Provides sum of elements.
}
\usage{
sum(x, na.rm=FALSE)
}
\arguments{
  \item{x}{numeric or logical to be summed}
  \item{na.rm}{logical scalar. If \code{TRUE}, remove missing values}
}
\details{
  This is an S4 generic function.
}
\value{
  Return value depends on argument \code{x}. If vector, returns the same as
  \code{\link[base]{sum}}. If matrix, returns vector containing the sum of
  each column.
}
\examples{
sum(1:9)
sum(matrix(1:9, 3, 3))
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\keyword{arith}

