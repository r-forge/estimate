% $Id$

\name{rot90}
\alias{rot90}
\title{MATLAB rot90 function}
\description{
  Rotates matrix counterclockwise \code{k*90} degrees.
}
\usage{
rot90(A, k=1)
}
\arguments{
  \item{A}{matrix to be rotated}
  \item{k}{numeric scalar specifying the number of times to rotate (1..4)}
}
\details{
  Rotating 4 times (360 degrees) returns the original matrix unchanged.
}
\value{
  Returns matrix corresponding to argument \code{A} having been rotated
  argument \code{k} number of times.
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\seealso{ 
  \code{\link{fliplr}},
  \code{\link{flipud}}
}   
\examples{
rot90(matrix(1:4, 2, 2))
}
\keyword{array}

