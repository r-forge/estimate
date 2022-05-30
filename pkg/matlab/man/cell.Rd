% $Id$

\name{cell}
\alias{cell}
\title{MATLAB cell function}
\description{
  Create cell array.
}
\usage{
cell(\dots)
}
\arguments{
  \item{\dots}{numeric dimensions for the result}
}
\value{
  Returns list consisting of empty matrices.
  Defaults to square if dimension argument resolves to a single value.
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\seealso{ 
  \code{\link{ones}},
  \code{\link{zeros}}
}   
\examples{
cell(3)
cell(c(3, 3))			# same thing
cell(3, 3)			# same thing
cell(size(matrix(NA, 3, 3)))	# same thing
}
\keyword{array}
\keyword{list}

