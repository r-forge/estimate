% $Id$

\name{eye}
\alias{eye}
\title{MATLAB eye function}
\description{
  Create an identity matrix.
}
\usage{
eye(m, n)
}
\arguments{
  \item{m, n}{numeric scalar specifying dimensions for the result}
}
\value{
  Returns matrix of order 1.
  Defaults to square if second dimension argument \code{n} not provided.
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\seealso{ 
  \code{\link{ones}},
  \code{\link{zeros}}
}   
\examples{
eye(3)
}
\keyword{array}

