% $Id$

\name{mod}
\alias{mod}
\alias{rem}
\title{MATLAB mod/rem functions}
\description{
  Provides modulus and remainder after division.
}
\usage{
mod(x, y)
rem(x, y)
}
\arguments{
  \item{x, y}{numeric vectors or objects}
}
\value{
  Returns vector containing result of the element by element operations.
}
\examples{
## same results with x, y having the same sign
mod(5, 3)
rem(5, 3)
## same results with x, y having different signs
mod(-5, 3)
rem(-5, 3)
}
\author{
  P. Roebuck \email{proebuck@mdanderson.org}
}
\keyword{arith}

