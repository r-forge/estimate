% $Id$

\name{isprime}
\alias{isprime}
\title{MATLAB isprime function}
\description{
  Array elements that are prime numbers.
}
\usage{
isprime(x)
}
\arguments{
  \item{x}{numeric vector or matrix containing nonnegative integer values}
}
\value{
  Returns an array (or vector) the same size as \code{x} containing logical
  1 (\code{true}) for the elements of \code{x} which are prime, and logical
  0 (\code{false}) otherwise.
}
\seealso{
  \code{\link{factors}},
  \code{\link{primes}}
}
\examples{
x <- c(2, 3, 0, 6, 10)
ans  <- isprime(x)	## 1, 1, 0, 0, 0
as.logical(ans)         ## true, true, false, false, false
}
\author{
  H. Borchers \email{hwborchers@googlemail.com},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\keyword{arith}

