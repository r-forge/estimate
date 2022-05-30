% $Id$

\name{factors}
\alias{factors}
\title{MATLAB factor function}
\description{
  Performs prime factorization.
}
\usage{
factors(n)
}
\arguments{
  \item{n}{numeric scalar specifying composite number to be factored}
}
\details{
  Computes the prime factors of \code{n} in ascending order, each one as often
  as its multiplicity requires, such that \code{n == prod(factors(n))}.
}
\value{
  Returns vector containing the prime factors of \code{n}.
}
\note{
  The corresponding MATLAB function is called 'factor', but was renamed here
  to avoid conflict with R's compound object class.
}
\seealso{
  \code{\link{isprime}},
  \code{\link{primes}}
}
\examples{
factors(1002001)       # 7  7  11  11  13  13
factors(65537)         # is prime
## Euler's calculation
factors(2^32 + 1)      # 641  6700417
}
\author{
  H. Borchers \email{hwborchers@googlemail.com},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\keyword{arith}

