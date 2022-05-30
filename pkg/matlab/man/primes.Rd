% $Id$

\name{primes}
\alias{primes}
\title{MATLAB primes function}
\description{
  Generate a list of prime numbers.
}
\usage{
primes(n)
}
\arguments{
  \item{n}{scalar numeric specifying largest prime number desired.}
}
\details{
  Generates the list of prime numbers less than or equal to \code{n} using a
  variant of the basic "Sieve of Eratosthenes" algorithm. This approach is
  reasonably fast, but requires a copious amount of memory when \code{n} is
  large. A prime number is one that has no other factors other than \code{1}
  and itself.
}
\value{
  Returns numeric vector containing prime numbers less than or equal to
  argument \code{n}.
}
\seealso{
  \code{\link{isprime}},
  \code{\link{factors}}
}
\examples{
primes(1000)
length(primes(1e6))    #   78498 prime numbers less than one million
\dontrun{
  length(primes(1e7))  #  664579 prime numbers less than ten million
  length(primes(1e8))  # 5761455 prime numbers less than one hundred million
}
}
\author{
  H. Borchers \email{hwborchers@googlemail.com},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\keyword{arith}

