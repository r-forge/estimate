% $Id$

\name{vander}
\alias{vander}
\title{MATLAB vander function}
\description{
  Generate Vandermonde matrix from a vector of numbers.
}
\usage{
vander(v)
}
\arguments{
  \item{v}{numeric or complex vector of values}
}
\details{
  Generates the Vandermonde matrix whose columns are powers of the vector
  \code{v} (of length \code{n}) using the formula \preformatted{
    A[i, j] = v[i]^(n-j)
}
  Used when fitting a polynomial to given points.
}
\value{
  Returns an \code{n}-by-\code{n} matrix constructed as described above.
}
\examples{
vander(1:5)
}
\author{
  H. Borchers \email{hwborchers@googlemail.com},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\keyword{array}

